import numpy as np
import pandas as pd
import torch

from src.data.stock_data_downloader import StockDataDownloader
from src.features.feature_engineer import FeatureEngineer
from src.utils.model_loader import ModelLoader


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class ForecastPredictor:
    def __init__(self, model_path):
        self.model_path = model_path
        self.package, self.model = ModelLoader.load_package(model_path)

    @staticmethod
    def enable_dropout(model):
        for module in model.modules():
            if isinstance(module, torch.nn.Dropout):
                module.train()

    def prepare_latest_sequence(self):
        ticker = self.package["ticker"]
        seq_length = self.package["seq_length"]
        feature_names = self.package["feature_names"]
        scaler = self.package["scaler"]
        period = self.package.get("prediction_download_period", "7y")

        df = StockDataDownloader.download_stock_data(ticker=ticker, period=period)
        if df is None:
            raise RuntimeError(f"Could not download latest data for {ticker}")

        df = StockDataDownloader.calculate_technical_indicators(df)
        X = FeatureEngineer.select_features_by_names(df, feature_names)

        if len(X) < seq_length:
            raise RuntimeError(f"Need {seq_length} rows after indicators, got {len(X)}")

        latest_sequence = X[-seq_length:]
        latest_sequence_scaled = scaler.transform(latest_sequence)
        current_price = float(df["Close"].iloc[-1])
        current_date = df.index[-1]

        return latest_sequence_scaled, current_price, current_date

    def predict_latest_path(self, n_simulations=None):
        if n_simulations is None:
            n_simulations = self.package.get("mc_simulations", 100)

        latest_sequence, current_price, current_date = self.prepare_latest_sequence()
        return self.predict_path_with_confidence(latest_sequence, current_price, current_date, n_simulations)

    def predict_path_with_confidence(self, latest_sequence, current_price, current_date, n_simulations=100):
        horizons = self.package["horizons"]
        X = torch.FloatTensor(latest_sequence).unsqueeze(0).to(device)

        predictions = []
        self.model.eval()
        ForecastPredictor.enable_dropout(self.model)

        with torch.no_grad():
            for _ in range(n_simulations):
                pred = self.model(X).cpu().numpy()[0]
                predictions.append(pred)

        predictions = np.array(predictions)
        mean_log = predictions.mean(axis=0)
        std_log = predictions.std(axis=0)
        lower_log = mean_log - 1.96 * std_log
        upper_log = mean_log + 1.96 * std_log

        names = list(horizons.keys())
        days = list(horizons.values())
        forecast_dates = [pd.Timestamp(current_date) + pd.tseries.offsets.BDay(day) for day in days]

        df = pd.DataFrame({
            "horizon_name": names,
            "trading_day_ahead": days,
            "forecast_date": forecast_dates,
            "expected_price": current_price * np.exp(mean_log),
            "lower_price": current_price * np.exp(lower_log),
            "upper_price": current_price * np.exp(upper_log),
            "expected_return_pct": (np.exp(mean_log) - 1) * 100,
            "std_log_return": std_log
        })

        return {
            "ticker": self.package["ticker"],
            "model_type": self.package["model_type"],
            "current_date": str(pd.Timestamp(current_date).date()),
            "current_price": float(current_price),
            "forecast_path": df
        }

    @staticmethod
    def get_milestones_from_forecast_df(forecast_df):
        milestones = {
            "1d": 1,
            "1m": 21,
            "1y": 252,
            "5y": 1260
        }
        results = {}
        for label, day in milestones.items():
            row = forecast_df[forecast_df["trading_day_ahead"] == day]
            if len(row) == 0:
                continue
            row = row.iloc[0]
            results[label] = {
                "forecast_date": str(pd.Timestamp(row["forecast_date"]).date()),
                "expected_price": float(row["expected_price"]),
                "lower_price": float(row["lower_price"]),
                "upper_price": float(row["upper_price"]),
                "expected_return_pct": float(row["expected_return_pct"]),
                "std_log_return": float(row["std_log_return"])
            }
        return results
