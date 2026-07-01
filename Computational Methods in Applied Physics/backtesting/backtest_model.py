import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from src.data.stock_data_downloader import StockDataDownloader
from src.features.feature_engineer import FeatureEngineer
from src.utils.model_loader import ModelLoader
from src.plotting import ForecastPlotter
from training.train_models import predict_daily_path_from_origin_date


def backtest_saved_model(model_path):
    package, model = ModelLoader.load_package(model_path)

    ticker = package["ticker"]
    horizons = package["horizons"]
    seq_length = package["seq_length"]
    feature_names = package["feature_names"]
    scaler = package["scaler"]

    output_dir = Path("outputs") / ticker
    plot_dir = output_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    df_raw = StockDataDownloader.download_stock_data(ticker=ticker, start=package["download_start"], end=package["download_end"])
    if df_raw is None:
        raise RuntimeError(f"No data downloaded for {ticker}")

    df_features = StockDataDownloader.calculate_technical_indicators(df_raw)

    forecast_df, origin_price, origin_date = predict_daily_path_from_origin_date(
        model=model,
        df_features=df_features,
        feature_names=feature_names,
        scaler=scaler,
        seq_length=seq_length,
        origin_date=package["train_end"],
        horizons=horizons,
        n_simulations=package.get("mc_simulations", 50)
    )

    csv_path = output_dir / "daily_forecast_path_from_training_end_from_saved_model.csv"
    forecast_df.to_csv(csv_path, index=False)

    plot_path = ForecastPlotter.plot_daily_forecast_curve_from_origin(
        full_dates=df_features.index,
        full_close_prices=df_features["Close"].values,
        forecast_df=forecast_df,
        origin_date=origin_date,
        origin_price=origin_price,
        ticker=ticker,
        train_start=package["train_start"],
        train_end=package["train_end"],
        save_dir=str(plot_dir)
    )

    print("Backtest plot regenerated.")
    print(f"CSV: {csv_path}")
    print(f"Plot: {plot_path}")


def parse_args():
    parser = argparse.ArgumentParser(description="Regenerate daily-path plot from a saved model")
    parser.add_argument("--model-path", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    backtest_saved_model(args.model_path)


if __name__ == "__main__":
    main()
