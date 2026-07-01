import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from torch.utils.data import DataLoader

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from src.data.stock_data_downloader import StockDataDownloader
from src.features.feature_engineer import FeatureEngineer
from src.models.neural_models import CNNRegressor, TDNNRegressor, CNNLSTMRegressor
from src.plotting import ForecastPlotter
from src.predictor import ForecastPredictor

from training.dataset import StockDataset
from training.trainer import Trainer
from training.evaluator import Evaluator
from training.baselines import BaselineModels


def build_daily_horizons(max_days=1512):
    return {f"d{day:04d}": day for day in range(1, max_days + 1)}


def build_daily_horizon_weights(horizons):
    days = np.array(list(horizons.values()), dtype=float)
    weights = 1.0 / np.sqrt(days)
    weights = weights / weights.max()
    return weights.tolist()


def predict_daily_path_from_origin_date(model, df_features, feature_names, scaler, seq_length, origin_date, horizons, n_simulations=100):
    import torch
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    origin_date = pd.Timestamp(origin_date)
    df_until_origin = df_features[df_features.index <= origin_date].copy()

    if len(df_until_origin) < seq_length:
        raise RuntimeError(f"Not enough data before origin date {origin_date.date()}. Need {seq_length}, got {len(df_until_origin)}")

    X_all = FeatureEngineer.select_features_by_names(df_until_origin, feature_names)
    latest_sequence = X_all[-seq_length:]
    latest_sequence_norm = scaler.transform(latest_sequence)

    origin_price = float(df_until_origin["Close"].iloc[-1])
    actual_origin_date = df_until_origin.index[-1]

    X = torch.FloatTensor(latest_sequence_norm).unsqueeze(0).to(device)
    predictions = []

    model = model.to(device)
    model.eval()
    ForecastPredictor.enable_dropout(model)

    with torch.no_grad():
        for _ in range(n_simulations):
            pred = model(X).cpu().numpy()[0]
            predictions.append(pred)

    predictions = np.array(predictions)
    mean_log = predictions.mean(axis=0)
    std_log = predictions.std(axis=0)
    lower_log = mean_log - 1.96 * std_log
    upper_log = mean_log + 1.96 * std_log

    names = list(horizons.keys())
    days = list(horizons.values())
    forecast_dates = [pd.Timestamp(actual_origin_date) + pd.tseries.offsets.BDay(day) for day in days]

    forecast_df = pd.DataFrame({
        "horizon_name": names,
        "trading_day_ahead": days,
        "forecast_date": forecast_dates,
        "expected_price": origin_price * np.exp(mean_log),
        "lower_price": origin_price * np.exp(lower_log),
        "upper_price": origin_price * np.exp(upper_log),
        "expected_return_pct": (np.exp(mean_log) - 1) * 100,
        "std_log_return": std_log
    })

    return forecast_df, origin_price, actual_origin_date


def build_milestone_backtest_csv(ticker, dates_test, prices_test, y_true, y_pred, horizon_names, horizons):
    milestone_days = {"1d": 1, "1m": 21, "1y": 252, "5y": 1260, "6y": 1512}
    rows = []
    dates_test = pd.to_datetime(dates_test)

    horizon_to_index = {h: i for i, h in enumerate(horizon_names)}
    day_to_horizon = {v: k for k, v in horizons.items()}

    for label, day in milestone_days.items():
        if day not in day_to_horizon:
            continue
        h = day_to_horizon[day]
        i = horizon_to_index[h]
        true_log = y_true[:, i]
        pred_log = y_pred[:, i]
        actual_price = prices_test * np.exp(true_log)
        predicted_price = prices_test * np.exp(pred_log)
        forecast_dates = [d + pd.tseries.offsets.BDay(day) for d in dates_test]

        for j in range(len(dates_test)):
            rows.append({
                "ticker": ticker,
                "milestone": label,
                "trading_day_ahead": day,
                "prediction_date": str(dates_test[j].date()),
                "forecast_date": str(pd.Timestamp(forecast_dates[j]).date()),
                "current_price": float(prices_test[j]),
                "actual_future_price": float(actual_price[j]),
                "predicted_future_price": float(predicted_price[j]),
                "price_error": float(predicted_price[j] - actual_price[j])
            })

    return pd.DataFrame(rows)


def train_one_ticker(args, ticker):
    print("\n" + "=" * 100)
    print(f"TRAINING: {ticker}")
    print("=" * 100)

    horizons = build_daily_horizons(args.max_horizon_days)
    horizon_names = list(horizons.keys())
    horizon_weights = build_daily_horizon_weights(horizons)

    output_dir = Path("outputs") / ticker
    plot_dir = output_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    df_raw = StockDataDownloader.download_stock_data(ticker=ticker, start=args.download_start, end=args.download_end)
    if df_raw is None:
        raise RuntimeError(f"No data downloaded for {ticker}")

    df_features = StockDataDownloader.calculate_technical_indicators(df_raw)
    df_targets = StockDataDownloader.create_forecast_targets(df_features, horizons)

    print(f"Usable supervised data: {df_targets.index[0].date()} to {df_targets.index[-1].date()}")
    print(f"Target columns: {len(horizon_names)} daily horizons")

    X, feature_names = FeatureEngineer.select_features(df_targets)
    y, _ = FeatureEngineer.extract_targets(df_targets, horizon_names)
    dates = df_targets.index.to_numpy()
    close_prices = df_targets["Close"].values

    X_seq, y_seq, date_seq = FeatureEngineer.create_sequences(X, y, dates, seq_length=args.seq_length, stride=1)
    current_prices_seq = close_prices[args.seq_length - 1: args.seq_length - 1 + len(X_seq)]

    (X_train, X_val, X_test, y_train, y_val, y_test, dates_train, dates_val, dates_test, prices_train, prices_val, prices_test) = FeatureEngineer.split_by_date_masks(
        X_seq, y_seq, date_seq, current_prices_seq,
        train_start=args.train_start,
        train_end=args.train_end,
        val_end=args.val_end,
        test_end=args.test_end
    )

    if len(X_train) == 0 or len(X_val) == 0 or len(X_test) == 0:
        raise RuntimeError("Train/validation/test split produced an empty set. Check dates and max horizon.")

    print(f"Train: {X_train.shape}, {dates_train[0]} to {dates_train[-1]}")
    print(f"Val:   {X_val.shape}, {dates_val[0]} to {dates_val[-1]}")
    print(f"Test:  {X_test.shape}, {dates_test[0]} to {dates_test[-1]}")

    X_train_norm, X_val_norm, X_test_norm, scaler = FeatureEngineer.normalize_features(
        X_train.reshape(-1, X_train.shape[-1]),
        X_val.reshape(-1, X_val.shape[-1]),
        X_test.reshape(-1, X_test.shape[-1])
    )
    X_train_norm = X_train_norm.reshape(X_train.shape)
    X_val_norm = X_val_norm.reshape(X_val.shape)
    X_test_norm = X_test_norm.reshape(X_test.shape)

    train_loader = DataLoader(StockDataset(X_train_norm, y_train), batch_size=args.batch_size, shuffle=False)
    val_loader = DataLoader(StockDataset(X_val_norm, y_val), batch_size=args.batch_size, shuffle=False)
    test_loader = DataLoader(StockDataset(X_test_norm, y_test), batch_size=args.batch_size, shuffle=False)

    n_features = X_train_norm.shape[2]
    n_outputs = len(horizon_names)

    baseline_results = BaselineModels.evaluate_baselines(y_train, y_test, horizon_names)

    candidates = {
        "CNN": CNNRegressor(n_features=n_features, seq_length=args.seq_length, n_outputs=n_outputs),
        "TDNN": TDNNRegressor(n_features=n_features, seq_length=args.seq_length, n_outputs=n_outputs),
        "CNN_LSTM": CNNLSTMRegressor(n_features=n_features, seq_length=args.seq_length, n_outputs=n_outputs)
    }

    neural_results = {}
    trained_models = {}

    for model_type, model in candidates.items():
        trained_model, history, best_val = Trainer.train_model(
            model=model,
            train_loader=train_loader,
            val_loader=val_loader,
            epochs=args.epochs,
            model_name=f"{ticker} {model_type}",
            lr=args.learning_rate,
            weight_decay=args.weight_decay,
            horizon_weights=horizon_weights,
            patience=args.patience
        )
        result = Evaluator.evaluate_model(trained_model, test_loader, horizon_names, model_name=f"{ticker} {model_type}")
        result["history"] = history
        result["best_val_loss"] = best_val
        neural_results[model_type] = result
        trained_models[model_type] = trained_model

    best_model_type = min(neural_results.keys(), key=lambda m: neural_results[m]["metrics"]["summary"]["mean_rmse_log_return"])
    best_model = trained_models[best_model_type]
    best_result = neural_results[best_model_type]

    print("\n" + "=" * 80)
    print(f"BEST MODEL: {best_model_type}")
    print("=" * 80)

    forecast_df, origin_price, origin_date = predict_daily_path_from_origin_date(
        model=best_model,
        df_features=df_features,
        feature_names=feature_names,
        scaler=scaler,
        seq_length=args.seq_length,
        origin_date=args.train_end,
        horizons=horizons,
        n_simulations=args.mc_simulations
    )

    forecast_csv = output_dir / "daily_forecast_path_from_training_end.csv"
    forecast_df.to_csv(forecast_csv, index=False)

    plot_path = ForecastPlotter.plot_daily_forecast_curve_from_origin(
        full_dates=df_features.index,
        full_close_prices=df_features["Close"].values,
        forecast_df=forecast_df,
        origin_date=origin_date,
        origin_price=origin_price,
        ticker=ticker,
        train_start=args.train_start,
        train_end=args.train_end,
        save_dir=str(plot_dir)
    )

    milestone_csv = output_dir / "backtest_predictions_milestones.csv"
    milestone_df = build_milestone_backtest_csv(ticker, dates_test, prices_test, best_result["targets"], best_result["predictions"], horizon_names, horizons)
    milestone_df.to_csv(milestone_csv, index=False)

    metrics_package = {
        "ticker": ticker,
        "best_model_type": best_model_type,
        "max_horizon_days": args.max_horizon_days,
        "baseline_summary": {k: v["summary"] for k, v in baseline_results.items()},
        "neural_summary": {k: v["metrics"]["summary"] for k, v in neural_results.items()},
        "milestone_plot": plot_path,
        "forecast_csv": str(forecast_csv),
        "milestone_backtest_csv": str(milestone_csv)
    }

    metrics_path = output_dir / "metrics.json"
    with open(metrics_path, "w", encoding="utf-8") as f:
        json.dump(metrics_package, f, indent=4)

    model_path = Path("saved_models") / f"{ticker}_best_model.pth"
    Trainer.save_model_package(
        model=best_model,
        save_path=str(model_path),
        model_type=best_model_type,
        ticker=ticker,
        download_start=args.download_start,
        download_end=args.download_end,
        train_start=args.train_start,
        train_end=args.train_end,
        val_end=args.val_end,
        test_end=args.test_end,
        seq_length=args.seq_length,
        horizons=horizons,
        feature_names=feature_names,
        scaler=scaler,
        n_features=n_features,
        n_outputs=n_outputs,
        evaluation_results=best_result["metrics"]["summary"],
        prediction_download_period=args.prediction_download_period,
        mc_simulations=args.mc_simulations
    )

    print("\nSaved outputs:")
    print(f"  Model: {model_path}")
    print(f"  Main plot: {plot_path}")
    print(f"  Forecast CSV: {forecast_csv}")
    print(f"  Backtest milestone CSV: {milestone_csv}")
    print(f"  Metrics: {metrics_path}")

    return metrics_package


def parse_args():
    parser = argparse.ArgumentParser(description="Train daily-path stock forecasting models")
    parser.add_argument("--tickers", nargs="+", default=["SPY"])
    parser.add_argument("--download-start", default="2010-01-01")
    parser.add_argument("--download-end", default=None)
    parser.add_argument("--train-start", default="2012-01-01")
    parser.add_argument("--train-end", default="2018-12-31")
    parser.add_argument("--val-end", default="2019-12-31")
    parser.add_argument("--test-end", default="2020-06-30")
    parser.add_argument("--seq-length", type=int, default=60)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--learning-rate", type=float, default=0.001)
    parser.add_argument("--weight-decay", type=float, default=1e-5)
    parser.add_argument("--patience", type=int, default=15)
    parser.add_argument("--mc-simulations", type=int, default=50)
    parser.add_argument("--max-horizon-days", type=int, default=1512)
    parser.add_argument("--prediction-download-period", default="7y")
    return parser.parse_args()


def main():
    args = parse_args()
    results = {}
    for ticker in args.tickers:
        try:
            results[ticker.upper()] = train_one_ticker(args, ticker.upper())
        except Exception as e:
            print("\n" + "=" * 100)
            print(f"ERROR TRAINING {ticker}")
            print("=" * 100)
            print(e)
    print("\nTraining complete.")
    return results


if __name__ == "__main__":
    main()
