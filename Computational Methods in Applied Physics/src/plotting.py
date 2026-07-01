import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class ForecastPlotter:
    @staticmethod
    def plot_daily_forecast_curve_from_origin(
        full_dates,
        full_close_prices,
        forecast_df,
        origin_date,
        origin_price,
        ticker,
        train_start,
        train_end,
        save_dir="outputs/plots"
    ):
        os.makedirs(save_dir, exist_ok=True)

        full_dates = pd.to_datetime(full_dates)
        full_close_prices = np.array(full_close_prices)
        origin_date = pd.Timestamp(origin_date)
        train_start = pd.Timestamp(train_start)
        train_end = pd.Timestamp(train_end)

        forecast_df = forecast_df.copy()
        forecast_df["forecast_date"] = pd.to_datetime(forecast_df["forecast_date"])

        plt.figure(figsize=(17, 8))

        plt.plot(full_dates, full_close_prices, label="Actual stock price", linewidth=2)

        plt.axvspan(train_start, train_end, alpha=0.12, label="Training period")
        plt.axvline(origin_date, linestyle="--", linewidth=2, label="Last training date / forecast origin")

        path_dates = [origin_date] + forecast_df["forecast_date"].tolist()
        path_prices = [origin_price] + forecast_df["expected_price"].tolist()

        plt.plot(path_dates, path_prices, linestyle="--", linewidth=2.2, label="Model predicted daily future path")

        plt.fill_between(
            forecast_df["forecast_date"],
            forecast_df["lower_price"],
            forecast_df["upper_price"],
            alpha=0.15,
            label="95% forecast range"
        )

        milestones = {"1d": 1, "1m": 21, "1y": 252, "5y": 1260, "6y": 1512}
        for label, day in milestones.items():
            row = forecast_df[forecast_df["trading_day_ahead"] == day]
            if len(row) == 0:
                continue
            row = row.iloc[0]
            plt.scatter(row["forecast_date"], row["expected_price"], s=110, marker=".")
            plt.annotate(
                label,
                xy=(row["forecast_date"], row["expected_price"]),
                xytext=(8, 8),
                textcoords="offset points",
                fontsize=11
            )

        plt.title(f"{ticker} - Daily Forecast Path from Last Training Date")
        plt.xlabel("Date")
        plt.ylabel("Price")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()

        predict_start_date = pd.Timestamp(origin_date).strftime("%Y%m%d")

        predict_end_date = pd.to_datetime(
            forecast_df["forecast_date"].iloc[-1]
        ).strftime("%Y%m%d")

        file_path = os.path.join(
            save_dir,
            f"{ticker}_forecast_from_{predict_start_date}_to_{predict_end_date}_daily_path.png"
)
        plt.savefig(file_path, dpi=150)
        plt.close()
        print(f"Saved plot: {file_path}")
        return file_path
