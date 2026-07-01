import numpy as np
import pandas as pd
import yfinance as yf
import pandas_ta as ta


class StockDataDownloader:
    """
    Downloads stock/ETF data and creates technical indicators.
    """

    @staticmethod
    def download_stock_data(ticker, period="10y", start=None, end=None, progress=False):
        try:
            if start is not None or end is not None:
                df = yf.download(
                    ticker,
                    start=start,
                    end=end,
                    auto_adjust=False,
                    progress=progress
                )
            else:
                df = yf.download(
                    ticker,
                    period=period,
                    auto_adjust=False,
                    progress=progress
                )

            if df is None or len(df) == 0:
                raise ValueError("No data downloaded")

            if isinstance(df.columns, pd.MultiIndex):
                df.columns = df.columns.get_level_values(0)

            print(f"{len(df)} trading days downloaded for {ticker}")
            print(f"Date range: {df.index[0].date()} to {df.index[-1].date()}")

            return df

        except Exception as e:
            print(f"Error downloading data for {ticker}: {e}")
            return None

    @staticmethod
    def calculate_technical_indicators(df):
        df = df.copy()

        if isinstance(df.columns, pd.MultiIndex):
            df.columns = df.columns.get_level_values(0)

        required = ["Open", "High", "Low", "Close", "Volume"]
        for col in required:
            if col not in df.columns:
                raise ValueError(f"Missing required column: {col}")

        # Technical indicators from pandas_ta
        df.ta.rsi(length=14, append=True)
        df.ta.macd(fast=12, slow=26, signal=9, append=True)
        df.ta.bbands(length=20, std=2, append=True)
        df.ta.stoch(length=14, append=True)
        df.ta.atr(length=14, append=True)

        # Custom features
        df["Daily_Return"] = df["Close"].pct_change()
        df["Log_Return"] = np.log(df["Close"] / df["Close"].shift(1))
        df["Price_Change"] = df["Close"].diff()
        df["Volume_Change"] = df["Volume"].pct_change()

        df["High_Low_Ratio"] = df["High"] / (df["Low"] + 1e-8)
        df["Close_Open_Ratio"] = df["Close"] / (df["Open"] + 1e-8)

        df["Price_MA10"] = df["Close"].rolling(window=10).mean()
        df["Price_MA20"] = df["Close"].rolling(window=20).mean()
        df["Price_MA50"] = df["Close"].rolling(window=50).mean()
        df["Volume_MA20"] = df["Volume"].rolling(window=20).mean()

        df["Close_vs_MA10"] = df["Close"] / (df["Price_MA10"] + 1e-8) - 1
        df["Close_vs_MA20"] = df["Close"] / (df["Price_MA20"] + 1e-8) - 1
        df["Close_vs_MA50"] = df["Close"] / (df["Price_MA50"] + 1e-8) - 1

        df = df.replace([np.inf, -np.inf], np.nan)
        return df.dropna()

    @staticmethod
    def create_forecast_targets(df, horizons):
        """
        Creates one target column per future trading day.

        Target_d0001 = log(Close[t+1] / Close[t])
        Target_d0021 = log(Close[t+21] / Close[t])
        Target_d1260 = log(Close[t+1260] / Close[t])
        """
        df = df.copy()

        target_data = {}
        close = df["Close"]

        for name, days in horizons.items():
            future_close = close.shift(-days)
            target_data[f"Target_{name}"] = np.log(future_close / close)

        target_df = pd.DataFrame(target_data, index=df.index)
        df = pd.concat([df, target_df], axis=1)
        df = df.replace([np.inf, -np.inf], np.nan)

        return df.dropna()
