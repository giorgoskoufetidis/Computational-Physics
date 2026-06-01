import yfinance as yf 
import pandas_ta as ta
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import numpy as np 
import pandas as pd 
class StockDataDownloader:
    @staticmethod    
    def download_stock_data(ticker, period="5y", progress=False):
        try:
            df = yf.download(ticker, period=period)
            print(f" {len(df)} trading days of data downloaded for {ticker} with period {period}.")
            print(f"  Date range: {df.index[0].date()} to {df.index[-1].date()}")
            return df
        except Exception as e:
            print(f"Error downloading data for {ticker}: {e}")
            return None
    @staticmethod
    def calculate_technical_indicators(df):
        df = df.copy()
        if isinstance(df.columns, pd.MultiIndex):
            df.columns = df.columns.get_level_values(0)
        # Add technical indicators using pandas_ta    
        df.ta.rsi(length=14, append=True)
        df.ta.macd(fast=12, slow=26, signal=9, append=True)
        df.ta.bbands(length=20, std=2, append=True)
        df.ta.stoch(length=14, append=True)
        df.ta.atr(length=14, append=True)
       # Calculate custom features
        df['Daily_Return'] = df['Close'].pct_change()
        df['ProcessLookupError'] = df['Close'].diff()
        df['Volume_Change'] = df['Volume'].pct_change()
        df['High_Low_Ratio'] = df['High'] / (df['Low']+ 1e-8)
        df['Close_Open_Ratio'] = df['Close'] / (df['Open'] + 1e-8)
        
        #Rolling statistics
        df['Price_MA10'] = df['Close'].rolling(window=10).mean()
        df['Price_MA20'] = df['Close'].rolling(window=20).mean()
        df['Volume_MA'] = df['Volume'].rolling(window=20).mean()
        return df.dropna()
    def create_target_variable(df, shift_days=1):
        df = df.copy()
        next_return = df['Daily_Return'].shift(-shift_days)
        target = np.zeros(len(df))
        target[next_return < -0.01] = 0 # Down
        target[(next_return >= -0.01) & (next_return <= 0.01)] = 1 # Neutral
        target[next_return > 0.01] = 2 # Up
        df['Target'] = target 
        return df
    
class FeatureEngineer:
    @staticmethod
    def select_features(df):
        exclude_cols = ['Target', 'Dividens', 'Stock Splits']
        features = [col for col in df.columns if col not in exclude_cols and 
                    df[col].dtype in ['float64', 'float32', 'int64', 'int32']]
        if len(features) < 10:
            if 'Close' in df.columns:
                return df[['Open', 'High', 'Low', 'Close', 'Volume']].values
        return df[features].values
    
    @staticmethod
    def normalize_features(X_train, X_val, X_test):
        scaler = MinMaxScaler()
        X_trained_scaled = scaler.fit_transform(X_train)
        X_val_scaled = scaler.transform(X_val)
        X_test_scaled = scaler.transform(X_test)
        return X_trained_scaled, X_val_scaled, X_test_scaled
    
    @staticmethod
    def create_sequences(X, y, seq_length=20, stride=1):
        """
        Create sequences for time series models
        
        Args:
            X: Feature matrix
            y: Target vector
            seq_length: Sequence length (days)
            stride: Stride between sequences
        
        Returns:
            X_seq: (n_samples, seq_length, n_features)
            y_seq: Target values
        """
        
        X_seq, y_seq = [], []
        
        for i in range(0, len(X) - seq_length, stride):
            X_seq.append(X[i:i + seq_length])
            y_seq.append(y[i + seq_length - 1])  
        return np.array(X_seq), np.array(y_seq)