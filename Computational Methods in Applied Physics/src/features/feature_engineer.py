import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler


class FeatureEngineer:
    @staticmethod
    def select_features(df):
        exclude_cols = ["Dividends", "Dividens", "Stock Splits"]
        exclude_cols += [col for col in df.columns if col.startswith("Target_")]

        features = [
            col for col in df.columns
            if col not in exclude_cols
            and df[col].dtype in ["float64", "float32", "int64", "int32"]
        ]

        return df[features].values, features

    @staticmethod
    def select_features_by_names(df, feature_names):
        missing = [feature for feature in feature_names if feature not in df.columns]
        if missing:
            raise ValueError(f"Missing features in new data: {missing}")
        return df[feature_names].values

    @staticmethod
    def extract_targets(df, horizon_names):
        target_cols = [f"Target_{h}" for h in horizon_names]
        missing = [col for col in target_cols if col not in df.columns]
        if missing:
            raise ValueError(f"Missing target columns: {missing}")
        return df[target_cols].values, target_cols

    @staticmethod
    def normalize_features(X_train, X_val, X_test):
        scaler = MinMaxScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_val_scaled = scaler.transform(X_val)
        X_test_scaled = scaler.transform(X_test)
        return X_train_scaled, X_val_scaled, X_test_scaled, scaler

    @staticmethod
    def create_sequences(X, y, dates, seq_length=60, stride=1):
        X_seq = []
        y_seq = []
        date_seq = []

        for i in range(0, len(X) - seq_length + 1, stride):
            end_idx = i + seq_length - 1
            X_seq.append(X[i:i + seq_length])
            y_seq.append(y[end_idx])
            date_seq.append(dates[end_idx])

        return np.array(X_seq), np.array(y_seq), np.array(date_seq)

    @staticmethod
    def split_by_date_masks(X, y, dates, current_prices, train_start, train_end, val_end, test_end=None):
        dates = pd.to_datetime(dates)
        train_start = pd.Timestamp(train_start)
        train_end = pd.Timestamp(train_end)
        val_end = pd.Timestamp(val_end)
        test_end = pd.Timestamp(test_end) if test_end is not None else None

        train_mask = (dates >= train_start) & (dates <= train_end)
        val_mask = (dates > train_end) & (dates <= val_end)

        if test_end is not None:
            test_mask = (dates > val_end) & (dates <= test_end)
        else:
            test_mask = dates > val_end

        return (
            X[train_mask], X[val_mask], X[test_mask],
            y[train_mask], y[val_mask], y[test_mask],
            dates[train_mask], dates[val_mask], dates[test_mask],
            current_prices[train_mask], current_prices[val_mask], current_prices[test_mask]
        )
