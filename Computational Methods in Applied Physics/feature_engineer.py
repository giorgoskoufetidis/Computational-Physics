import numpy as np
from sklearn.preprocessing import MinMaxScaler

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
        return X_trained_scaled, X_val_scaled, X_test_scaled, scaler
    
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
            y_seq.append(y[i + seq_length - 1])  # Target at end of sequence
        
        return np.array(X_seq), np.array(y_seq)