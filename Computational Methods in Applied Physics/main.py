
from StockDataDownloader import StockDataDownloader
from feature_engineer import FeatureEngineer
from models import CNNClassifier, TDNNClassifier
from train_classifier import TrainClassifier
from sklearn.model_selection import train_test_split
import torch
from torch.utils.data import DataLoader
from stock_dataset import StockDataset
import numpy as np
def main():
    print("\n" + "="*80)
    print("PHASE 1: DOWNLOAD REAL STOCK DATA")
    print("="*80)
    
    tickers = 'SPY'
    all_data = {}
    
    for ticker in tickers:
        df = StockDataDownloader.download_stock_data(ticker, period="5y")
        if df is not None:
            df = StockDataDownloader.calculate_technical_indicators(df)
            df = StockDataDownloader.create_target_variable(df)
            all_data[ticker] = df
            print(f"   Features calculated: {df.shape[1]} columns, {df.shape[0]} rows")
    print("\n" + "="*80)
    print("PHASE 2: FEATURE ENGINEERING")
    print("="*80)
    
    for ticker, df in all_data.items():
        X = FeatureEngineer.select_features(df)
        y = all_data[ticker].values
    print(f"Features shape: {X.shape}")
    
    seq_length = 20
    X_seq, y_seq = FeatureEngineer.create_sequences(X, y, seq_length=seq_length, stride=1)
    print(f"Sequence shape: {X_seq.shape}, Target shape: {y_seq.shape}")
    X_train, X_test, y_train, y_test = train_test_split(X_seq, y_seq, test_size=0.2, random_state=42)
    X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.2, random_state=42)
    print(f"Train shape: {X_train.shape}, Val shape: {X_val.shape}, Test shape: {X_test.shape}")
  
    

    X_train_norm, X_val_norm, X_test_norm, scaler = FeatureEngineer.normalize_features(
        X_train.reshape(-1, X_train.shape[-1]),
        X_val.reshape(-1, X_val.shape[-1]),
        X_test.reshape(-1, X_test.shape[-1])
    )
    print(f"Normalized shape: Train {X_train_norm.shape}, Val {X_val_norm.shape}, Test {X_test_norm.shape}")
    X_train_norm = X_train_norm.reshape(X_train.shape)
    X_val_norm = X_val_norm.reshape(X_val.shape)
    X_test_norm = X_test_norm.reshape(X_test.shape)
    
    print(f"Normalized Train shape: {X_train_norm.shape}, Val shape: {X_val_norm.shape}, Test shape: {X_test_norm.shape}")
    
    print("\n" + "="*80)
    print("PHASE 3: CREATE DATALOADERS")
    print("="*80)
    batch_size = 32
    train_ds = StockDataset(X_train_norm, y_train)
    val_ds = StockDataset(X_val_norm, y_val)
    test_ds = StockDataset(X_test_norm, y_test)
    
    
    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False)
    test_loader = DataLoader(test_ds, batch_size=batch_size, shuffle=False)
    print(f"✓ DataLoaders created (batch_size={batch_size})")
    
    print("\n" + "="*80)
    print("PHASE 4: TRAIN MODELS")
    print("="*80)
    n_features = X_train_norm.shape[2]
    print(f"✓ Ready to train models with {n_features} features")
    #trainer 
    trainer = TrainClassifier()
    # CNN
    print("\n[4.1] CNN Classifier")
    model_cnn = CNNClassifier(n_features=n_features, seq_length=seq_length, n_classes=3)
    trainer.train_classifier(model_cnn, train_loader, val_loader, epochs=50, model_name="CNN Classifier")
    
    
    return all_data

if __name__ == "__main__":
    results = main()
    # print(results)
    print("\n" + "="*80)
    print("PROJECT COMPLETE")
    print("="*80)