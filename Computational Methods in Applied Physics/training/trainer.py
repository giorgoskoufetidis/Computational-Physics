import copy
import os
from datetime import datetime

import torch
import torch.nn as nn
import torch.optim as optim


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class WeightedHuberLoss(nn.Module):
    def __init__(self, horizon_weights=None, delta=1.0):
        super().__init__()
        self.delta = delta
        self.horizon_weights = torch.FloatTensor(horizon_weights) if horizon_weights is not None else None

    def forward(self, predictions, targets):
        error = predictions - targets
        abs_error = torch.abs(error)
        quadratic = torch.minimum(abs_error, torch.tensor(self.delta, device=predictions.device))
        linear = abs_error - quadratic
        loss = 0.5 * quadratic ** 2 + self.delta * linear
        if self.horizon_weights is not None:
            loss = loss * self.horizon_weights.to(predictions.device)
        return loss.mean()


class Trainer:
    @staticmethod
    def train_model(model, train_loader, val_loader, epochs=100, model_name="Model", lr=0.001, weight_decay=1e-5, horizon_weights=None, patience=20):
        print("\n" + "=" * 80)
        print(f"Training {model_name}")
        print("=" * 80)

        model = model.to(device)
        criterion = WeightedHuberLoss(horizon_weights=horizon_weights)
        optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="min", factor=0.5, patience=5)

        best_val = float("inf")
        best_state = copy.deepcopy(model.state_dict())
        wait = 0
        history = {"train_loss": [], "val_loss": []}

        for epoch in range(epochs):
            model.train()
            train_loss = 0.0
            for X_batch, y_batch in train_loader:
                X_batch = X_batch.to(device)
                y_batch = y_batch.to(device)
                optimizer.zero_grad()
                out = model(X_batch)
                loss = criterion(out, y_batch)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()
                train_loss += loss.item()
            train_loss /= len(train_loader)

            model.eval()
            val_loss = 0.0
            with torch.no_grad():
                for X_batch, y_batch in val_loader:
                    X_batch = X_batch.to(device)
                    y_batch = y_batch.to(device)
                    out = model(X_batch)
                    loss = criterion(out, y_batch)
                    val_loss += loss.item()
            val_loss /= len(val_loader)

            history["train_loss"].append(train_loss)
            history["val_loss"].append(val_loss)
            scheduler.step(val_loss)

            if val_loss < best_val:
                best_val = val_loss
                best_state = copy.deepcopy(model.state_dict())
                wait = 0
            else:
                wait += 1

            if epoch == 0 or (epoch + 1) % 10 == 0:
                print(f"Epoch {epoch + 1}/{epochs} | Train {train_loss:.6f} | Val {val_loss:.6f} | Best {best_val:.6f}")

            if wait >= patience:
                print(f"Early stopping at epoch {epoch + 1}")
                break

        model.load_state_dict(best_state)
        print(f"Restored best model with val loss: {best_val:.6f}")
        return model, history, best_val

    @staticmethod
    def save_model_package(model, save_path, model_type, ticker, download_start, download_end, train_start, train_end, val_end, test_end, seq_length, horizons, feature_names, scaler, n_features, n_outputs, evaluation_results=None, prediction_download_period="7y", mc_simulations=100):
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        package = {
            "model_state_dict": model.state_dict(),
            "model_type": model_type,
            "ticker": ticker,
            "download_start": download_start,
            "download_end": download_end,
            "train_start": train_start,
            "train_end": train_end,
            "val_end": val_end,
            "test_end": test_end,
            "seq_length": seq_length,
            "horizons": horizons,
            "feature_names": feature_names,
            "scaler": scaler,
            "n_features": n_features,
            "n_outputs": n_outputs,
            "evaluation_results": evaluation_results,
            "prediction_download_period": prediction_download_period,
            "mc_simulations": mc_simulations,
            "created_at": datetime.now().isoformat()
        }
        torch.save(package, save_path)
        print(f"Saved model package: {save_path}")
