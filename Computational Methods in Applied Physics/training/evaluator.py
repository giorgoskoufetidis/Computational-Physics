import numpy as np
import torch
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class Evaluator:
    @staticmethod
    def collect_predictions(model, loader):
        model = model.to(device)
        model.eval()
        preds = []
        targets = []
        with torch.no_grad():
            for X_batch, y_batch in loader:
                X_batch = X_batch.to(device)
                out = model(X_batch)
                preds.append(out.cpu().numpy())
                targets.append(y_batch.numpy())
        return np.concatenate(targets, axis=0), np.concatenate(preds, axis=0)

    @staticmethod
    def calculate_metrics(y_true, y_pred, horizon_names, print_all=False):
        per_horizon = {}
        maes, rmses, dirs, corrs, r2s = [], [], [], [], []

        for i, h in enumerate(horizon_names):
            true_h = y_true[:, i]
            pred_h = y_pred[:, i]
            mae = mean_absolute_error(true_h, pred_h)
            rmse = np.sqrt(mean_squared_error(true_h, pred_h))
            direction = np.mean(np.sign(true_h) == np.sign(pred_h))
            corr = np.corrcoef(true_h, pred_h)[0, 1] if np.std(true_h) > 0 and np.std(pred_h) > 0 else 0.0
            try:
                r2 = r2_score(true_h, pred_h)
            except Exception:
                r2 = 0.0

            per_horizon[h] = {
                "mae_log_return": float(mae),
                "rmse_log_return": float(rmse),
                "directional_accuracy": float(direction),
                "correlation": float(corr),
                "r2": float(r2)
            }
            maes.append(mae); rmses.append(rmse); dirs.append(direction); corrs.append(corr); r2s.append(r2)

        return {
            "per_horizon": per_horizon,
            "summary": {
                "mean_mae_log_return": float(np.mean(maes)),
                "mean_rmse_log_return": float(np.mean(rmses)),
                "mean_directional_accuracy": float(np.mean(dirs)),
                "mean_correlation": float(np.mean(corrs)),
                "mean_r2": float(np.mean(r2s))
            }
        }

    @staticmethod
    def evaluate_model(model, loader, horizon_names, model_name="Model"):
        y_true, y_pred = Evaluator.collect_predictions(model, loader)
        metrics = Evaluator.calculate_metrics(y_true, y_pred, horizon_names)

        print("\n" + "=" * 80)
        print(f"Evaluation: {model_name}")
        print("=" * 80)
        print("Summary:")
        for k, v in metrics["summary"].items():
            print(f"  {k}: {v:.6f}")

        milestones = {"d0001": "1d", "d0021": "1m", "d0252": "1y", "d1260": "5y", "d1512": "6y"}
        print("Milestone horizons:")
        for h, label in milestones.items():
            if h in metrics["per_horizon"]:
                m = metrics["per_horizon"][h]
                print(f"  {label} ({h}): RMSE={m['rmse_log_return']:.6f}, Direction={m['directional_accuracy']:.4f}, Corr={m['correlation']:.4f}")

        return {"metrics": metrics, "targets": y_true, "predictions": y_pred}
