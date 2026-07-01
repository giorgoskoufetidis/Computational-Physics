import numpy as np
from training.evaluator import Evaluator


class BaselineModels:
    @staticmethod
    def zero_return_baseline(y_test):
        return np.zeros_like(y_test)

    @staticmethod
    def mean_return_baseline(y_train, y_test):
        mean_returns = np.mean(y_train, axis=0)
        return np.tile(mean_returns, (len(y_test), 1))

    @staticmethod
    def evaluate_baselines(y_train, y_test, horizon_names):
        zero = BaselineModels.zero_return_baseline(y_test)
        mean = BaselineModels.mean_return_baseline(y_train, y_test)
        results = {
            "Zero Return Baseline": Evaluator.calculate_metrics(y_test, zero, horizon_names),
            "Mean Return Baseline": Evaluator.calculate_metrics(y_test, mean, horizon_names)
        }
        print("\n" + "=" * 80)
        print("Baseline comparison")
        print("=" * 80)
        for name, res in results.items():
            print(f"{name}: RMSE={res['summary']['mean_rmse_log_return']:.6f}, Direction={res['summary']['mean_directional_accuracy']:.4f}")
        return results
