import torch
from src.models.neural_models import CNNRegressor, TDNNRegressor, CNNLSTMRegressor


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class ModelLoader:
    @staticmethod
    def build_model_from_package(package):
        model_type = package["model_type"]
        n_features = package["n_features"]
        seq_length = package["seq_length"]
        n_outputs = package["n_outputs"]

        if model_type == "CNN":
            model = CNNRegressor(n_features=n_features, seq_length=seq_length, n_outputs=n_outputs)
        elif model_type == "TDNN":
            model = TDNNRegressor(n_features=n_features, seq_length=seq_length, n_outputs=n_outputs)
        elif model_type == "CNN_LSTM":
            model = CNNLSTMRegressor(n_features=n_features, seq_length=seq_length, n_outputs=n_outputs)
        else:
            raise ValueError(f"Unknown model type: {model_type}")

        model.load_state_dict(package["model_state_dict"])
        model = model.to(device)
        model.eval()
        return model

    @staticmethod
    def load_package(path):
        package = torch.load(path, map_location=device, weights_only=False)
        model = ModelLoader.build_model_from_package(package)
        return package, model
