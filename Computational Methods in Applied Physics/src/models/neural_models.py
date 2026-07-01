import torch
import torch.nn as nn
import torch.nn.functional as F


class CNNRegressor(nn.Module):
    def __init__(self, n_features, seq_length, n_outputs):
        super().__init__()
        self.conv1 = nn.Conv1d(n_features, 32, kernel_size=3, padding=1)
        self.bn1 = nn.BatchNorm1d(32)
        self.conv2 = nn.Conv1d(32, 64, kernel_size=3, padding=1)
        self.bn2 = nn.BatchNorm1d(64)
        self.conv3 = nn.Conv1d(64, 128, kernel_size=3, padding=1)
        self.bn3 = nn.BatchNorm1d(128)
        self.pool = nn.MaxPool1d(2)
        self.dropout = nn.Dropout(0.30)

        size_after_conv = (seq_length // 8) * 128
        if size_after_conv <= 0:
            raise ValueError("seq_length too small for this CNN architecture")

        self.fc1 = nn.Linear(size_after_conv, 256)
        self.fc2 = nn.Linear(256, 128)
        self.fc3 = nn.Linear(128, n_outputs)

    def forward(self, x):
        x = x.transpose(1, 2)
        x = self.pool(F.relu(self.bn1(self.conv1(x))))
        x = self.dropout(x)
        x = self.pool(F.relu(self.bn2(self.conv2(x))))
        x = self.dropout(x)
        x = self.pool(F.relu(self.bn3(self.conv3(x))))
        x = self.dropout(x)
        x = x.reshape(x.size(0), -1)
        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        x = F.relu(self.fc2(x))
        return self.fc3(x)


class TDNNRegressor(nn.Module):
    def __init__(self, n_features, seq_length, n_outputs, time_delays=(1, 2, 4, 8)):
        super().__init__()
        self.time_delays = time_delays
        n_versions = 1 + len(time_delays)
        input_size = seq_length * n_features * n_versions

        self.fc1 = nn.Linear(input_size, 512)
        self.norm1 = nn.LayerNorm(512)
        self.drop1 = nn.Dropout(0.30)
        self.fc2 = nn.Linear(512, 256)
        self.norm2 = nn.LayerNorm(256)
        self.drop2 = nn.Dropout(0.30)
        self.fc3 = nn.Linear(256, 128)
        self.norm3 = nn.LayerNorm(128)
        self.drop3 = nn.Dropout(0.20)
        self.fc4 = nn.Linear(128, n_outputs)

    def forward(self, x):
        batch, seq_len, features = x.shape
        delayed_features = [x]

        for delay in self.time_delays:
            if delay >= seq_len:
                continue
            delayed = torch.cat([
                torch.zeros(batch, delay, features, device=x.device, dtype=x.dtype),
                x[:, :-delay, :]
            ], dim=1)
            delayed_features.append(delayed)

        x = torch.cat(delayed_features, dim=2)
        x = x.reshape(batch, -1)
        x = F.relu(self.norm1(self.fc1(x)))
        x = self.drop1(x)
        x = F.relu(self.norm2(self.fc2(x)))
        x = self.drop2(x)
        x = F.relu(self.norm3(self.fc3(x)))
        x = self.drop3(x)
        return self.fc4(x)


class CNNLSTMRegressor(nn.Module):
    def __init__(self, n_features, seq_length, n_outputs, lstm_hidden=64, lstm_layers=1):
        super().__init__()
        self.conv1 = nn.Conv1d(n_features, 32, kernel_size=3, padding=1)
        self.bn1 = nn.BatchNorm1d(32)
        self.conv2 = nn.Conv1d(32, 64, kernel_size=3, padding=1)
        self.bn2 = nn.BatchNorm1d(64)
        self.dropout = nn.Dropout(0.30)

        self.lstm = nn.LSTM(
            input_size=64,
            hidden_size=lstm_hidden,
            num_layers=lstm_layers,
            batch_first=True
        )

        self.fc1 = nn.Linear(lstm_hidden, 256)
        self.fc2 = nn.Linear(256, n_outputs)

    def forward(self, x):
        x = x.transpose(1, 2)
        x = F.relu(self.bn1(self.conv1(x)))
        x = self.dropout(x)
        x = F.relu(self.bn2(self.conv2(x)))
        x = self.dropout(x)
        x = x.transpose(1, 2)
        out, _ = self.lstm(x)
        x = out[:, -1, :]
        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        return self.fc2(x)
