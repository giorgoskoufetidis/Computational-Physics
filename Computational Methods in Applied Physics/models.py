import torch
import torch.nn as nn
import torch.nn.functional as F


class CNNClassifier(nn.Module):
    def __init__(self, n_features=16, seq_length=20, n_classes=3):
        super(CNNClassifier, self).__init__()
        
        # Encoder: Extract features
        self.conv1 = nn.Conv1d(n_features, 32, kernel_size=3, padding=1)
        self.bn1 = nn.BatchNorm1d(32)
        self.conv2 = nn.Conv1d(32, 64, kernel_size=3, padding=1)
        self.bn2 = nn.BatchNorm1d(64)
        self.conv3 = nn.Conv1d(64, 128, kernel_size=3, padding=1)
        self.bn3 = nn.BatchNorm1d(128)
        
        self.pool = nn.MaxPool1d(2)
        self.dropout = nn.Dropout(0.3)
        
        # Calculate size after convolutions
        size_after_conv = seq_length // 8 * 128
        
        # Classifier
        self.fc1 = nn.Linear(size_after_conv, 64)
        self.fc2 = nn.Linear(64, 32)
        self.fc3 = nn.Linear(32, n_classes)
    
    def forward(self, x):
        # (batch, seq_len, features) -> (batch, features, seq_len)
        x = x.transpose(1, 2)
        
        # Convolutional layers
        x = F.relu(self.bn1(self.conv1(x)))
        x = self.pool(x)
        x = self.dropout(x)
        
        x = F.relu(self.bn2(self.conv2(x)))
        x = self.pool(x)
        x = self.dropout(x)
        
        x = F.relu(self.bn3(self.conv3(x)))
        x = self.pool(x)
        x = self.dropout(x)
        
        # Flatten
        x = x.view(x.size(0), -1)
        
        # Classification
        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        x = F.relu(self.fc2(x))
        x = self.fc3(x)
        
        return x
    
class TDNNClassifier(nn.Module):
    """TDNN for stock direction classification"""
    
    def __init__(self, n_features=16, seq_length=20, n_classes=3, time_delays=(1, 2, 4, 8)):
        super(TDNNClassifier, self).__init__()
        
        self.seq_length = seq_length
        self.time_delays = time_delays
        
        # Number of delayed features
        n_delayed = 1 + len(time_delays)
        
        # Dense layers
        self.fc1 = nn.Linear(seq_length * n_features * n_delayed, 128)
        self.bn1 = nn.BatchNorm1d(128)
        self.dropout1 = nn.Dropout(0.3)
        
        self.fc2 = nn.Linear(128, 64)
        self.bn2 = nn.BatchNorm1d(64)
        self.dropout2 = nn.Dropout(0.3)
        
        self.fc3 = nn.Linear(64, 32)
        self.fc4 = nn.Linear(32, n_classes)
    
    def forward(self, x):
        # x: (batch, seq_len, features)
        batch_size = x.size(0)
        n_features = x.size(2)
        
        # Create time-delayed versions
        delayed_features = [x]
        
        for delay in self.time_delays:
            if delay > 0:
                # Pad and shift
                delayed = torch.cat([
                    torch.zeros(batch_size, delay, n_features, device=x.device),
                    x[:, :-delay, :]
                ], dim=1)
                delayed_features.append(delayed)
        
        # Concatenate all
        x_delayed = torch.cat(delayed_features, dim=2)  # (batch, seq_len, features*(1+delays))
        
        # Flatten
        x = x_delayed.view(batch_size, -1)
        
        # Dense layers
        x = F.relu(self.bn1(self.fc1(x)))
        x = self.dropout1(x)
        
        x = F.relu(self.bn2(self.fc2(x)))
        x = self.dropout2(x)
        
        x = F.relu(self.fc3(x))
        x = self.fc4(x)
        
        return x