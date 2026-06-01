import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix, classification_report
import numpy as np 
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
class TrainClassifier:
    
    @staticmethod
    def train_classifier(model, train_loader, val_loader, epochs=50, model_name='Model'):
        print(f'\n{'='*80}')
        print(f'Training {model_name}')
        print(f'{'='*80}\n')
        
        model = model.to(device)
        optimizer = optim.Adam(model.parameters(), lr=0.001)
        criterion = nn.CrossEntropyLoss()
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=5)
        best_val_loss = float('inf')
        patience = 10
        patience_counter = 0
        train_losses, val_losses = [], []
        for epoch in range(epochs):
            model.train()
            train_loss = 0.0
            for X_batch, y_batch in train_loader:
                X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                optimizer.zero_grad()
                outputs = model(X_batch)
                print(y_batch)
                loss = criterion(outputs, y_batch)
                loss.backward()
                optimizer.step()
                train_loss += loss.item() 
            train_loss /= len(train_loader)
            train_losses.append(train_loss)
            model.eval()
            val_loss = 0.0
            with torch.no_grad():
                for X_batch, y_batch in val_loader:
                    X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                    outputs = model(X_batch)
                    loss = criterion(outputs, y_batch)
                    val_loss += loss.item()
            val_loss /= len(val_loader)
            val_losses.append(val_loss)
            scheduler.step(val_loss)
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                patience_counter = 0
            else: 
                patience_counter += 1
            
            if patience_counter >= patience:
                print(f'Early stopping at epoch {epoch+1}')
                break
            if (epoch + 1) % 10 == 0 or epoch == 0:
                print(f'Epoch {epoch+1}/{epochs} - Train Loss: {train_loss:.4f} - Val Loss: {val_loss:.4f}')
        return train_losses, val_losses
    
    @staticmethod
    def evaluate_classifier(model, test_loader, model_name="Model"):
        """Evaluate the model on the test set and print metrics"""
        print(f"\n{'='*60}")
        print(f" Evaluating {model_name} on Test Set ")
        print(f"{'='*60}\n")
        model.eval()
        
        all_preds, all_targets, all_probs = [], [], []
        with torch.no_grad():
            for X_batch, y_batch in test_loader:
                X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                outputs = model(X_batch)
                probs = F.softmax(outputs, dim=1)
                preds = torch.argmax(outputs, dim=1)
                
                all_preds.append(preds.cpu().numpy())
                all_targets.append(y_batch.cpu().numpy())
                all_probs.append(probs.cpu().numpy())
        y_pred = np.concatenate(all_preds)
        y_true = np.concatenate(all_targets)
        y_probs = np.concatenate(all_probs) 
            
        accuracy = accuracy_score(y_true, y_pred)
        precision = precision_score(y_true, y_pred, average='weighted', zero_division=0)
        recall = recall_score(y_true, y_pred, average='weighted', zero_division=0)
        f1 = f1_score(y_true, y_pred, average='weighted', zero_division=0)
        print(f"Accuracy:  {accuracy:.4f}")
        print(f"Precision: {precision:.4f}")
        print(f"Recall:    {recall:.4f}")
        print(f"F1-Score:  {f1:.4f}")
        
        print(f"\nConfusion Matrix:")
        cm = confusion_matrix(y_true, y_pred)
        print(cm)
        
        print(f"\nClassification Report:")
        class_names = ['Down', 'Neutral', 'Up']
        print(classification_report(y_true, y_pred, target_names=class_names, zero_division=0))
        
        return {
            'accuracy': accuracy,
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'confusion_matrix': cm,
            'predictions': y_pred,
            'targets': y_true,
            'probabilities': y_probs
        }