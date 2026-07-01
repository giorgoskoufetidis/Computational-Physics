import os
import sys
import subprocess
import threading
import traceback
import tkinter as tk
from pathlib import Path
from tkinter import ttk, messagebox

from src.predictor import ForecastPredictor


class PredictionGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Stock / ETF Forecast Predictor")
        self.root.geometry("920x720")

        self.project_root = Path(__file__).resolve().parents[1]
        self.saved_models_folder = self.project_root / "saved_models"
        self.available_models = {}
        self.predictor = None

        self.create_widgets()
        self.scan_saved_models()

    def create_widgets(self):
        main = ttk.Frame(self.root, padding=15)
        main.pack(fill="both", expand=True)

        ttk.Label(main, text="Stock / ETF Forecast Predictor", font=("Arial", 20, "bold")).pack(pady=10)
        ttk.Label(main, text="Loads saved models and shows milestone predictions from the latest market date.").pack(pady=5)

        model_frame = ttk.LabelFrame(main, text="Available Trained Models", padding=10)
        model_frame.pack(fill="x", pady=10)

        ttk.Label(model_frame, text="Ticker:").grid(row=0, column=0, padx=5, pady=5)
        self.ticker_var = tk.StringVar()
        self.ticker_combo = ttk.Combobox(model_frame, textvariable=self.ticker_var, values=[], width=25, state="readonly")
        self.ticker_combo.grid(row=0, column=1, padx=5, pady=5)

        ttk.Button(model_frame, text="Refresh Models", command=self.scan_saved_models).grid(row=0, column=2, padx=5, pady=5)
        ttk.Button(model_frame, text="Load Model", command=self.load_selected_model).grid(row=0, column=3, padx=5, pady=5)

        info_frame = ttk.LabelFrame(main, text="Model Information", padding=10)
        info_frame.pack(fill="x", pady=10)
        self.model_info_label = ttk.Label(info_frame, text="No model loaded.", justify="left")
        self.model_info_label.pack(anchor="w")

        button_frame = ttk.Frame(main)
        button_frame.pack(fill="x", pady=10)

        self.predict_button = ttk.Button(button_frame, text="Predict Latest", command=self.start_prediction, state="disabled")
        self.predict_button.pack(side="left", padx=5)
        ttk.Button(button_frame, text="Open Plots Folder", command=self.open_plots_folder).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Clear Output", command=self.clear_output).pack(side="left", padx=5)

        self.status_label = ttk.Label(button_frame, text="Status: No model loaded", font=("Arial", 10, "bold"))
        self.status_label.pack(side="right", padx=5)

        output_frame = ttk.LabelFrame(main, text="Prediction Output", padding=10)
        output_frame.pack(fill="both", expand=True, pady=10)
        self.output_text = tk.Text(output_frame, wrap="word", height=24, font=("Consolas", 10))
        self.output_text.pack(fill="both", expand=True)

    def scan_saved_models(self):
        self.available_models = {}
        self.saved_models_folder.mkdir(parents=True, exist_ok=True)

        for model_path in self.saved_models_folder.glob("*.pth"):
            name = model_path.name
            ticker = name.replace("_best_model.pth", "") if name.endswith("_best_model.pth") else model_path.stem
            self.available_models[ticker.upper()] = str(model_path)

        tickers = sorted(self.available_models.keys())
        self.ticker_combo["values"] = tickers

        if tickers:
            self.ticker_var.set(tickers[0])
            self.status_label.config(text=f"Status: Found {len(tickers)} model(s)")
            self.model_info_label.config(text=f"Found models in:\n{self.saved_models_folder}\n\nAvailable: {', '.join(tickers)}")
        else:
            self.ticker_var.set("")
            self.status_label.config(text="Status: No trained models found")
            self.model_info_label.config(text=f"No models found. Expected folder:\n{self.saved_models_folder}")

        self.predictor = None
        self.predict_button.config(state="disabled")

    def load_selected_model(self):
        try:
            ticker = self.ticker_var.get().strip().upper()
            if not ticker:
                raise ValueError("No ticker selected")

            model_path = self.available_models[ticker]
            self.predictor = ForecastPredictor(model_path)
            package = self.predictor.package

            info = (
                f"Model loaded successfully\n"
                f"File: {model_path}\n"
                f"Ticker: {package.get('ticker', ticker)}\n"
                f"Model type: {package.get('model_type', 'N/A')}\n"
                f"Train period: {package.get('train_start', 'N/A')} to {package.get('train_end', 'N/A')}\n"
                f"Validation end: {package.get('val_end', 'N/A')}\n"
                f"Test end: {package.get('test_end', 'N/A')}\n"
                f"Sequence length: {package.get('seq_length', 'N/A')}\n"
                f"Outputs: {package.get('n_outputs', 'N/A')} daily horizons\n"
                f"Features: {len(package.get('feature_names', []))}"
            )
            self.model_info_label.config(text=info)
            self.status_label.config(text=f"Status: Loaded {ticker}")
            self.predict_button.config(state="normal")
            self.clear_output()
            self.write_output(info)
        except Exception as e:
            messagebox.showerror("Load Error", str(e))

    def start_prediction(self):
        threading.Thread(target=self.predict, daemon=True).start()

    def predict(self):
        try:
            if self.predictor is None:
                raise RuntimeError("Load a model first")

            self.predict_button.config(state="disabled")
            self.status_label.config(text="Status: Predicting...")
            self.clear_output()

            result = self.predictor.predict_latest_path()
            milestones = ForecastPredictor.get_milestones_from_forecast_df(result["forecast_path"])

            self.write_output("=" * 80)
            self.write_output(f"{result['ticker']} LATEST FORECAST")
            self.write_output("=" * 80)
            self.write_output(f"Model type:    {result['model_type']}")
            self.write_output(f"Current date:  {result['current_date']}")
            self.write_output(f"Current price: ${result['current_price']:.2f}")
            self.write_output("")

            for label, row in milestones.items():
                self.write_output(f"{label} forecast")
                self.write_output("-" * 60)
                self.write_output(f"Forecast date:   {row['forecast_date']}")
                self.write_output(f"Expected price:  ${row['expected_price']:.2f}")
                self.write_output(f"Expected return: {row['expected_return_pct']:.2f}%")
                self.write_output(f"95% price range: ${row['lower_price']:.2f} to ${row['upper_price']:.2f}")
                self.write_output("")

            self.status_label.config(text="Status: Prediction complete")
        except Exception:
            self.write_output(traceback.format_exc())
            self.status_label.config(text="Status: Error")
        finally:
            self.predict_button.config(state="normal")

    def open_plots_folder(self):
        ticker = self.ticker_var.get().strip().upper()
        if not ticker:
            messagebox.showinfo("No ticker", "Select a ticker first")
            return

        plots_folder = self.project_root / "outputs" / ticker / "plots"
        if not plots_folder.exists():
            messagebox.showinfo("Plots not found", f"No plots found in:\n{plots_folder}\n\nRun training or backtesting first.")
            return

        if sys.platform.startswith("win"):
            os.startfile(plots_folder)
        elif sys.platform == "darwin":
            subprocess.call(["open", str(plots_folder)])
        else:
            subprocess.call(["xdg-open", str(plots_folder)])

    def write_output(self, text):
        self.output_text.insert(tk.END, str(text) + "\n")
        self.output_text.see(tk.END)

    def clear_output(self):
        self.output_text.delete("1.0", tk.END)


def main():
    root = tk.Tk()
    PredictionGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
