# Stock / ETF Forecasting App — Daily Future Path Version

This project trains deep-learning models for stock/ETF forecasting and creates a plot that matches this idea:

```text
Train the model on a selected historical period.
Take the last training date as the forecast origin.
Predict a daily future price path from that date.
Highlight 1 day, 1 month, 1 year, and 5 years on the same path.
Compare the path visually against the real stock price when historical data exists.
```

The important difference from the previous version is this:

```text
Old version:
    predicted only 4 points: 1d, 1m, 1y, 5y

New version:
    predicts every trading day from day 1 to day 1260
    then highlights:
        day 1    = 1d
        day 21   = 1m
        day 252  = 1y
        day 1260 = 5y
```

This is a **direct multi-output model**. It does not recursively predict tomorrow and then use tomorrow's prediction to predict the next day. Instead, the model directly learns:

```text
Input:
    last 60 trading days before date t

Output:
    future log returns for:
    t+1, t+2, t+3, ..., t+1260
```

---

## Project structure

```text
stock_forecast_daily_path_project/
│
├── app.py
├── README.md
├── requirements.txt
├── saved_models/
├── outputs/
│
├── src/
│   ├── gui.py
│   ├── predictor.py
│   ├── plotting.py
│   ├── data/
│   │   └── stock_data_downloader.py
│   ├── features/
│   │   └── feature_engineer.py
│   ├── models/
│   │   └── neural_models.py
│   └── utils/
│       └── model_loader.py
│
├── training/
│   ├── train_models.py
│   ├── trainer.py
│   ├── evaluator.py
│   ├── dataset.py
│   └── baselines.py
│
└── backtesting/
    └── backtest_model.py
```

---

## Install

```bash
pip install -r requirements.txt
```

On Linux, Tkinter may need:

```bash
sudo apt-get install python3-tk
```

---

## Train one model

```bash
python training/train_models.py --tickers SPY
```

This trains CNN, TDNN, and CNN-LSTM, compares them, saves the best one, and creates the daily forecast path plot.

---

## Train multiple tickers

```bash
python training/train_models.py --tickers SPY QQQ AAPL MSFT
```

---

## Use custom dates

```bash
python training/train_models.py \
  --tickers SPY \
  --download-start 2010-01-01 \
  --download-end 2026-01-01 \
  --train-start 2012-01-01 \
  --train-end 2018-12-31 \
  --val-end 2019-12-31 \
  --test-end 2020-12-31
```

For Windows CMD, use `^` instead of `\`.

---

## Important date logic

If the maximum forecast horizon is 1260 trading days, the model needs future data to create training labels.

Example:

```text
If train_end = 2018-12-31
and max horizon = 1260 trading days,
then the data should extend to around 2024.
```

This is necessary because the model must know the real price 1260 trading days after each training date during supervised learning.

---

## Outputs

After training SPY, you get:

```text
saved_models/SPY_best_model.pth
outputs/SPY/metrics.json
outputs/SPY/daily_forecast_path_from_training_end.csv
outputs/SPY/backtest_predictions_milestones.csv
outputs/SPY/plots/SPY_daily_forecast_path_from_training_end.png
```

The main plot is:

```text
outputs/SPY/plots/SPY_daily_forecast_path_from_training_end.png
```

It shows:

- actual stock price
- training period shading
- last training date
- model daily predicted future path
- confidence band
- markers for 1d, 1m, 1y, and 5y

---

## Run prediction-only GUI

```bash
python app.py
```

The GUI scans:

```text
saved_models/
```

and shows only tickers that already have saved trained models.

---

## Backtest saved model again

```bash
python backtesting/backtest_model.py --model-path saved_models/SPY_best_model.pth
```

This reloads the saved model and recreates the daily forecast path plot.

---


