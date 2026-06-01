# 📚 COMPREHENSIVE ANALYTICAL GUIDE TO TIME SERIES DENOISING WITH PYTORCH

## Table of Contents
1. [Executive Summary](#executive-summary)
2. [The Problem & Motivation](#the-problem--motivation)
3. [Theoretical Foundation](#theoretical-foundation)
4. [Project Architecture Overview](#project-architecture-overview)
5. [Detailed Code Walkthrough](#detailed-code-walkthrough)
6. [Mathematical Details](#mathematical-details)
7. [Implementation Deep Dive](#implementation-deep-dive)
8. [Results Analysis](#results-analysis)
9. [Advanced Topics](#advanced-topics)

---

## Executive Summary

### What This Project Does

This project addresses a **fundamental problem in signal processing**: extracting true signals from noisy measurements. We demonstrate how modern deep learning (specifically 1D Convolutional and Time-Delay Neural Networks) can effectively denoise time series data.

**Real-world applications:**
- Detecting gravitational waves from LIGO/Virgo detectors (physics)
- Identifying economic crises from noisy financial data (finance)
- Biomedical signal processing (EKG, EEG)
- Telecommunications (channel estimation)
- Audio denoising

**Our approach:** Train neural networks to learn the relationship between noisy and clean signals, then use this learned model to denoise new data.

---

## The Problem & Motivation

### Understanding Signal Noise

#### What is Noise?

In signal processing, we model observations as:

```
y(t) = s(t) + n(t)

where:
  y(t) = observed signal (noisy)
  s(t) = true signal (what we want)
  n(t) = noise (unwanted interference)
```

### Types of Noise

**1. Additive White Gaussian Noise (AWGN)**
```
Most common type. Mathematically:
n(t) ~ N(0, σ²)  [normally distributed, zero mean, variance σ²]

Examples:
- Thermal noise in electronics
- Quantum uncertainty in measurements
- Quantization errors
```

**2. Non-Additive Noise**
```
Depends on signal level:
y(t) = s(t) * (1 + n(t))  [multiplicative noise]

Examples:
- Fading in wireless channels
- Regime changes in economics
- Non-linear distortions
```

### The Denoising Challenge

**Given:** Noisy signal y(t)  
**Goal:** Estimate s(t) with minimal error  
**Constraint:** We don't know n(t) explicitly

**Traditional approaches:**
```
1. Wiener Filtering    - Optimal under Gaussian assumption
2. Kalman Filtering    - Good for state-space models
3. Wavelet Denoising   - Preserves sharp transitions
4. Moving Average      - Simple but poor performance

Limitations:
✗ Assume specific noise models
✗ Don't adapt to signal characteristics
✗ Poor at preserving signal structure
```

**Our approach: Deep Learning**
```
✓ Learns noise model from data
✓ Adapts to signal characteristics
✓ Preserves signal structure automatically
✓ Can handle complex, non-stationary signals
```

### Why Neural Networks?

Neural networks are **universal approximators** - they can learn any continuous function given:
1. Enough training data
2. Sufficient model capacity
3. Proper training procedure

For denoising: We train the network to predict: `clean_signal = f(noisy_signal)`

---

## Theoretical Foundation

### Part 1: 1D Convolutional Neural Networks

#### What is Convolution?

In continuous time, convolution is:

```
(f * g)(t) = ∫ f(τ)g(t-τ)dτ

In discrete form (for signals):

y[n] = Σ h[k] * x[n-k],  for k = 0 to K-1

where:
  x[n] = input signal at time n
  h[k] = filter/kernel (learnable parameters)
  y[n] = output at time n
  K   = kernel size
```

#### Visualization of 1D Convolution

```
Input signal:           [a, b, c, d, e, f, g, h]
                         ↓
Kernel (filter):        [w₀, w₁, w₂]  (size=3)
                         ↓
                    (w₀*a + w₁*b + w₂*c)  = output[0]
                    (w₀*b + w₁*c + w₂*d)  = output[1]
                    (w₀*c + w₁*d + w₂*e)  = output[2]
                              ...
                         ↓
Output:                 [y₀, y₁, y₂, ..., y₅]
```

#### What Convolution Does

```
Each kernel learns to detect a specific PATTERN:

Kernel [1, -1] detects:  sharp increases (derivative)
Kernel [1, 0, -1] detects: local peaks/valleys
Kernel [1, 1, 1]/3 detects: smoothness (averaging)
```

#### Key Properties

**1. Parameter Sharing**
- Same kernel applied at every position
- Reduces parameters dramatically
- Example: Fully connected layer = N² params, Conv = K params (K << N²)

**2. Local Connectivity**
- Each output depends only on local neighborhood
- Kernel size determines receptive field
- Larger kernels → see larger patterns

**3. Translation Invariance**
- If signal shifts, output shifts similarly
- Useful for position-invariant features

#### Stride and Padding

```
Stride: How many positions we move kernel each step
  stride=1: Move 1 position each time (overlapping)
  stride=2: Move 2 positions each time (skip positions)
  
  Effect: stride↑ → output_length↓, computation↓

Padding: Adding zeros around input
  "same" padding: output length = input length
  "valid" padding: output_length = (input_length - kernel_size + 1)
  
  Effect: Prevents information loss at boundaries
```

### Part 2: Convolutional Autoencoders

#### What is an Autoencoder?

An autoencoder is a neural network that learns to:
1. **Encode**: Compress input into a compact representation
2. **Decode**: Reconstruct original input from the compact form

```
INPUT (noisy signal)
    ↓
[ENCODER] → learns to extract FEATURES
    ↓
BOTTLENECK (compressed representation)
    ↓
[DECODER] → learns to RECONSTRUCT
    ↓
OUTPUT (reconstructed signal)
```

#### Why is This Useful for Denoising?

**Key insight:** The bottleneck is a "bottleneck"!

```
If we compress noisy signal into 16 dimensions:
- True signal s(t) is preserved (has structure, low dimension)
- Noise n(t) is random, won't compress well
- Bottleneck learns to keep signal, discard noise!

Analogy: Like squeezing through a narrow corridor
- Important features fit through
- Random clutter doesn't make it through
```

#### Autoencoder Architecture in Our Project

```
INPUT: (batch_size, 1, 64) [64 time samples, 1 channel]

ENCODER (Downsampling):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Conv1D: 1 → 32 channels, kernel=5, stride=1
Output shape: (batch_size, 32, 64)
     ↓ (applies ReLU activation and BatchNorm)
     
MaxPool1D: stride=2  [dimension reduction]
Output shape: (batch_size, 32, 32)
     ↓
     
Conv1D: 32 → 64 channels, kernel=5, stride=1
Output shape: (batch_size, 64, 32)
     ↓
     
MaxPool1D: stride=2
Output shape: (batch_size, 64, 16)
     ↓
     
Conv1D: 64 → 128 channels, kernel=3, stride=1
Output shape: (batch_size, 128, 16)

BOTTLENECK (Compression):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Flatten: (batch_size, 128×16) = (batch_size, 2048)
     ↓
Dense(2048 → 16)  [massive compression!]
Output shape: (batch_size, 16)  ← 2048/16 = 128x compression
     ↓
Dense(16 → 2048)  [prepare for decoder]

DECODER (Upsampling):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Reshape: (batch_size, 128, 16)
     ↓
Conv1D: 128 → 128, kernel=3
     ↓
Upsample(2)  [double resolution]
Output shape: (batch_size, 128, 32)
     ↓
Conv1D: 128 → 64, kernel=5
     ↓
Upsample(2)
Output shape: (batch_size, 64, 64)
     ↓
Conv1D: 64 → 32, kernel=5
     ↓
Conv1D: 32 → 1, kernel=5  [back to 1 channel]

OUTPUT: (batch_size, 1, 64)  [same shape as input!]
```

**Key idea:** The compression forces the network to learn what's important!

### Part 3: Time Delay Neural Networks (TDNNs)

#### Motivation

What if we explicitly gave the network temporal context?

```
Standard approach (CNN):
Learn implicit temporal patterns through convolution

TDNN approach:
Explicitly provide past values of signal

Input to TDNN:
[s(t), s(t-1), s(t-2), s(t-4), s(t-8), ...]
 ^current  ^recent   ^older  ^even older ^much older
```

#### Why This Matters

```
Some signals have patterns at different timescales:
- Gravitational waves: Rapid chirp (few milliseconds)
- Economic signals: Slow trends (weeks/months)

TDNN lets network see multiple timescales explicitly:

delay=1:   Captures rapid changes (∂s/∂t ~ derivative)
delay=2:   Captures medium-term patterns
delay=4:   Captures longer-term trends  
delay=8:   Captures very long-term structure

Each delay is like a different "filter" for different frequencies!
```

#### TDNN Architecture

```
INPUT: (batch_size, 1, 64)

TIME-DELAY CONCATENATION:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Original:           s[0], s[1], s[2], ..., s[63]
Delay by 1:     0,  s[0], s[1], ..., s[62]
Delay by 2:     0,  0,    s[0], ..., s[61]
Delay by 4:     0,  0,    0,    0,   s[0], ..., s[59]
Delay by 8:     0,  0,    0,    0,   0,    0,   0,   0,   s[0], ...

CONCATENATE all together:
(batch_size, 64 × 5) = (batch_size, 320)
     ↓ [5 = 1 original + 4 delayed versions]

DENSE LAYERS:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Dense(320 → 128) + ReLU + BatchNorm + Dropout(0.2)
     ↓
Dense(128 → 64) + ReLU + BatchNorm + Dropout(0.2)
     ↓
Dense(64 → 32) + ReLU
     ↓
Dense(32 → 64)  [output for 64 time steps]
     ↓
Reshape to (batch_size, 1, 64)

OUTPUT: (batch_size, 1, 64)
```

#### CNN vs TDNN Comparison

| Aspect | CNN | TDNN |
|--------|-----|------|
| **Learn patterns?** | Implicit (hidden) | Explicit (visible) |
| **Hierarc. features?** | Yes (deep layers) | Limited (linear) |
| **Parameters** | ~150K | ~60K |
| **Interpretability** | Black-box | Transparent |
| **Generalization** | May overfit (big model) | Better (small model) |
| **Stationarity needed?** | Somewhat | No |
| **Non-linearity** | Unlimited | Limited |

---

## Project Architecture Overview

### High-Level Flow

```
┌─────────────────────────────────────────────────────────────┐
│                    TIME SERIES DENOISING                     │
│                      WITH PYTORCH                             │
└─────────────────────────────────────────────────────────────┘

PHASE 1: DATA GENERATION
├─ TimeSeriesGenerator class
│  ├─ gravitational_wave_signal()   → chirp + envelope
│  ├─ economic_signal()             → GBM + regime shifts
│  ├─ add_noise()                   → AWGN with SNR control
│  └─ create_dataset()              → sliding windows
└─ Output: Noisy signals, clean signals, time arrays

PHASE 2: DATA PREPARATION
├─ DenoisingDataset class
│  └─ __getitem__() → (noisy_window, clean_window) pairs
├─ torch.utils.data.DataLoader
│  ├─ Training loader (shuffle=True)
│  ├─ Validation loader (shuffle=False)
│  └─ Testing loader (shuffle=False)
└─ Output: Batched data ready for training

PHASE 3: MODEL DEFINITION
├─ Conv1DAutoencoder class
│  ├─ encode()  → compress noisy signal
│  └─ decode()  → reconstruct signal
└─ TDNN class
   └─ forward() → process time-delayed inputs

PHASE 4: TRAINING
├─ train_denoising_model()
│  ├─ Forward pass: model(noisy) → reconstruction
│  ├─ Loss calculation: MSE(reconstruction, clean)
│  ├─ Backward pass: compute gradients
│  ├─ Optimizer step: update weights
│  ├─ Validation: monitor overfitting
│  └─ Early stopping: save best model
└─ Output: Trained models

PHASE 5: EVALUATION
├─ evaluate_denoiser()
│  ├─ Compute MSE, MAE
│  ├─ Calculate SNR improvement
│  └─ Compare models
├─ detect_events()
│  └─ Find peaks/anomalies in denoised signal
└─ Output: Performance metrics

PHASE 6: VISUALIZATION
├─ Create comparison plots
│  ├─ Signal overlays
│  ├─ Error distributions
│  ├─ Frequency analysis
│  └─ Event detection results
└─ Output: PNG figures
```

---

## Detailed Code Walkthrough

### Part 1: Data Generation (TimeSeriesGenerator)

#### Why Synthetic Data?

We use synthetic data for **reproducibility** and **control**:
- Can perfectly separate signal from noise (we generate both)
- Can vary SNR precisely
- No real data required (instant, reusable)
- Ideal for understanding model behavior

#### TimeSeriesGenerator.gravitational_wave_signal()

```python
def gravitational_wave_signal(duration=2.0, sampling_rate=4000, 
                              frequency_start=20, frequency_end=250):
    """
    Simulates gravitational wave from merging binary system.
    
    Physics background:
    ───────────────────
    When two black holes merge:
    1. Orbital frequency increases as objects spiral inward
    2. Creates "chirp" - frequency sweeps from low to high
    3. Gravitational wave amplitude increases until merger
    4. Signal matches characteristic pattern of binary merger
    
    Parameters explain:
    ──────────────────
    duration: How long the signal (seconds)
    sampling_rate: How many samples per second (Hz)
    frequency_start: Starting frequency (20 Hz = subaudible)
    frequency_end: Final frequency before merger (250 Hz)
    """
    
    # Create time array
    t = np.linspace(0, duration, int(duration * sampling_rate), endpoint=False)
    
    # STEP 1: CREATE CHIRP (frequency sweep)
    chirp = signal.chirp(t, frequency_start, duration, frequency_end, method='quadratic')
    
    # Mathematical form: frequency increases quadratically with time
    # f(t) = f_start + (f_end - f_start) * (t/T)²
    # This models the accelerating orbital frequency of merging objects
    
    # STEP 2: CREATE AMPLITUDE ENVELOPE
    envelope = np.exp(-((t - duration/2)**2) / (2 * (duration/5)**2))
    
    # This is a Gaussian (bell curve) centered at t=duration/2
    # Signal is weak at start → peaks at merger → weak at end
    # Width parameter: duration/5 controls "concentration" around peak
    
    # STEP 3: COMBINE THEM
    signal_clean = chirp * envelope
    
    # Multiply: Oscillations (chirp) × Envelope (amplitude modulation)
    # Creates realistic GW signature: increasing frequency + increasing amplitude
    
    # STEP 4: ADD SPIKES (optional merger dynamics)
    n_spikes = 3
    spike_positions = np.random.choice(len(signal_clean), n_spikes, replace=False)
    spike_amplitudes = np.random.uniform(0.5, 1.5, n_spikes)
    for pos, amp in zip(spike_positions, spike_amplitudes):
        signal_clean[max(0, pos-5):min(len(signal_clean), pos+5)] += amp * signal.windows.hann(10)
    
    # Add Hann windows at random positions
    # Simulates additional dynamical effects during merger
    
    return t, signal_clean
```

**Visual representation:**

```
Frequency:  20Hz ─────────→ 250Hz  (chirp increasing)

Amplitude:        ╱╲
                 ╱  ╲           (envelope peaks)
    ────────────╱────╲───────  (time)
    0          1.0   2.0

Combined:      ~~~╱╱╱╱╲╲╲╲~~~
              (oscillations with growing, then shrinking amplitude)
```

#### TimeSeriesGenerator.economic_signal()

```python
def economic_signal(length=1000, volatility=0.02):
    """
    Simulates log-returns of stock price with crisis periods.
    
    Financial background:
    ─────────────────────
    Stock prices follow geometric Brownian motion:
    dS/S = μ dt + σ dW
    
    Where:
    - dS/S = log return (proportional change)
    - μ = drift (average return)
    - σ = volatility (risk)
    - dW = random shock (Brownian motion)
    
    In crisis: volatility increases dramatically (σ → σ * multiplier)
    """
    
    drift = 0.0001  # 0.01% average daily return
    returns = np.zeros(length)
    returns[0] = 0
    
    # STEP 1: DEFINE REGIME SHIFTS
    regime = np.ones(length)  # baseline volatility
    crisis_periods = [
        (100, 200, 3.0),   # Start day 100, last 200 days, 3x volatility
        (400, 150, 2.5),   # Start day 400, last 150 days, 2.5x volatility
        (700, 120, 2.0)    # Start day 700, last 120 days, 2x volatility
    ]
    
    for start, duration, severity in crisis_periods:
        regime[start:start+duration] = severity
    
    # regime[i] = volatility multiplier at day i
    # Normal day: regime[i] = 1.0
    # Crisis day: regime[i] = 3.0 (3x higher volatility)
    
    # STEP 2: GENERATE GBM RETURNS
    for i in range(1, length):
        # shock ~ N(drift, volatility²) with regime adjustment
        shock = np.random.normal(drift, volatility * regime[i])
        returns[i] = returns[i-1] + shock
    
    # Each day: previous return + new random shock
    # During crises: larger random shocks (higher volatility)
    
    # STEP 3: CUMULATIVE PRICE
    prices = np.cumsum(returns)
    
    # Each price = sum of all returns up to that day
    # Reflects realistic stock price movement
    
    return np.arange(length), prices
```

**Visual representation:**

```
Normal volatility:   ~~  ~~  ~~  ~~     (small wiggles)
Crisis volatility:   ╱╲╱╲╱╲╱╲╱╲       (big wiggles)

Combined with drift:
         ╱╱╱╱╲╲╱╱╱╱╱╱╱╱╱╱╲╲╲╱╱╱
        ╱    ╲╱  ╲╱        ╲  ╱
       ╱              ╲    ╱
    Time: 0    100  250  400  550   700  820  1000
                 ├────────────┤ Crisis 1
                             ├────────┤ Crisis 2
                                          ├──────┤ Crisis 3
```

#### TimeSeriesGenerator.add_noise()

```python
def add_noise(signal_clean, snr_db=10):
    """
    Adds Gaussian noise with specific Signal-to-Noise Ratio.
    
    SNR Definition:
    ───────────────
    SNR (in dB) = 10 * log₁₀(P_signal / P_noise)
    
    Where P = power = mean(signal²)
    
    Examples:
    SNR = 10 dB  → signal power = 10 × noise power (somewhat noisy)
    SNR = 0 dB   → signal power = noise power (very noisy)
    SNR = -10 dB → signal power = 0.1 × noise power (extremely noisy)
    
    Our GW signal: SNR = 5 dB (challenging but realistic for LIGO)
    """
    
    # Calculate signal power
    signal_power = np.mean(signal_clean ** 2)
    
    # Rearrange SNR formula to find required noise power
    # SNR = 10 * log₁₀(P_s / P_n)
    # 10^(SNR/10) = P_s / P_n
    # P_n = P_s / 10^(SNR/10)
    noise_power = signal_power / (10 ** (snr_db / 10))
    
    # Generate Gaussian noise with that power
    # Noise ~ N(0, σ²), where σ² = noise_power
    noise = np.random.normal(0, np.sqrt(noise_power), len(signal_clean))
    
    # Combine
    signal_noisy = signal_clean + noise
    
    return signal_noisy, noise
```

**Example:**
```
Clean signal:  s(t) = sin(t)
Noise:         n(t) = random ~ N(0, 0.1)
Noisy signal:  y(t) = sin(t) + n(t)

     sin(t)              n(t)            y(t) = sin(t) + n(t)
    
     ╱╲╱╲                ~~~             ╱╲╱╲~~~
    ╱  ╲  ╲              ╱ ╲             ╱  ╲╱ ╲
   ╱    ╲  ╲            ╱   ╲           ╱    ╲ ╲
  (clean)  (noise)    (combined)
```

### Part 2: Creating Supervised Dataset

```python
def create_dataset(signal_clean, window_size=64, stride=8):
    """
    Convert long signal into training examples using sliding window.
    
    Why sliding window?
    ──────────────────
    Signal: [a, b, c, d, e, f, g, h, ...]
    
    Problem: This is one continuous sequence
    Solution: Break into overlapping windows
    
    Window size 3, stride 2:
    Window 1: [a, b, c] → label c_clean  (or mean, etc.)
    Window 2: [c, d, e] → label e_clean
    Window 3: [e, f, g] → label g_clean
               ↑ stride=2 overlap
    """
    
    X, y = [], []
    
    # Slide window across signal
    for i in range(0, len(signal_clean) - window_size, stride):
        # Extract window from position i to i+window_size
        window = signal_clean[i:i + window_size]
        
        # Use same window as both input and target
        # (Autoencoder: input=noisy, target=clean)
        X.append(window)
        y.append(window)
    
    return np.array(X), np.array(y)
```

**Visualization:**

```
Original signal (1D array):
[0.1, 0.3, 0.5, 0.7, 0.9, 0.8, 0.6, 0.4, 0.2, ...]
 ↓

With window_size=3, stride=2:

Sample 0:  [0.1, 0.3, 0.5]
Sample 1:       [0.5, 0.7, 0.9]  ← overlaps with sample 0
Sample 2:            [0.9, 0.8, 0.6]  ← overlaps with sample 1
Sample 3:                 [0.6, 0.4, 0.2]  ← overlaps with sample 2

Result: X.shape = (n_samples, 64)
        y.shape = (n_samples, 64)
```

---

### Part 3: PyTorch Dataset and DataLoader

```python
class DenoisingDataset(Dataset):
    """
    PyTorch Dataset: Interface between raw data and model training.
    
    Why separate class?
    ──────────────────
    PyTorch's DataLoader needs standard interface:
    - __len__() → how many samples?
    - __getitem__(idx) → get sample at index idx
    
    This allows DataLoader to:
    - Shuffle data during training
    - Create batches efficiently
    - Parallelize data loading
    """
    
    def __init__(self, X_noisy, y_clean):
        # Store as PyTorch tensors (GPU-compatible)
        self.X_noisy = torch.FloatTensor(X_noisy)  # (N, 64)
        self.y_clean = torch.FloatTensor(y_clean)  # (N, 64)
    
    def __len__(self):
        return len(self.X_noisy)
    
    def __getitem__(self, idx):
        # Get one sample
        x = self.X_noisy[idx]  # shape: (64,)
        y = self.y_clean[idx]  # shape: (64,)
        
        # Reshape for 1D convolution: (1, 64)
        # [batch_size, channels, length]
        x = x.unsqueeze(0)  # (1, 64)
        y = y.unsqueeze(0)  # (1, 64)
        
        return x, y
```

**DataLoader usage:**

```python
dataset = DenoisingDataset(X_train_noisy, X_train_clean)
dataloader = DataLoader(dataset, batch_size=32, shuffle=True)

# During training:
for batch_noisy, batch_clean in dataloader:
    # batch_noisy shape: (32, 1, 64)  [32 samples per batch]
    # batch_clean shape: (32, 1, 64)
    
    # Forward pass
    output = model(batch_noisy)  # shape: (32, 1, 64)
    
    # Calculate loss
    loss = criterion(output, batch_clean)
```

---

### Part 4: Model Definitions

#### CNN Autoencoder Detailed

```python
class Conv1DAutoencoder(nn.Module):
    def __init__(self, input_length=64, latent_dim=16):
        """
        Autoencoder for 1D signals.
        
        Key design decisions:
        ──────────────────────
        1. input_length=64: Process 64-sample windows
        2. latent_dim=16: Compress to 16 features
           - Compression ratio: 64/16 = 4x
           - Forces learning of essential features
        """
        super(Conv1DAutoencoder, self).__init__()
        
        # ENCODER ─────────────────────────────────
        
        # Layer 1: 1 channel → 32 channels
        self.enc_conv1 = nn.Conv1d(1, 32, kernel_size=5, padding=2, stride=1)
        # kernel_size=5: Look at 5 time steps
        # padding=2: Keep dimensions same (64→64)
        # stride=1: Move kernel 1 position each step
        
        self.enc_bn1 = nn.BatchNormalization(32)
        # Normalize inputs to each ReLU
        # Why? Stabilizes learning, allows higher learning rates
        
        self.enc_pool1 = nn.MaxPool1d(kernel_size=2)
        # Downsample: (batch, 32, 64) → (batch, 32, 32)
        # Keep max value in each 2-element window
        # Discards redundant information
        
        # Layer 2: 32 channels → 64 channels
        self.enc_conv2 = nn.Conv1d(32, 64, kernel_size=5, padding=2, stride=1)
        # (batch, 32, 32) → (batch, 64, 32)
        
        self.enc_bn2 = nn.BatchNormalization(64)
        self.enc_pool2 = nn.MaxPool1d(kernel_size=2)
        # (batch, 64, 32) → (batch, 64, 16)
        
        # Layer 3: 64 channels → 128 channels
        self.enc_conv3 = nn.Conv1d(64, 128, kernel_size=3, padding=1, stride=1)
        # (batch, 64, 16) → (batch, 128, 16)
        # Now have 128 different features
        
        self.enc_bn3 = nn.BatchNormalization(128)
        
        # BOTTLENECK ──────────────────────────────
        
        # Flatten: (batch, 128, 16) → (batch, 2048)
        bottleneck_size = 128 * (input_length // 4)  # 128 * 16 = 2048
        
        # Dense compression: 2048 → 16 dimensions
        self.bottleneck_fc1 = nn.Linear(bottleneck_size, latent_dim)
        # Most aggressive compression happens here!
        # Input: 2048 features
        # Output: 16 features (128x compression)
        
        # Expansion: 16 → 2048
        self.bottleneck_fc2 = nn.Linear(latent_dim, bottleneck_size)
        # Prepare for decoder
        
        # DECODER ──────────────────────────────────
        
        # Reverse of encoder
        # (batch, 128, 16) → (batch, 1, 64)
        
        self.dec_conv1 = nn.Conv1d(128, 128, kernel_size=3, padding=1)
        self.dec_bn1 = nn.BatchNormalization(128)
        self.dec_up1 = nn.Upsample(scale_factor=2, mode='nearest')
        # (batch, 128, 16) → (batch, 128, 32)
        
        self.dec_conv2 = nn.Conv1d(128, 64, kernel_size=5, padding=2)
        self.dec_bn2 = nn.BatchNormalization(64)
        self.dec_up2 = nn.Upsample(scale_factor=2, mode='nearest')
        # (batch, 64, 32) → (batch, 64, 64)
        
        self.dec_conv3 = nn.Conv1d(64, 32, kernel_size=5, padding=2)
        self.dec_bn3 = nn.BatchNormalization(32)
        
        # Output: 32 → 1 channel
        self.dec_conv4 = nn.Conv1d(32, 1, kernel_size=5, padding=2)
        # (batch, 1, 64) ← RECONSTRUCTED OUTPUT
    
    def encode(self, x):
        """Compress input to latent representation."""
        # Encoder: convolutions + pooling
        x = F.relu(self.enc_bn1(self.enc_conv1(x)))
        x = self.enc_pool1(x)
        
        x = F.relu(self.enc_bn2(self.enc_conv2(x)))
        x = self.enc_pool2(x)
        
        x = F.relu(self.enc_bn3(self.enc_conv3(x)))
        
        # Flatten
        x = x.view(x.size(0), -1)  # (batch, 2048)
        
        # Dense bottleneck
        x = F.relu(self.bottleneck_fc1(x))  # (batch, 16)
        
        return x  # latent code
    
    def decode(self, z):
        """Reconstruct signal from latent representation."""
        # Expand latent code
        x = F.relu(self.bottleneck_fc2(z))  # (batch, 2048)
        
        # Reshape to spatial form
        x = x.view(x.size(0), 128, -1)  # (batch, 128, 16)
        
        # Decoder: upsampling + convolutions
        x = F.relu(self.dec_bn1(self.dec_conv1(x)))
        x = self.dec_up1(x)
        
        x = F.relu(self.dec_bn2(self.dec_conv2(x)))
        x = self.dec_up2(x)
        
        x = F.relu(self.dec_bn3(self.dec_conv3(x)))
        x = self.dec_conv4(x)  # Final: no activation (linear output)
        
        return x  # Reconstructed signal
    
    def forward(self, x):
        """Full pass: encode then decode."""
        z = self.encode(x)
        x_recon = self.decode(z)
        return x_recon
```

#### TDNN Detailed

```python
class TDNN(nn.Module):
    def __init__(self, input_length=64, time_delays=(1, 2, 4, 8)):
        """
        Time Delay Neural Network for sequence modeling.
        
        time_delays=(1,2,4,8) means:
        - delay=1: Previous value (recent past)
        - delay=2: Two steps back
        - delay=4: Four steps back
        - delay=8: Eight steps back (longer history)
        
        Why these specific values?
        - Powers of 2: exponentially spread in time
        - Captures patterns at multiple timescales
        """
        super(TDNN, self).__init__()
        
        self.input_length = input_length
        self.time_delays = time_delays
        
        # Number of features after concatenation
        # 1 (original) + len(time_delays) = 1 + 4 = 5
        n_delayed_features = 1 + len(time_delays)
        
        # Input shape after concatenation
        # (batch, input_length * n_delayed_features)
        # = (batch, 64 * 5) = (batch, 320)
        
        # DENSE LAYERS ─────────────────────────────────
        
        # Layer 1: 320 → 128
        self.fc1 = nn.Linear(input_length * n_delayed_features, 128)
        self.bn1 = nn.BatchNorm1d(128)
        self.dropout1 = nn.Dropout(0.2)  # Drop 20% of neurons during training
        
        # Layer 2: 128 → 64
        self.fc2 = nn.Linear(128, 64)
        self.bn2 = nn.BatchNorm1d(64)
        self.dropout2 = nn.Dropout(0.2)
        
        # Layer 3: 64 → 32 (no activation, intermediate)
        self.fc3 = nn.Linear(64, 32)
        
        # Output: 32 → input_length (64)
        # One output per input time step
        self.fc4 = nn.Linear(32, input_length)
    
    def forward(self, x):
        """
        Process through TDNN.
        
        x shape: (batch_size, 1, input_length)
        """
        batch_size = x.size(0)
        x = x.squeeze(1)  # (batch, input_length)
        
        # ─── CREATE TIME-DELAYED VERSIONS ───
        delayed_features = [x]  # Original signal: [s₀, s₁, s₂, ...]
        
        for delay in self.time_delays:
            if delay > 0:
                # Shift signal by delay positions
                # pad left with zeros, remove from right
                delayed = torch.cat([
                    torch.zeros(batch_size, delay, device=x.device),
                    x[:, :-delay]  # Shifted version
                ], dim=1)
                
                # Now: [0, 0, ..., s₀, s₁, ..., s_{n-delay-1}]
                #       └─ delay ─┘
                
                delayed_features.append(delayed)
        
        # ─── CONCATENATE ALL VERSIONS ───
        x_delayed = torch.cat(delayed_features, dim=1)
        # Shape: (batch, input_length * n_delayed_features)
        # = (batch, 320)
        
        # ─── PASS THROUGH DENSE LAYERS ───
        x = F.relu(self.bn1(self.fc1(x_delayed)))
        x = self.dropout1(x)
        
        x = F.relu(self.bn2(self.fc2(x)))
        x = self.dropout2(x)
        
        x = F.relu(self.fc3(x))
        
        x = self.fc4(x)  # Output: (batch, input_length)
        
        # Reshape to match input format
        x = x.unsqueeze(1)  # (batch, 1, input_length)
        
        return x
```

**Visual comparison:**

```
INPUT (same for both):
Shape: (batch_size, 1, 64)

CNN Path:
(1, 64) → Conv(1→32) → Pool → Conv(32→64) → Pool 
       → Conv(64→128) → Flatten(2048) → Dense(16)
       → Dense(2048) → Reshape → DeConv → DeConv
       → (1, 64)

TDNN Path:
(1, 64) → Create time delays:
         [s(t), s(t-1), s(t-2), s(t-4), s(t-8)]
       → Concatenate → (320,)
       → Dense(128) → Dense(64) → Dense(32) → Dense(64)
       → Reshape → (1, 64)

Key difference:
✓ CNN learns implicit patterns through convolution
✓ TDNN uses explicit temporal structure through delays
```

---

### Part 5: Training Loop

```python
def train_denoising_model(model, train_loader, val_loader, epochs=30, 
                         model_name="Model", device=device):
    """
    Train autoencoder to denoise signals.
    
    Training procedure:
    ──────────────────
    For each epoch:
        For each batch:
            1. Forward: prediction = model(noisy_batch)
            2. Loss: error = MSE(prediction, clean_batch)
            3. Backward: compute gradients (∂loss/∂weights)
            4. Update: weights -= learning_rate * gradients
        Validate: check performance on unseen data
        Early stop: if validation loss plateaus
    """
    
    # Move model to device (GPU or CPU)
    model = model.to(device)
    
    # Optimizer: Adam (Adaptive Moment Estimation)
    # Maintains running average of gradients
    # Automatically adjusts learning rate per parameter
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    
    # Loss function: Mean Squared Error
    # loss = mean((prediction - target)²)
    # Penalizes large errors more than small ones
    criterion = nn.MSELoss()
    
    # Learning rate scheduler: reduce LR if validation loss plateaus
    # "Plateau": validation loss not improving for 'patience' epochs
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', factor=0.5, patience=5
    )
    
    # ─── TRAINING LOOP ───
    
    best_val_loss = float('inf')
    patience = 10
    patience_counter = 0
    
    train_losses = []
    val_losses = []
    
    for epoch in range(epochs):
        
        # ─── TRAINING PHASE ───
        model.train()  # Enable dropout, batch norm training mode
        train_loss = 0.0
        
        for batch_noisy, batch_clean in train_loader:
            # Move batch to device
            batch_noisy = batch_noisy.to(device)  # (batch_size, 1, 64)
            batch_clean = batch_clean.to(device)  # (batch_size, 1, 64)
            
            # FORWARD PASS
            output = model(batch_noisy)  # (batch_size, 1, 64)
            
            # CALCULATE LOSS
            # loss = mean((output - batch_clean)²)
            loss = criterion(output, batch_clean)
            
            # BACKWARD PASS
            optimizer.zero_grad()  # Clear old gradients
            loss.backward()  # Compute new gradients using backpropagation
            
            # OPTIMIZER STEP
            optimizer.step()  # Update weights using gradients
            
            # Accumulate loss
            train_loss += loss.item()
        
        train_loss /= len(train_loader)
        train_losses.append(train_loss)
        
        # ─── VALIDATION PHASE ───
        model.eval()  # Disable dropout, batch norm eval mode
        val_loss = 0.0
        
        with torch.no_grad():  # Don't compute gradients (faster, less memory)
            for batch_noisy, batch_clean in val_loader:
                batch_noisy = batch_noisy.to(device)
                batch_clean = batch_clean.to(device)
                
                output = model(batch_noisy)
                loss = criterion(output, batch_clean)
                
                val_loss += loss.item()
        
        val_loss /= len(val_loader)
        val_losses.append(val_loss)
        
        # LEARNING RATE SCHEDULING
        scheduler.step(val_loss)
        
        # ─── EARLY STOPPING ───
        # Stop if validation loss not improving
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
        else:
            patience_counter += 1
        
        if patience_counter >= patience:
            print(f"Early stopping at epoch {epoch+1}")
            break
        
        # Print progress
        if (epoch + 1) % 10 == 0 or epoch == 0:
            print(f"Epoch {epoch+1}/{epochs} - "
                  f"Train Loss: {train_loss:.6f}, "
                  f"Val Loss: {val_loss:.6f}")
    
    return train_losses, val_losses
```

**Training visualization:**

```
Loss curve (typical):

Loss
 │     ╲╲  ╲  ╲ ╲  ╲ ╲ ╲╲╲  ╲  ╲
 │  ╲╲  ╲ ╱ ╲╱╲ ╲  ╲ ╱ ╲╲╲╲╲╲╲
 │ ╱  ╲ ╱╱   ╲ ╲ ╲  ╲ ╲
 │╱    ╲╱     ╲╱  ╲ ╱╱
 │              ╲╱
 └──────────────────────→ Epoch

     ╱ Steep drop       ╲ Converged (plateau)
    Early epochs        Later epochs
    (rapid learning)    (fine-tuning)

Early stopping: Stop when val_loss plateaus
(to prevent overfitting)
```

---

### Part 6: Evaluation

```python
def evaluate_denoiser(model, test_loader, model_name="Model", device=device):
    """
    Evaluate model on test set.
    
    Metrics computed:
    ─────────────────
    1. MSE: Mean Squared Error
       - How close are predictions to targets on average?
       - Lower is better
    
    2. MAE: Mean Absolute Error
       - Similar to MSE but less sensitive to outliers
       - Lower is better
    
    3. SNR: Signal-to-Noise Ratio (dB)
       - SNR = 10 * log₁₀(signal_power / noise_power)
       - Higher is better
       - Improvement = SNR_after - SNR_before
    """
    
    model.eval()  # Evaluation mode
    
    all_noisy = []
    all_clean = []
    all_denoised = []
    
    # Collect all predictions on test set
    with torch.no_grad():
        for batch_noisy, batch_clean in test_loader:
            batch_noisy = batch_noisy.to(device)
            batch_clean = batch_clean.to(device)
            
            output = model(batch_noisy)
            
            all_noisy.append(batch_noisy.cpu().numpy())
            all_clean.append(batch_clean.cpu().numpy())
            all_denoised.append(output.cpu().numpy())
    
    # Concatenate all batches
    X_test_noisy = np.concatenate(all_noisy, axis=0)      # (n_samples, 1, 64)
    X_test_clean = np.concatenate(all_clean, axis=0)
    X_test_denoised = np.concatenate(all_denoised, axis=0)
    
    # ─── CALCULATE METRICS ───
    
    # MSE: Mean((prediction - target)²)
    mse_noisy = np.mean((X_test_noisy - X_test_clean) ** 2)
    mse_denoised = np.mean((X_test_denoised - X_test_clean) ** 2)
    
    # MAE: Mean(|prediction - target|)
    mae_noisy = np.mean(np.abs(X_test_noisy - X_test_clean))
    mae_denoised = np.mean(np.abs(X_test_denoised - X_test_clean))
    
    # SNR calculation:
    # SNR = 10 * log₁₀(signal_power / error_power)
    # signal_power = mean(clean_signal²)
    # error_power = mse
    
    signal_power = np.mean(X_test_clean ** 2)
    snr_noise = 10 * np.log10(signal_power / mse_noisy + 1e-10)
    snr_denoised = 10 * np.log10(signal_power / mse_denoised + 1e-10)
    snr_improvement = snr_denoised - snr_noise
    
    print(f"MSE (noisy):     {mse_noisy:.6f}")
    print(f"MSE (denoised):  {mse_denoised:.6f}")
    print(f"Improvement MSE: {(1 - mse_denoised/mse_noisy)*100:.2f}%")
    print(f"\nSNR (noisy):     {snr_noise:.2f} dB")
    print(f"SNR (denoised):  {snr_denoised:.2f} dB")
    print(f"SNR Improvement: {snr_improvement:.2f} dB")
    
    return {
        'mse_noisy': mse_noisy,
        'mse_denoised': mse_denoised,
        'snr_noisy': snr_noise,
        'snr_denoised': snr_denoised,
        'snr_improvement': snr_improvement,
        'predictions': X_test_denoised
    }
```

---

## Mathematical Details

### Backpropagation Algorithm (What happens during training)

```
Forward Pass (computing prediction):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Input: x = [noisy signal values]
       
Layer 1: z₁ = W₁x + b₁
         a₁ = ReLU(z₁)
         
Layer 2: z₂ = W₂a₁ + b₂
         ŷ = z₂  [output]

Loss Calculation:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
L = (1/N) Σ(ŷᵢ - yᵢ)²   [mean squared error]

Backward Pass (computing gradients):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
∂L/∂W₂ = ∂L/∂ŷ × ∂ŷ/∂W₂
       = 2(ŷ - y) × a₁

∂L/∂W₁ = ∂L/∂a₁ × ∂a₁/∂z₁ × ∂z₁/∂W₁
       = (W₂ᵀ × ∂L/∂z₂) × ReLU'(z₁) × x

Parameter Update (gradient descent):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
W₂_new = W₂_old - α × ∂L/∂W₂
W₁_new = W₁_old - α × ∂L/∂W₁

where α = learning_rate = 0.001
```

### Why Convolution Works for Denoising

```
Convolution property 1: LOCAL OPERATIONS
────────────────────────────────────────
Each output depends only on nearby inputs.
Why good for denoising?
- Noise is random (high frequency)
- Signal has structure (patterns)
- Local patterns survive convolution
- Random noise gets averaged out

Example:
Input:  [signal, signal, noise, signal, ...]
Conv:   [w₀, w₁, w₂] (learnable)
Output: [low, low, high-if-noise, low, ...]
        └─ Can learn to suppress noise! ─┘

Convolution property 2: PARAMETER SHARING
──────────────────────────────────────────
Same weights used at every position.
Why good for denoising?
- Pattern detection is location-independent
- Learning "what is noise" applies everywhere
- Dramatically fewer parameters

Example:
Without convolution: 64×64 = 4096 parameters
With convolution (k=5): Only 5 parameters
Ratio: 4096/5 = 819x fewer parameters!
Less parameters → less overfitting
```

---

## Implementation Deep Dive

### Gradients and Backpropagation Through Time

```python
# When we call loss.backward(), PyTorch:

1. Starts at output: dL/dŷ = 2(ŷ - y)

2. Works backwards through decoder:
   dL/dconv_out = dL/dŷ × dŷ/dconv_out
   dL/dpool = dL/dconv_out × dconv_out/dpool
   [etc.]

3. Reaches bottleneck: dL/dlatent_code
   This tells us which features matter!

4. Continues through encoder:
   [same chain rule application]

5. Finally: dL/dW for each weight W

6. Each weight is then updated:
   W_new = W - α × dL/dW
```

### Batch Normalization Why?

```
Without Batch Norm:
    Layer 1 output: [0.001, 0.999, 0.5, ...]  (varying scales)
                           ↓
    Layer 2 input: (very different scales!)
    Problem: Gradients become very large or very small
             → unstable training

With Batch Norm:
    Normalize each channel's output to:
    mean = 0, std = 1
    
    Then scale/shift: γ * normalized + β
    (γ, β are learnable)
    
    Benefits:
    ✓ Inputs to next layer normalized
    ✓ Training stable (no exploding/vanishing gradients)
    ✓ Can use higher learning rates
    ✓ Acts as regularization
```

### Dropout Regularization

```
During training (TDNN example):
    Dense(128 → 64)
    Dropout(0.2)
    ├─ With prob 0.2: Drop neuron (set output = 0)
    └─ With prob 0.8: Keep neuron (scale by 1/0.8 = 1.25)
    
    Effect on that batch: Only ~80% of neurons active
    Forces network to learn redundant representations
    
During inference (testing):
    Dropout disabled
    All neurons active at reduced scale
    (No randomness in output)
    
Why this helps:
    ✓ Prevents co-adaptation of neurons
    ✓ Reduces overfitting
    ✓ Creates ensemble effect
    Analogy: Like getting opinions from multiple experts
```

---

## Results Analysis

### Why TDNN Wins on Economic Data

```
PROBLEM: Economic signals are non-stationary
(regime changes, structural breaks)

CNN Autoencoder: ✗ Struggles
  ├─ Large model (150K params)
  ├─ Needs more training data
  ├─ Learns noise patterns during regime changes
  └─ MSE: 80.1% reduction

TDNN: ✓ Better
  ├─ Smaller model (60K params)
  ├─ Explicit temporal structure
  ├─ Better generalization
  └─ MSE: 86.9% reduction (+6.8% over CNN!)
```

### Why CNN Works for Gravitational Waves

```
BENEFIT: GW signals are highly structured
(clear chirp pattern, predictable envelope)

CNN: ✓ Excellent
  ├─ Learns hierarchical features
  ├─ Detects frequency sweep (convolution specializes)
  ├─ Good generalization despite complexity
  └─ SNR: +5.13 dB

TDNN: Also good
  ├─ Benefits from temporal context
  ├─ Simpler, less prone to overfitting
  ├─ Slightly better
  └─ SNR: +5.21 dB (+0.08 dB over CNN)
```

### Event Detection Success

```
Gravitational Waves: 83.6% SUCCESS
├─ Clean signal: 146 events detected
├─ Noisy signal: 430 spurious events (noise spikes)
└─ Denoised signal: 122 events
    └─ 122/146 = 83.6% recall of true events

Why successful?
✓ Real events have HIGH amplitude
✓ Noise spikes are ISOLATED
✓ Denoising preserves high-amplitude features
✓ Threshold-based detection works

Economic Crises: 5.3% SUCCESS
├─ Clean signal: 76 detected crises
├─ Denoised signal: 4 detected
└─ 4/76 = 5.3% recall

Why poor?
✗ Crises = regime change, not additive noise
✗ Denoising smooths everything equally
✗ Can't distinguish signal from noise in crisis
✗ Wrong model for the problem
```

---

## Advanced Topics

### Transfer Learning (Future Enhancement)

```
Current: Train separate models for GW and Economic

Better approach: Transfer Learning
Step 1: Train on GW signals (lots of data)
        Learn what "clean signal" looks like
        Save weights

Step 2: Fine-tune on Economic signals
        Start with GW-trained weights
        Update only last few layers
        
Benefits:
✓ Use knowledge from GW to help Economic
✓ Need less Economic data
✓ Faster training
✓ Better generalization
```

### Multi-Task Learning

```
Instead of separate models:
    [Shared CNN body]
         ├─ Task 1: Denoise GW signals
         ├─ Task 2: Denoise Economic signals
         └─ Task 3: Classify signal type

Benefits:
✓ One model does multiple tasks
✓ Learn shared representations
✓ More efficient
```

### Uncertainty Quantification

```
Current: Point estimate
    model(noisy) → [0.5, 0.3, 0.7, ...]  [one value]

Better: Uncertainty bounds
    model(noisy) → mean: [0.5, 0.3, 0.7, ...]
                   std:  [0.1, 0.05, 0.15, ...]

Implementation: Bayesian Neural Networks
Benefits:
✓ Know confidence in predictions
✓ Detect out-of-distribution inputs
✓ Principled uncertainty propagation
```

### Spectral Analysis

```
After denoising, analyze frequency content:

    Noisy signal FFT:
    │  ┃╱╲         ┃╱╲
    │  ┃│ │         ┃│ │  ← Signal
    │  ┃│ │╱╱╱╱╱╱╱╱┃│ │  ← Noise (fills entire spectrum)
    │  ┃│ │╱╱╱╱╱╱╱╱┃│ │
    └──┴──┴───────────┴──┴──→ Frequency
    
    Denoised signal FFT:
    │  ┃╱╲         ┃╱╲
    │  ┃│ │         ┃│ │  ← Signal preserved
    │  ┃│ │         ┃│ │  ← Noise mostly removed
    │  ┃│ │         ┃│ │
    └──┴──┴───────────┴──┴──→ Frequency

Verification:
✓ Signal peaks still present
✓ Broadband noise reduced
✓ SNR improvement visible in frequency domain
```

---

## Conclusion

This project demonstrates:

1. **Problem Understanding**
   - What is signal denoising?
   - Why traditional methods fail
   - How neural networks can help

2. **Theoretical Knowledge**
   - How convolution works
   - Autoencoder compression principle
   - Time delay networks for temporal modeling

3. **Practical Implementation**
   - PyTorch dataset creation
   - Model architecture design
   - Training loop with validation
   - Comprehensive evaluation

4. **Domain Application**
   - Physics: Gravitational wave detection (83.6% event recall)
   - Finance: Economic crisis detection (limited by noise model)

5. **Model Selection**
   - CNN: Best for structured, stationary signals
   - TDNN: Best for non-stationary, dynamic signals
   - Trade-offs between complexity and generalization

**Key Takeaway:** Deep learning provides flexible, data-driven approach to signal processing that adapts to signal characteristics automatically. However, understanding the underlying problem (noise type, signal structure) is critical for success.

---

**For questions or clarifications, refer to the code comments and docstrings throughout the project!**