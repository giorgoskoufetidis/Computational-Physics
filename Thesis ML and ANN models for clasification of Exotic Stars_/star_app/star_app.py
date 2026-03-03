# ==========================================================
# Star Classifier — Final Hybrid App (Scikit + PyTorch)
# ==========================================================
# Run: streamlit run star_classifier_app_final.py
# ==========================================================

import os, time, joblib, torch, torch.nn as nn
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import streamlit as st
import itertools
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    confusion_matrix, ConfusionMatrixDisplay, classification_report
)
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
import matplotlib.cm as cm
# -----------------------------
# Page Config
# -----------------------------
st.set_page_config(page_title="Neutron vs Quark Star Classifier", layout="wide")
st.title("🌟 Neutron vs Quark Star Classifier")

# -----------------------------
# Paths & Globals
# -----------------------------
PRETRAINED_PATH = "pretrained_models"
TRAINED_PATH = "trained_models"
os.makedirs(PRETRAINED_PATH, exist_ok=True)
os.makedirs(TRAINED_PATH, exist_ok=True)

LABEL_MAP = {0: "Neutron Star", 1: "Quark Star"}

# -----------------------------
# Utilities
# -----------------------------
def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [c.strip().lower() for c in df.columns]
    return df

def metric_table(y_true, y_pred):
    return {
        "Accuracy": accuracy_score(y_true, y_pred),
        "Precision": precision_score(y_true, y_pred, average="macro", zero_division=0),
        "Recall": recall_score(y_true, y_pred, average="macro", zero_division=0),
        "F1": f1_score(y_true, y_pred, average="macro", zero_division=0),
    }

def plot_confusion(y_true, y_pred, title):
    cm = confusion_matrix(y_true, y_pred)
    fig, ax = plt.subplots()
    ConfusionMatrixDisplay(cm, display_labels=[LABEL_MAP[0], LABEL_MAP[1]]).plot(ax=ax)
    plt.title(title)
    st.pyplot(fig)

def load_dataframe(upload):
    if upload is None:
        return None
    try:
        if upload.name.endswith(".csv"):
            return pd.read_csv(upload)
        else:
            return pd.read_excel(upload)
    except Exception as e:
        st.error(f"❌ Failed to read file: {e}")
        return None

def encode_k2_inplace(df: pd.DataFrame):
    """Encode k2 dynamically; if it's already numeric it's left as-is."""
    if "k2" not in df.columns:
        st.error("❌ Missing column 'k2'")
        st.stop()
    # If k2 is object/string -> encode


def prepare_features_for_model(df: pd.DataFrame):
    """
    - Ensures column order: mass, radius, k2
    - Encodes k2 if needed
    - Applies scaling only for KNN / Neural_Network (scikit) / any PyTorch NN
    """
    df = df.copy()
    encode_k2_inplace(df)
    df = df[["mass", "radius", "k2"]]
    return df.values


def load_or_fit_scaler_if_needed(scaler_path: str, X_sample: np.ndarray | None):
    """
    - If scaler exists and is fitted -> return it
    - If exists but not fitted and X_sample provided -> fit, save, return
    - If missing and X_sample provided -> fit new, save, return
    - Otherwise return None (prediction can continue for unscaled models)
    """
    try:
        if os.path.exists(scaler_path):
            scaler = joblib.load(scaler_path)
            if hasattr(scaler, "mean_") and hasattr(scaler, "scale_"):
                return scaler
            if X_sample is not None:
                st.info("🔄 Refitting existing scaler with sample data...")
                scaler.fit(X_sample)
                joblib.dump(scaler, scaler_path)
                st.success("✅ Scaler re-fitted and saved.")
                return scaler
            st.warning("⚠️ Scaler not fitted and no sample data provided.")
            return None
        else:
            if X_sample is not None:
                st.info("🆕 Creating and fitting a new scaler with sample data...")
                scaler = StandardScaler().fit(X_sample)
                joblib.dump(scaler, scaler_path)
                st.success("✅ New scaler saved.")
                return scaler
            return None
    except Exception as e:
        st.error(f"❌ Failed to load/fit scaler: {e}")
        return None

# -----------------------------
# PyTorch model class
# -----------------------------
class FlexibleModel(nn.Module):
    def __init__(self, in_features, hidden_layers, out_features):
        super(FlexibleModel, self).__init__()
        layers = []
        prev_units = in_features
        for h in hidden_layers:
            layers.append(nn.Linear(prev_units, h))
            layers.append(nn.ReLU())
            prev_units = h
        layers.append(nn.Linear(prev_units, out_features))
        self.network = nn.Sequential(*layers)
    def forward(self, x):
        return self.network(x)

# -----------------------------
# Sidebar: Model Source + Banner + Info
# -----------------------------
st.sidebar.header("⚙️ Model Source")
model_source = st.sidebar.radio(
    "Select which models to use:",
    ["Pretrained Models", "Newly Trained Models"],
    index=0
)
ACTIVE_PATH = PRETRAINED_PATH if model_source == "Pretrained Models" else TRAINED_PATH
os.makedirs(ACTIVE_PATH, exist_ok=True)

# Visual banner
color = "#22c55e" if model_source == "Pretrained Models" else "#8b5cf6"
label = "🟢 Pretrained Models Active" if model_source == "Pretrained Models" else "🟣 Trained Models Active"
st.markdown(
    f"<div style='background-color:{color};padding:10px;border-radius:10px;margin-bottom:10px;text-align:center;color:white;'>{label}</div>",
    unsafe_allow_html=True,
)

# Show folder and models with timestamps
st.sidebar.markdown(f"📂 **Active Folder:** `{ACTIVE_PATH}/`")
model_files_in_folder = [
    f for f in os.listdir(ACTIVE_PATH)
    if (f.endswith(".joblib") or f.endswith(".pth"))
    and all(x not in f for x in ["scaler", "encoder", "arch"])
]
if model_files_in_folder:
    st.sidebar.success(f"✅ Models found ({len(model_files_in_folder)}):")
    for f in model_files_in_folder:
        full = os.path.join(ACTIVE_PATH, f)
        ts = time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(full)))
        st.sidebar.markdown(f"- {f} _(updated {ts})_")
else:
    st.sidebar.warning("⚠️ No model files found in this folder.")

# Scaler path is tied to the ACTIVE folder
ACTIVE_SCALER_PATH_NN = os.path.join(ACTIVE_PATH, "scaler.joblib")
ACTIVE_SCALER_PATH_ML = os.path.join(ACTIVE_PATH, "scaler_ml.joblib")
# -----------------------------
# Load models from ACTIVE_PATH
# -----------------------------
loaded_models: dict[str, dict] = {}

for file in model_files_in_folder:
    name = file.replace(".joblib", "").replace(".pth", "")
    path = os.path.join(ACTIVE_PATH, file)
    try:
        # PyTorch model by .pth + architecture sidecar
        if file.endswith(".pth"):
            arch_path = os.path.join(ACTIVE_PATH, f"{name}_arch.joblib")
            if not os.path.exists(arch_path):
                st.sidebar.warning(f"⚠️ Missing architecture for {name}, skipping.")
                continue
            params = joblib.load(arch_path)
            nn_model = FlexibleModel(**params)
            nn_model.load_state_dict(torch.load(path, map_location="cpu"))
            nn_model.eval()
            scaler = joblib.load(ACTIVE_SCALER_PATH_NN) if os.path.exists(ACTIVE_SCALER_PATH_NN) else None
            loaded_models[name] = {"type": "torch", "model": nn_model, "scaler": scaler}
        else:
            # scikit-learn model (.joblib)
            obj = joblib.load(path)
            if isinstance(obj, dict) and "model_state_dict" in obj:
                st.sidebar.warning(f"⚠️ {name} looks like a raw torch state dict; skipping.")
                continue
            scaler = joblib.load(ACTIVE_SCALER_PATH_ML) if os.path.exists(ACTIVE_SCALER_PATH_ML) else None
            loaded_models[name] = {"type": "sklearn", "model": obj, "scaler": scaler}
    except Exception as e:
        st.sidebar.error(f"❌ Could not load {file}: {e}")

if loaded_models:
    st.sidebar.success(f"✅ Loaded Models: {', '.join(loaded_models.keys())}")
else:
    st.sidebar.warning("⚠️ No models loaded.")

# -----------------------------
# TABS
# -----------------------------
tab_predict, tab_train = st.tabs(["🔮 Predict", "🧠 Train / Retrain"])

# ==========================================================
# 🔮 PREDICT TAB
# ==========================================================
with tab_predict:
    st.subheader("Manual Prediction (Single Value)")

    # Model selector for single value
    model_choice_single = st.selectbox(
        "Choose Model for Single Prediction",
        list(loaded_models.keys()) if loaded_models else [],
        key="model_choice_single"
    )

    c1, c2, c3 = st.columns(3)
    with c1: mass = st.number_input("Mass", value=1.0, format="%.6f")
    with c2: radius = st.number_input("Radius", value=10.0, format="%.6f")
    with c3: k2_val = st.number_input("K2", value=0.1, format="%.6f")

    if st.button("Predict Single Value"):
        if not loaded_models or not model_choice_single:
            st.error("❌ No model available.")
        else:
            sample = pd.DataFrame({"mass": [mass], "radius": [radius], "k2": [k2_val]})
            X_prepared = prepare_features_for_model(sample)
            info = loaded_models[model_choice_single]
            if info["type"] == "sklearn":
                X_scaled = info["scaler"].transform(X_prepared)
                probs = info["model"].predict_proba(X_scaled)[0]
            else:
                nn_model = info["model"]
                scaler = info["scaler"]
                X_scaled = scaler.transform(X_prepared)
                with torch.no_grad():
                    logits = nn_model(torch.tensor(X_scaled, dtype=torch.float32))
                    p = torch.sigmoid(logits).numpy().flatten()
                    probs = np.array([1 - p[0], p[0]])
            pred = int(np.argmax(probs))
            # ----------------------------------------------------------
            # Bar chart for single prediction probabilities
            # ----------------------------------------------------------
            st.markdown("### 📊 Prediction Probabilities")

            labels = ["Neutron Star", "Quark Star"]
            values = [probs[0], probs[1]]

            fig_prob, ax_prob = plt.subplots(figsize=(3, 2.5))

            bars = ax_prob.bar(labels, values, color=["#3b82f6", "#ef4444"], width=0.5)

            # Add percentage labels inside bars
            for bar, val in zip(bars, values):
                ax_prob.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() / 2,
                    f"{val*100:.1f}%",
                    ha="center",
                    va="center",
                    color="black",
                    fontsize=5,
                    fontweight="bold"
                )

            ax_prob.set_ylim(0, 1)
            ax_prob.set_ylabel("Probability")
            ax_prob.set_title("Single Prediction Probabilities")
            ax_prob.grid(axis="y", linestyle="--", alpha=0.4)

            st.pyplot(fig_prob, use_container_width=False)

            st.success(f"🌠 Predicted: **{LABEL_MAP[pred]}** (using {model_choice_single})")

    # -----------------------------
    # Range Prediction Mode
    # -----------------------------
    st.markdown("---")
    st.subheader("📊 Range Prediction (Multiple Combinations)")

    model_choice_range = st.selectbox(
        "Choose Model for Range Prediction",
        list(loaded_models.keys()) if loaded_models else [],
        key="model_choice_range"
    )

    col1, col2, col3 = st.columns(3)
    with col1:
        mass_min = st.number_input("Mass min", value=1.0)
        mass_max = st.number_input("Mass max", value=1.5)
        step_mass = st.number_input("Step for Mass", value=0.1)
    with col2:
        radius_min = st.number_input("Radius min", value=9.0)
        radius_max = st.number_input("Radius max", value=12.0)
        step_radius = st.number_input("Step for Radius", value=0.1)
    with col3:
        k2_min = st.number_input("K2 min", value=0.05)
        k2_max = st.number_input("K2 max", value=0.15)
        step_k2 = st.number_input("Step for K2", value=0.05)

    if st.button("🚀 Run Range Prediction"):
        masses = np.arange(mass_min, mass_max + step_mass / 2, step_mass)
        radii = np.arange(radius_min, radius_max + step_radius / 2, step_radius)
        k2s = np.arange(k2_min, k2_max + step_k2 / 2, step_k2)

        grid = list(itertools.product(masses, radii, k2s))
        df_grid = pd.DataFrame(grid, columns=["mass", "radius", "k2"])
        st.info(f"🧮 Generated {len(df_grid)} combinations.")

        fig = plt.figure(figsize=(4, 3))
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(df_grid["mass"], df_grid["radius"], df_grid["k2"], color="#4fa3f7", s=10, alpha=0.7, edgecolors="k")
        ax.set_xlabel("Mass")
        ax.set_ylabel("Radius")
        ax.set_zlabel("k2")
        ax.set_title("Generated Parameter Grid")
        st.pyplot(fig, use_container_width=False)

        if not model_choice_range:
            st.error("❌ Please select a model for range prediction.")
            st.stop()

        info = loaded_models[model_choice_range]
        X_prepared = prepare_features_for_model(df_grid)

        if info["type"] == "sklearn":
            X_scaled = info["scaler"].transform(X_prepared)
            probs = info["model"].predict_proba(X_scaled)
        else:
            nn_model = info["model"]
            scaler = info["scaler"]
            X_scaled = scaler.transform(X_prepared)
            with torch.no_grad():
                logits = nn_model(torch.tensor(X_scaled, dtype=torch.float32))
                p = torch.sigmoid(logits).numpy().flatten()
                probs = np.column_stack([1 - p, p])

        preds = np.argmax(probs, axis=1)
        df_grid["Predicted_Type"] = [LABEL_MAP[p] for p in preds]
        df_grid["Prob_Neutron"] = probs[:, 0]
        df_grid["Prob_Quark"] = probs[:, 1]

        st.success("✅ Range predictions complete!")
        st.dataframe(df_grid.head())

        # 3D Plot of Predicted Classes
        col_plot3d, col_barchart = st.columns([2, 1])
        color_map = {"Neutron Star": "#3b82f6", "Quark Star": "#ef4444"}
        colors = df_grid["Predicted_Type"].map(color_map)
        with col_plot3d:
            st.markdown("### 🌌 3D Visualization of Predicted Classes")
            fig2 = plt.figure(figsize=(4, 3))
            ax2 = fig2.add_subplot(111, projection="3d")
            ax2.scatter(df_grid["mass"], df_grid["radius"], df_grid["k2"], c=colors, s=15, alpha=0.8, edgecolors="k")
            for label, color in color_map.items():
                ax2.scatter([], [], [], color=color, label=label, s=25)
            ax2.legend(loc="upper right", fontsize=7)
            ax2.set_xlabel("Mass")
            ax2.set_ylabel("Radius")
            ax2.set_zlabel("k2")
            ax2.set_title("Predicted Classes (Neutron vs Quark)")
            st.pyplot(fig2, use_container_width=False)

        with col_barchart:
            st.markdown("### 📊 Prediction Totals")
            counts = df_grid["Predicted_Type"].value_counts().reindex(["Neutron Star", "Quark Star"], fill_value=0)
            fig_bar, ax_bar = plt.subplots(figsize=(4, 3))
            bars = ax_bar.bar(counts.index, counts.values, color=["#3b82f6", "#ef4444"], width=0.5)
            for bar in bars:
                ax_bar.text(bar.get_x() + bar.get_width() / 2, bar.get_height() / 2,
                            f"{int(bar.get_height())}", ha="center", va="center", color="white", fontsize=8, fontweight="bold")
            ax_bar.set_ylabel("Count")
            ax_bar.set_title("Total Predictions")
            ax_bar.grid(axis="y", linestyle="--", alpha=0.4)
            st.pyplot(fig_bar, use_container_width=True)

        csv = df_grid.to_csv(index=False).encode("utf-8")
        st.download_button("📥 Download Predictions", csv, "grid_predictions.csv", "text/csv")





    st.markdown("---")
    st.subheader("📂 Batch Prediction")
    upload = st.file_uploader("Upload CSV/XLSX with columns: mass, radius, k2", type=["csv", "xlsx", "xls"])
    batch_model = st.selectbox("Select Model for Batch Prediction", list(loaded_models.keys()), key="batch_model")

    if st.button("🚀 Predict Batch File"):
        if upload is None:
            st.error("❌ Please upload a file.")
        else:
            df = load_dataframe(upload)
            if df is None: st.stop()
            df = normalize_columns(df)
            required = {"mass", "radius", "k2"}
            if not required.issubset(df.columns):
                st.error("❌ Missing required columns: mass, radius, k2")
            else:
                # Ensure scaler for models that need it
                needs_scaling = (loaded_models[batch_model]["type"] == "torch") or (batch_model in ["Neural_Network"])
                # if needs_scaling:
                #     _ = load_or_fit_scaler_if_needed(ACTIVE_SCALER_PATH, X_sample=df.assign(k2=[0]*len(df)).values)

                X_prepared = prepare_features_for_model(df)
                info = loaded_models[batch_model]

                if info["type"] == "sklearn":
                    X_scaled = info["scaler"].transform(X_prepared)
                    probs = info["model"].predict_proba(X_scaled)
                else:
                    nn_model = info["model"]
                    scaler = info["scaler"]
                    if scaler is None or not hasattr(scaler, "mean_"):
                        st.error("❌ Scaler missing/unfitted for the neural network in this folder.")
                        st.stop()
                    X_scaled = scaler.transform(X_prepared)
                    with torch.no_grad():
                        logits = nn_model(torch.tensor(X_scaled, dtype=torch.float32))
                        p = torch.sigmoid(logits).numpy().reshape(-1, 1)
                        probs = np.hstack([1 - p, p])
                preds = np.argmax(probs, axis=1)
                out = df.copy()
                out["Predicted_Type"] = [LABEL_MAP[p] for p in preds]
                out["Prob_Neutron"] = probs[:, 0]
                out["Prob_Quark"] = probs[:, 1]
                st.success(f"✅ Predictions completed using {batch_model}")
                st.dataframe(out.head())
                st.download_button("📥 Download Predictions", out.to_csv(index=False).encode("utf-8"),
                                   "batch_predictions.csv", "text/csv")

# ==========================================================
# 🧠 TRAIN / RETRAIN TAB
# ==========================================================
with tab_train:
    st.subheader("Upload Labeled Data for Training / Retraining")
    file = st.file_uploader("Upload CSV/XLSX (mass, radius, k2, star_type)", type=["csv", "xlsx", "xls"])

    if file:
        df = load_dataframe(file)
        if df is not None:
            df = normalize_columns(df)
            st.dataframe(df.head())
            required = {"mass", "radius", "k2", "star_type"}

            if not required.issubset(df.columns):
                st.error("❌ Missing required columns")
            else:
                # Encode target and k2
                y = LabelEncoder().fit_transform(df["star_type"].astype(str))

                # Dynamic k2 encoding (no persistent encoder dependency)
                encode_k2_inplace(df)

                X = df[["mass", "radius", "k2"]]
                X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, test_size=0.2, random_state=42)

                # Fit scaler for models that need scaling and SAVE IT TO TRAINED_PATH
                scaler_new = StandardScaler()
                X_train_s = scaler_new.fit_transform(X_train)
                X_test_s = scaler_new.transform(X_test)
                ACTIVE_TRAIN_SCALER = os.path.join(TRAINED_PATH, "scaler.joblib")
                joblib.dump(scaler_new, ACTIVE_TRAIN_SCALER)
                st.info("💾 Saved scaler for trained models.")

                st.info("🔄 Training models (GridSearchCV, 5-fold CV)...")
                cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

                grids = {
                    "Decision_Tree": {
                        "estimator": DecisionTreeClassifier(random_state=0),
                        "grid": {"max_depth": [None, 5, 10, 20]},
                        "use_scaled": False
                    },
                    "Random_Forest": {
                        "estimator": RandomForestClassifier(random_state=0),
                        "grid": {"n_estimators": [100, 200], "max_depth": [None, 10, 20]},
                        "use_scaled": False
                    },
                    "KNN": {
                        "estimator": KNeighborsClassifier(),
                        "grid": {"n_neighbors": [3, 5, 7], "weights": ["uniform", "distance"]},
                        "use_scaled": True
                    },
                    "Neural_Network": {
                        "estimator": MLPClassifier(max_iter=1000, random_state=0),
                        "grid": {"hidden_layer_sizes": [(50,), (100,)], "activation": ["relu", "tanh"]},
                        "use_scaled": True
                    }
                }

                results = []
                progress = st.progress(0)
                best_acc = -1.0
                best_model_name = None

                # Train & save into TRAINED_PATH (won't overwrite pretrained)
                for i, (name, spec) in enumerate(grids.items(), start=1):
                    Xtr = X_train_s if spec["use_scaled"] else X_train.values
                    Xte = X_test_s if spec["use_scaled"] else X_test.values

                    gs = GridSearchCV(spec["estimator"], spec["grid"], cv=cv, scoring="accuracy", n_jobs=-1)
                    gs.fit(Xtr, y_train)

                    best = gs.best_estimator_
                    y_pred = best.predict(Xte)
                    metrics = metric_table(y_test, y_pred)
                    results.append({"Model": name, **metrics})

                    joblib.dump(best, os.path.join(TRAINED_PATH, f"{name}.joblib"))
                    plot_confusion(y_test, y_pred, f"{name} Confusion Matrix (Test)")

                    if metrics["Accuracy"] > best_acc:
                        best_acc = metrics["Accuracy"]
                        best_model_name = name
                        joblib.dump(best, os.path.join(TRAINED_PATH, "best_model.joblib"))

                    progress.progress(i / len(grids))
                    time.sleep(0.2)

                st.success("✅ Training complete and models saved to `trained_models/`.")
                res_df = pd.DataFrame(results).sort_values("Accuracy", ascending=False)
                st.dataframe(res_df.style.format("{:.4f}"))
                if best_model_name:
                    st.success(f"🏆 Best Model: {best_model_name} (Acc={best_acc:.4f})")
