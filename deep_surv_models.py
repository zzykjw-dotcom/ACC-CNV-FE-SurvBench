############################################################
# deep_surv_models.py
# DeepSurv + Cox-Time under the SAME 5-fold splits and the
# SAME 59-dim input as Cox-AAE (56 CNV-driven genes + Age,
# Sex, Stage). Seed = 42 throughout.
#
# Usage (from project root or data folder):
#   python deep_surv_models.py
#
# Optional:
#   export ACC_DATA_DIR=/path/to/folder/with/csv
############################################################

import os
import json
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.model_selection import KFold
from lifelines.utils import concordance_index

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
DATA_DIR = os.environ.get("ACC_DATA_DIR", os.getcwd())

def _path(name):
    p = os.path.join(DATA_DIR, name)
    if os.path.isfile(p):
        return p
    p2 = os.path.join(DATA_DIR, "data", name)
    if os.path.isfile(p2):
        return p2
    return name


# ------------------------------------------------------------
# 1. Load data (identical logic to Cox_AAE.py)
# ------------------------------------------------------------
expr = pd.read_csv(_path("ACC_expression_data.csv"), index_col=0)
clin = pd.read_csv(_path("ACC_clinical.csv"))

core_path = _path("ME4_CNV_driven_core_genes.txt")
with open(core_path, "r") as f:
    core_genes = [
        line.strip() for line in f
        if line.strip() and not line.strip().startswith("#")
    ]

common_samples = list(set(expr.columns) & set(clin["sample_id"]))
print(f"[data] Intersected samples: {len(common_samples)}")

X_list = []
for g in core_genes:
    if g in expr.index:
        X_list.append(expr.loc[g, common_samples].values.astype(np.float64))
X_genes = np.array(X_list).T  # (n_samples, n_genes)

clin_sub = clin.set_index("sample_id").loc[common_samples]
X_clin = clin_sub[["Age", "Sex", "Stage"]].copy()
X_clin["Sex"] = X_clin["Sex"].apply(
    lambda x: 1 if str(x).upper() in ["MALE", "M"] else 0
)
X_clin["Stage"] = pd.factorize(X_clin["Stage"])[0]

X = np.hstack([X_genes, X_clin.values]).astype(np.float32)
T = clin_sub["OS_time"].values.astype(np.float32)
E = clin_sub["OS_status"].values.astype(np.float32)

mask = ~(np.isnan(X).any(axis=1) | np.isnan(T) | np.isnan(E))
X, T, E = X[mask], T[mask], E[mask]
print(f"[data] Final samples: {X.shape[0]}, features: {X.shape[1]}, events: {int(E.sum())}")

# Standardize features (fit on full set for simplicity; for stricter
# CV one could fit only on train fold — optional)
X_mean = X.mean(axis=0, keepdims=True)
X_std = X.std(axis=0, keepdims=True) + 1e-8
X = (X - X_mean) / X_std

# ------------------------------------------------------------
# 2. Fixed 5-fold splits (same seed as Cox-AAE)
# ------------------------------------------------------------
SEED = 42
np.random.seed(SEED)
torch.manual_seed(SEED)

kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
folds = []
for fold_idx, (tr, va) in enumerate(kf.split(X)):
    folds.append({"fold": fold_idx, "train": tr.tolist(), "val": va.tolist()})

with open(os.path.join(DATA_DIR, "folds_seed42.json"), "w") as f:
    json.dump(folds, f)
print(f"[folds] Saved folds_seed42.json ({len(folds)} folds)")


# ------------------------------------------------------------
# 3. Cox partial-likelihood loss (shared)
# ------------------------------------------------------------
def cox_ph_loss(risk, time, event):
    """
    risk: (N,) higher = higher hazard
    time, event: (N,)
    """
    order = torch.argsort(time, descending=True)
    risk = risk[order]
    event = event[order]
    log_cumsum = torch.logcumsumexp(risk, dim=0)
    return -torch.sum((risk - log_cumsum) * event) / (event.sum() + 1e-8)


# ------------------------------------------------------------
# 4. DeepSurv
# ------------------------------------------------------------
class DeepSurv(nn.Module):
    """Classic DeepSurv: MLP -> scalar risk (log-hazard)."""
    def __init__(self, input_dim, hidden=(32, 16), dropout=0.2):
        super().__init__()
        layers = []
        prev = input_dim
        for h in hidden:
            layers += [nn.Linear(prev, h), nn.ReLU(), nn.Dropout(dropout)]
            prev = h
        layers.append(nn.Linear(prev, 1, bias=False))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x).squeeze(-1)


def train_deepsurv_fold(X_tr, T_tr, E_tr, X_va, T_va, E_va,
                        epochs=100, lr=5e-4, device="cpu"):
    model = DeepSurv(X_tr.shape[1]).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)

    X_tr_t = torch.tensor(X_tr, device=device)
    T_tr_t = torch.tensor(T_tr, device=device)
    E_tr_t = torch.tensor(E_tr, device=device)
    X_va_t = torch.tensor(X_va, device=device)

    model.train()
    for _ in range(epochs):
        risk = model(X_tr_t)
        loss = cox_ph_loss(risk, T_tr_t, E_tr_t)
        opt.zero_grad()
        loss.backward()
        opt.step()

    model.eval()
    with torch.no_grad():
        risk_va = model(X_va_t).cpu().numpy()
    # concordance: higher risk -> shorter survival -> use -risk? 
    # lifelines concordance_index(event_times, predicted_scores, event_observed)
    # higher predicted_scores should mean longer survival for concordance,
    # so we pass -risk (lower risk = longer survival)
    ci = concordance_index(T_va, -risk_va, E_va)
    return model, ci


# ------------------------------------------------------------
# 5. Cox-Time style (time-dependent neural Cox)
#    Pure PyTorch: risk = f(x, t)  where t is scaled time.
#    Training uses discrete ranking over observed times (approximate
#    continuous-time Cox-Time). This does NOT require pycox.
# ------------------------------------------------------------
class CoxTimeNet(nn.Module):
    """
    Simplified Cox-Time: MLP takes concatenated [x, t_scaled] and
    outputs a time-dependent risk score.
    """
    def __init__(self, input_dim, hidden=(32, 16), dropout=0.2):
        super().__init__()
        layers = []
        prev = input_dim + 1  # + time
        for h in hidden:
            layers += [nn.Linear(prev, h), nn.ReLU(), nn.Dropout(dropout)]
            prev = h
        layers.append(nn.Linear(prev, 1, bias=False))
        self.net = nn.Sequential(*layers)

    def forward(self, x, t):
        # t: (N,) scaled to [0,1]
        xt = torch.cat([x, t.unsqueeze(-1)], dim=1)
        return self.net(xt).squeeze(-1)


def cox_time_loss(model, x, time, event, t_max):
    """
    Approximate continuous Cox partial likelihood with time-dependent
    risk: at each event time, risk is evaluated at that event's time.
    """
    order = torch.argsort(time, descending=True)
    x = x[order]
    time = time[order]
    event = event[order]
    t_scaled = (time / (t_max + 1e-8)).clamp(0, 1)

    risk = model(x, t_scaled)
    log_cumsum = torch.logcumsumexp(risk, dim=0)
    return -torch.sum((risk - log_cumsum) * event) / (event.sum() + 1e-8)


def train_coxtime_fold(X_tr, T_tr, E_tr, X_va, T_va, E_va,
                       epochs=100, lr=5e-4, device="cpu"):
    t_max = float(max(T_tr.max(), T_va.max()))
    model = CoxTimeNet(X_tr.shape[1]).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)

    X_tr_t = torch.tensor(X_tr, device=device)
    T_tr_t = torch.tensor(T_tr, device=device)
    E_tr_t = torch.tensor(E_tr, device=device)
    X_va_t = torch.tensor(X_va, device=device)
    T_va_t = torch.tensor(T_va, device=device)

    model.train()
    for _ in range(epochs):
        loss = cox_time_loss(model, X_tr_t, T_tr_t, E_tr_t, t_max)
        opt.zero_grad()
        loss.backward()
        opt.step()

    model.eval()
    with torch.no_grad():
        # For ranking at validation: evaluate risk at each sample's own time
        t_va_scaled = (T_va_t / (t_max + 1e-8)).clamp(0, 1)
        risk_va = model(X_va_t, t_va_scaled).cpu().numpy()
    ci = concordance_index(T_va, -risk_va, E_va)
    return model, ci


# ------------------------------------------------------------
# 6. Run 5-fold CV for both models
# ------------------------------------------------------------
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"[device] {device}")

results = {
    "DeepSurv": [],
    "CoxTime": [],
}

for fold in folds:
    tr, va = np.array(fold["train"]), np.array(fold["val"])
    X_tr, X_va = X[tr], X[va]
    T_tr, T_va = T[tr], T[va]
    E_tr, E_va = E[tr], E[va]

    # --- DeepSurv ---
    _, ci_ds = train_deepsurv_fold(
        X_tr, T_tr, E_tr, X_va, T_va, E_va,
        epochs=100, lr=5e-4, device=device
    )
    results["DeepSurv"].append(ci_ds)
    print(f"  Fold {fold['fold']+1} DeepSurv  C-index = {ci_ds:.4f}")

    # --- Cox-Time ---
    _, ci_ct = train_coxtime_fold(
        X_tr, T_tr, E_tr, X_va, T_va, E_va,
        epochs=100, lr=5e-4, device=device
    )
    results["CoxTime"].append(ci_ct)
    print(f"  Fold {fold['fold']+1} CoxTime   C-index = {ci_ct:.4f}")

# ------------------------------------------------------------
# 7. Summary
# ------------------------------------------------------------
print("\n========== SUMMARY ==========")
summary_rows = []
for name, scores in results.items():
    mean_ci = float(np.mean(scores))
    sd_ci = float(np.std(scores))
    print(f"{name:12s}: {mean_ci:.4f} ± {sd_ci:.4f}")
    summary_rows.append({
        "Algorithm": name,
        "Dim": 59,
        "C_index_mean": round(mean_ci, 4),
        "C_index_sd": round(sd_ci, 4),
        "fold_1": round(scores[0], 4),
        "fold_2": round(scores[1], 4),
        "fold_3": round(scores[2], 4),
        "fold_4": round(scores[3], 4),
        "fold_5": round(scores[4], 4),
        "Notes": "5-fold CV; same folds as Cox-AAE (seed=42); 59-dim input",
    })

# Reference line
print(f"\nReference: CoxPH optimism-corrected C-index = 0.882")
print(f"Reference: Cox-AAE (prior run)               = 0.831 ± 0.066")

out_csv = os.path.join(DATA_DIR, "deep_models_cindex_summary.csv")
pd.DataFrame(summary_rows).to_csv(out_csv, index=False)
print(f"\n[saved] {out_csv}")
print("[saved] folds_seed42.json")
print("Done.")
