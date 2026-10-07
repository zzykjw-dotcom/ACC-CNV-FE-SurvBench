############################################################
# Cox_AAE.py
# ALGORITHM TRACK – SurvBench deep model
# Cox adversarial autoencoder on 59-dimensional input
# (56 CNV-driven genes + Age, Sex, Stage)
# 5-fold cross-validated C-index
#
# Run from project root, or set environment variable:
#   export ACC_DATA_DIR=/path/to/folder/with/csv
############################################################

import os
import torch
import torch.nn as nn
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from lifelines.utils import concordance_index

# Data directory: ACC_DATA_DIR env, else current working directory
DATA_DIR = os.environ.get("ACC_DATA_DIR", os.getcwd())


def _path(name):
    p = os.path.join(DATA_DIR, name)
    if os.path.isfile(p):
        return p
    return name  # fallback to cwd


# ------------------------------------------------------------
# 1. Data loading
# ------------------------------------------------------------
expr = pd.read_csv(_path("ACC_expression_data.csv"), index_col=0)
clin = pd.read_csv(_path("ACC_clinical.csv"))

core_path = _path("ME4_CNV_driven_core_genes.txt")
if not os.path.isfile(core_path):
    core_path = _path(os.path.join("data", "ME4_CNV_driven_core_genes.txt"))
with open(core_path, "r") as f:
    core_genes = [
        line.strip() for line in f
        if line.strip() and not line.strip().startswith("#")
    ]

common_samples = list(set(expr.columns) & set(clin["sample_id"]))
print(f"Intersected samples: {len(common_samples)}")

X_list = []
for g in core_genes:
    if g in expr.index:
        X_list.append(expr.loc[g, common_samples].values)
X_genes = np.array(X_list).T

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
print(f"Final samples: {X.shape[0]}, features: {X.shape[1]}, events: {int(E.sum())}")


# ------------------------------------------------------------
# 2. Network
# ------------------------------------------------------------
class CoxAAE(nn.Module):
    def __init__(self, input_dim, latent_dim=16):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 32),
            nn.ReLU(),
            nn.Linear(32, latent_dim),
        )
        self.cox_head = nn.Linear(latent_dim, 1, bias=False)
        self.discriminator = nn.Sequential(
            nn.Linear(latent_dim, 16),
            nn.ReLU(),
            nn.Linear(16, 1),
            nn.Sigmoid(),
        )

    def forward(self, x):
        z = self.encoder(x)
        risk = self.cox_head(z)
        return risk, z


def cox_ph_loss(risk, time, event):
    order = torch.argsort(time, descending=True)
    risk = risk[order].squeeze()
    event = event[order]
    log_cumsum = torch.logcumsumexp(risk, dim=0)
    return -torch.sum((risk - log_cumsum) * event) / (event.sum() + 1e-8)


def adv_loss(discriminator, z, real=True):
    pred = discriminator(z)
    target = torch.ones_like(pred) if real else torch.zeros_like(pred)
    return nn.BCELoss()(pred, target)


# ------------------------------------------------------------
# 3. Training with 5-fold CV
# ------------------------------------------------------------
def train_cox_aae(X, T, E, n_splits=5, epochs=100, latent_dim=16,
                  lam=0.01, lr=5e-4):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    scores = []
    input_dim = X.shape[1]

    for fold, (tr, va) in enumerate(kf.split(X)):
        X_tr = torch.tensor(X[tr], device=device)
        T_tr = torch.tensor(T[tr], device=device)
        E_tr = torch.tensor(E[tr], device=device)
        X_va = torch.tensor(X[va], device=device)

        model = CoxAAE(input_dim, latent_dim).to(device)
        opt_enc = torch.optim.Adam(
            list(model.encoder.parameters()) + list(model.cox_head.parameters()),
            lr=lr,
        )
        opt_disc = torch.optim.Adam(model.discriminator.parameters(), lr=lr)

        for epoch in range(epochs):
            model.train()
            z_real = torch.randn(X_tr.size(0), latent_dim, device=device)
            with torch.no_grad():
                _, z_fake = model(X_tr)
            loss_d = adv_loss(model.discriminator, z_real, True) + \
                     adv_loss(model.discriminator, z_fake, False)
            opt_disc.zero_grad()
            loss_d.backward()
            opt_disc.step()

            risk, z = model(X_tr)
            loss = cox_ph_loss(risk, T_tr, E_tr) + \
                   lam * adv_loss(model.discriminator, z, True)
            opt_enc.zero_grad()
            loss.backward()
            opt_enc.step()

        model.eval()
        with torch.no_grad():
            risk_va, _ = model(X_va)
        ci = concordance_index(
            T[va], -risk_va.cpu().numpy().squeeze(), E[va]
        )
        scores.append(ci)
        print(f"Fold {fold + 1}: C-index = {ci:.4f}")

    print(f"\nMean C-index: {np.mean(scores):.4f} ± {np.std(scores):.4f}")
    return model, scores


if __name__ == "__main__":
    model, scores = train_cox_aae(
        X, T, E, n_splits=5, epochs=100, latent_dim=16, lam=0.01, lr=5e-4
    )
