"""Frozen Stage24 numerical functions, extracted without importing its driver.

Algorithm bodies are unchanged. The portable runner supplies input/output paths.
The original non-general aggregate/bootstrap drivers are replaced separately.
"""
from __future__ import annotations
import hashlib,math
from pathlib import Path
from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
ALPHAS=np.asarray([1e-4,1e-3,1e-2,1e-1,1.,10.,100.,1000.],float)
RIDGE_SOLVER="cholesky"
MIN_ROWS=50


def atomic_csv(frame: pd.DataFrame, path: Path, compress: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + (".tmp.gz" if compress else ".tmp"))
    compression = {"method": "gzip", "mtime": 0} if compress else None
    frame.to_csv(tmp, index=False, encoding="utf-8-sig", compression=compression)
    tmp.replace(path)

@dataclass
class Bundle:
    participant: np.ndarray
    record: np.ndarray
    proteins: list[str]
    x: np.ndarray
    position: dict[str, int]

def read_bundle(path: Path) -> Bundle:
    frame = pd.read_csv(path)
    if list(frame.columns[:2]) != ["participant_id", "record_id"]:
        raise RuntimeError(f"Identifier columns changed: {path}")
    proteins = frame.columns[2:].astype(str).tolist()
    x = frame.iloc[:, 2:].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    if np.isinf(x).any() or len(proteins) != len(set(proteins)):
        raise RuntimeError(f"Invalid quantitative matrix: {path}")
    return Bundle(
        participant=frame.participant_id.astype(str).to_numpy(),
        record=frame.record_id.astype(str).to_numpy(),
        proteins=proteins,
        x=x,
        position={p: i for i, p in enumerate(proteins)},
    )

def weighted_score(bundle: Bundle, weights: pd.Series, selected: list[str]) -> np.ndarray:
    chosen = [p for p in selected if p in weights.index and p in bundle.position]
    if len(chosen) < 2:
        return np.full(len(bundle.participant), np.nan)
    w = weights.loc[chosen].to_numpy(float)
    denom = float(np.abs(w).sum())
    if not np.isfinite(denom) or denom <= 0:
        return np.full(len(bundle.participant), np.nan)
    pos = [bundle.position[p] for p in chosen]
    return bundle.x[:, pos] @ w / denom

def feature_matrix(bundle: Bundle, selected: list[str]) -> np.ndarray:
    if not all(p in bundle.position for p in selected):
        raise RuntimeError("Observed protein absent from HMACS batch")
    return bundle.x[:, [bundle.position[p] for p in selected]]

def fold_seed(base: int, tag: str) -> int:
    value = int(hashlib.sha256(f"{base}|{tag}".encode("utf-8")).hexdigest()[:8], 16)
    return 1 + value % 2147483646

def select_alpha_one_se(x: np.ndarray, y: np.ndarray, seed: int) -> tuple[float, pd.DataFrame]:
    if len(y) < MIN_ROWS or x.shape[1] < 1:
        raise RuntimeError("Inadequate ridge training dimensions")
    splits = min(5, max(2, len(y) // 20))
    cv = KFold(n_splits=splits, shuffle=True, random_state=seed)
    errors = np.empty((len(ALPHAS), splits), float)
    for fold, (train, valid) in enumerate(cv.split(x)):
        for ai, alpha in enumerate(ALPHAS):
            model = Ridge(alpha=float(alpha), fit_intercept=True, solver=RIDGE_SOLVER)
            model.fit(x[train], y[train])
            prediction = model.predict(x[valid])
            errors[ai, fold] = float(np.mean((y[valid] - prediction) ** 2))
    means = errors.mean(axis=1)
    ses = errors.std(axis=1, ddof=1) / math.sqrt(splits)
    best = int(np.argmin(means))
    eligible = np.where(means <= means[best] + ses[best] + 1e-15)[0]
    chosen = int(eligible[-1])
    audit = pd.DataFrame({"alpha": ALPHAS, "mean_cv_mse": means, "se_cv_mse": ses, "selected": np.arange(len(ALPHAS)) == chosen})
    return float(ALPHAS[chosen]), audit

def within_crossfit(x: np.ndarray, residual: np.ndarray, valid: np.ndarray, seed: int) -> tuple[np.ndarray, list[float], list[dict[str, object]]]:
    prediction = np.full(len(residual), np.nan)
    indices = np.where(valid)[0]
    if len(indices) < MIN_ROWS:
        raise RuntimeError("Fewer than 50 target rows for within-batch cross-fitting")
    cv = KFold(n_splits=5, shuffle=True, random_state=seed)
    alphas, audits = [], []
    for outer_fold, (outer_train, outer_test) in enumerate(cv.split(indices), 1):
        train_idx, test_idx = indices[outer_train], indices[outer_test]
        alpha, _ = select_alpha_one_se(x[train_idx], residual[train_idx], fold_seed(seed, f"outer{outer_fold}"))
        model = Ridge(alpha=alpha, fit_intercept=True, solver=RIDGE_SOLVER)
        model.fit(x[train_idx], residual[train_idx])
        prediction[test_idx] = model.predict(x[test_idx])
        alphas.append(alpha)
        audits.append(
            {
                "outer_fold": outer_fold,
                "train_n": len(train_idx),
                "test_n": len(test_idx),
                "alpha": alpha,
                "train_test_overlap": len(set(train_idx).intersection(set(test_idx))),
            }
        )
    return prediction, alphas, audits

def fit_cross_batch_ridge(train_x: np.ndarray, train_residual: np.ndarray, train_valid: np.ndarray, test_x: np.ndarray, test_valid: np.ndarray, seed: int) -> tuple[np.ndarray, float, pd.DataFrame]:
    if train_valid.sum() < MIN_ROWS or test_valid.sum() < MIN_ROWS:
        raise RuntimeError("Fewer than 50 complete source or target rows")
    alpha, audit = select_alpha_one_se(train_x[train_valid], train_residual[train_valid], seed)
    model = Ridge(alpha=alpha, fit_intercept=True, solver=RIDGE_SOLVER)
    model.fit(train_x[train_valid], train_residual[train_valid])
    prediction = np.full(len(test_valid), np.nan)
    prediction[test_valid] = model.predict(test_x[test_valid])
    return prediction, alpha, audit

def pearson(x: np.ndarray, y: np.ndarray) -> float:
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 3 or np.std(x[keep]) <= 0 or np.std(y[keep]) <= 0:
        return math.nan
    return float(np.corrcoef(x[keep], y[keep])[0, 1])

def spearman(x: np.ndarray, y: np.ndarray) -> float:
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 3:
        return math.nan
    return pearson(rankdata(x[keep], method="average"), rankdata(y[keep], method="average"))

def ccc(x: np.ndarray, y: np.ndarray) -> float:
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 3:
        return math.nan
    a, b = x[keep], y[keep]
    va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
    cov = np.cov(a, b, ddof=1)[0, 1]
    den = va + vb + (np.mean(a) - np.mean(b)) ** 2
    return float(2 * cov / den) if den > 0 else math.nan

def metrics(reference: np.ndarray, direct: np.ndarray, projected: np.ndarray) -> dict[str, float | int]:
    keep = np.isfinite(reference) & np.isfinite(direct) & np.isfinite(projected)
    if keep.sum() < MIN_ROWS:
        return {k: math.nan for k in ("risk", "srmse", "pearson", "spearman", "ccc", "calibration_intercept", "calibration_slope", "alignment_benefit", "correction_cost", "net_gain_identity")} | {"n": int(keep.sum())}
    f, d, p = reference[keep], direct[keep], projected[keep]
    vf = float(np.var(f, ddof=1))
    if not np.isfinite(vf) or vf <= 0:
        return {k: math.nan for k in ("risk", "srmse", "pearson", "spearman", "ccc", "calibration_intercept", "calibration_slope", "alignment_benefit", "correction_cost", "net_gain_identity")} | {"n": int(keep.sum())}
    risk = float(np.mean((f - p) ** 2) / vf)
    design = np.column_stack([np.ones(len(p)), p])
    coef = np.linalg.lstsq(design, f, rcond=None)[0]
    e, correction = f - d, p - d
    benefit = float(2 * np.mean(e * correction) / vf)
    cost = float(np.mean(correction**2) / vf)
    return {
        "n": int(len(f)),
        "risk": risk,
        "srmse": math.sqrt(max(risk, 0)),
        "pearson": pearson(f, p),
        "spearman": spearman(f, p),
        "ccc": ccc(f, p),
        "calibration_intercept": float(coef[0]),
        "calibration_slope": float(coef[1]),
        "alignment_benefit": benefit,
        "correction_cost": cost,
        "net_gain_identity": benefit - cost,
    }

def compute_predictions(mode: str, out: Path, bundles: dict[str, Bundle], weights: pd.DataFrame, masks: pd.DataFrame, units: pd.DataFrame, seeds: pd.DataFrame, crnp_scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mask_lookup = {k: g.protein.astype(str).tolist() for k, g in masks.groupby("mask_id", sort=False)}
    seed_lookup = dict(zip(seeds.stream.astype(str), seeds.seed.astype(int)))
    crnp_lookup = {k: g for k, g in crnp_scores.groupby("direction_unit_id", sort=False)}
    predictions, failures, model_audit, fold_audit = [], [], [], []

    for number, unit in enumerate(units.itertuples(index=False), 1):
        try:
            wd = weights[weights.module.eq(unit.module)]
            w = wd.set_index("protein").weight.astype(float)
            full = wd.protein.astype(str).tolist()
            observed = [p for p in mask_lookup[unit.mask_unit_id] if p in w.index]
            train, test = bundles[unit.train_batch], bundles[unit.test_batch]
            train_full = weighted_score(train, w, full)
            train_direct = weighted_score(train, w, observed)
            test_full = weighted_score(test, w, full)
            test_direct = weighted_score(test, w, observed)
            train_x, test_x = feature_matrix(train, observed), feature_matrix(test, observed)
            train_residual = train_full - train_direct
            test_residual = test_full - test_direct
            train_valid = np.isfinite(train_residual) & np.isfinite(train_x).all(axis=1)
            test_valid = np.isfinite(test_residual) & np.isfinite(test_x).all(axis=1)
            unit_seed = int(seed_lookup[f"fold::{unit.direction_unit_id}"])

            within_correction, outer_alphas, outer_audits = within_crossfit(test_x, test_residual, test_valid, unit_seed)
            cross_correction, cross_alpha, cv_audit = fit_cross_batch_ridge(
                train_x, train_residual, train_valid, test_x, test_valid, unit_seed
            )
            within_score = test_direct + within_correction
            cross_score = test_direct + cross_correction

            if unit.direction_unit_id not in crnp_lookup:
                raise RuntimeError("CRNP score output absent")
            crnp = crnp_lookup[unit.direction_unit_id]
            crnp = crnp.set_index(["participant_id", "record_id"])
            crnp_score = np.asarray(
                [
                    pd.to_numeric(crnp.loc[(p, r), "CRNP_v0.2.0"], errors="coerce")
                    if (p, r) in crnp.index
                    else math.nan
                    for p, r in zip(test.participant, test.record)
                ],
                float,
            )

            frame = pd.DataFrame(
                {
                    "direction_unit_id": unit.direction_unit_id,
                    "mask_unit_id": unit.mask_unit_id,
                    "module": unit.module,
                    "platform_label": unit.platform_label,
                    "direction": unit.direction,
                    "train_batch": unit.train_batch,
                    "test_batch": unit.test_batch,
                    "participant_id": test.participant,
                    "record_id": test.record,
                    "reference_full": test_full,
                    "direct": test_direct,
                    "within_cf_ridge": within_score,
                    "cross_batch_ridge": cross_score,
                    "CRNP_v0.2.0": crnp_score,
                }
            )
            predictions.append(frame)
            model_audit.append(
                {
                    "direction_unit_id": unit.direction_unit_id,
                    "module": unit.module,
                    "direction": unit.direction,
                    "train_n_complete": int(train_valid.sum()),
                    "test_n_complete": int(test_valid.sum()),
                    "full_proteins": len(full),
                    "observed_proteins": len(observed),
                    "cross_batch_alpha": cross_alpha,
                    "within_outer_alphas": ";".join(str(x) for x in outer_alphas),
                    "status": "PASS",
                }
            )
            for audit in outer_audits:
                fold_audit.append({"direction_unit_id": unit.direction_unit_id, **audit})
            cv_audit.insert(0, "direction_unit_id", unit.direction_unit_id)
            atomic_csv(cv_audit, out / "00_preflight" / "ridge_cv" / f"{number:03d}.csv")
        except Exception as exc:
            failures.append(
                {
                    "direction_unit_id": unit.direction_unit_id,
                    "module": unit.module,
                    "direction": unit.direction,
                    "stage": "prediction_construction",
                    "detail": repr(exc),
                    "status": "FAILED_RETAINED",
                }
            )
            model_audit.append(
                {
                    "direction_unit_id": unit.direction_unit_id,
                    "module": unit.module,
                    "direction": unit.direction,
                    "status": "FAILED_RETAINED",
                }
            )
        if number == 1 or number % 12 == 0 or number == len(units):
            print(f"prediction units {number}/{len(units)}", flush=True)

    pred = pd.concat(predictions, ignore_index=True) if predictions else pd.DataFrame()
    fail = pd.DataFrame(failures, columns=["direction_unit_id", "module", "direction", "stage", "detail", "status"])
    return pred, fail, pd.DataFrame(model_audit), pd.DataFrame(fold_audit)

def cell_metrics(predictions: pd.DataFrame, units: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for unit_id, frame in predictions.groupby("direction_unit_id", sort=True):
        meta = units.set_index("direction_unit_id").loc[unit_id]
        direct_metric = metrics(frame.reference_full.to_numpy(float), frame.direct.to_numpy(float), frame.direct.to_numpy(float))
        method_metrics = {"direct": direct_metric}
        for method in ("within_cf_ridge", "cross_batch_ridge", "CRNP_v0.2.0"):
            method_metrics[method] = metrics(
                frame.reference_full.to_numpy(float),
                frame.direct.to_numpy(float),
                frame[method].to_numpy(float),
            )
        risk_direct = float(direct_metric["risk"])
        headroom = risk_direct - float(method_metrics["within_cf_ridge"]["risk"])
        for method, metric in method_metrics.items():
            gain = risk_direct - float(metric["risk"])
            row = {
                "direction_unit_id": unit_id,
                "mask_unit_id": meta.mask_unit_id,
                "module": meta.module,
                "platform_label": meta.platform_label,
                "direction": meta.direction,
                "train_batch": meta.train_batch,
                "test_batch": meta.test_batch,
                "method": method,
                **metric,
                "direct_risk": risk_direct,
                "empirical_headroom": headroom,
                "method_gain": gain,
                "relative_gain": gain / risk_direct if risk_direct > 0 else math.nan,
                "transport_loss": float(method_metrics["cross_batch_ridge"]["risk"] - method_metrics["within_cf_ridge"]["risk"]),
                "captured_headroom": gain / headroom if headroom > 0 and method != "direct" else math.nan,
                "headroom_status": "POSITIVE_EMPIRICAL_HEADROOM" if headroom > 0 else "NO_POSITIVE_EMPIRICAL_HEADROOM",
                "status": "PASS" if np.isfinite(metric["risk"]) else "FAILED_RETAINED",
            }
            identity_diff = gain - float(metric["net_gain_identity"])
            row["gain_identity_difference"] = identity_diff
            rows.append(row)
    return pd.DataFrame(rows)

def participant_bootstrap_indices(bundle: Bundle, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    unique = pd.unique(bundle.participant)
    sampled = rng.choice(unique, size=len(unique), replace=True)
    lookup = {p: np.where(bundle.participant == p)[0] for p in unique}
    return np.concatenate([lookup[p] for p in sampled])

def percentile_inference(point: pd.DataFrame, bootstrap: pd.DataFrame, out: Path) -> pd.DataFrame:
    rows = []
    for p in point.itertuples(index=False):
        data = bootstrap[(bootstrap.summary_id.eq(p.summary_id)) & (bootstrap.method.eq(p.method))]
        for metric in ("risk", "method_gain", "relative_gain", "empirical_headroom", "transport_loss", "captured_headroom"):
            values = pd.to_numeric(data[metric], errors="coerce").dropna().to_numpy(float)
            estimate = float(getattr(p, metric)) if np.isfinite(getattr(p, metric)) else math.nan
            lo, hi = np.quantile(values, [0.025, 0.975]) if len(values) else (math.nan, math.nan)
            rows.append(
                {
                    "summary_level": p.summary_level,
                    "summary_id": p.summary_id,
                    "method": p.method,
                    "metric": metric,
                    "estimate": estimate,
                    "valid_replicates": len(values),
                    "q025": float(lo),
                    "q975": float(hi),
                    "status": "PASS" if len(values) >= 0.95 * bootstrap.replicate.nunique() else "INDETERMINATE",
                }
            )
    frame = pd.DataFrame(rows)
    atomic_csv(frame, out / "04_aggregate" / "BOOTSTRAP_INFERENCE.csv")
    return frame
