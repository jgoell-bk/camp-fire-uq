"""The regularized EKI update (Tikhonov augmentation; plan sec. 10).

    F(theta) = [G(theta); B theta],   z = [y; 0],   Gamma = diag(R, I)
    theta_j  <- theta_j + C_thetaF (C_FF + alpha Gamma)^-1 (z - F(theta_j))

C_thetaF and C_FF are ensemble covariances over the J members (divided by J - 1).
Every member uses the same z (deterministic EKI, no perturbed observations). For a
linear G, one step's new mean is the Kalman analysis with the current ensemble as
prior; repeated steps approach the minimizer of 1/2 |y - G|^2/R + 1/2 |B theta|^2.

First controlled update (eki.objective):
- G = sqrt(w) d / L, y = 0, R = 1, so Gamma = I;
- B = regularization_matrix, so 1/2 |F|^2 = Phi.
controlled_step picks alpha by doubling until the step limits in config fit.step hold.
They are numerical limits on the step; they do not promise that Phi decreases.
"""
import numpy as np
import pandas as pd

from .common import StudyError, member_id, sha256, utc_now, write_json, write_text
from .corrections import corrected_weather, p_eps, theta_names
from .forcing import load_baseline
from .objective import penalty_scale_vector


def _check(theta, G, y, R, B, alpha):
    J, p = theta.shape
    if G.shape != (J, y.size) or R.shape != y.shape:
        raise StudyError(f"shapes: theta {theta.shape}, G {G.shape}, y {y.shape}, R {R.shape}")
    if not all(np.isfinite(a).all() for a in (theta, G, y, R)):
        raise StudyError("non-finite value in theta, G, y or R; NaN never goes into an update")
    if (R <= 0).any() or alpha <= 0:
        raise StudyError("R (variances) and alpha must be positive")
    if J <= p:
        raise StudyError(f"{J} members cannot span {p} coefficients")
    if B is not None and (B.ndim != 2 or B.shape[1] != p):
        raise StudyError(f"B must have {p} columns")


def eki_update(theta, G, y, R, B=None, alpha=1.0):
    """theta (J x p), G (J x N) model outputs, y (N), R (N) variances, B (q x p) or
    None. Returns the new theta (J x p)."""
    theta, G, y, R = (np.asarray(a, float) for a in (theta, G, y, R))
    B = None if B is None else np.asarray(B, float)
    _check(theta, G, y, R, B, alpha)
    J = len(theta)
    if B is None:
        F, z, gamma = G, y, R
    else:
        F = np.hstack([G, theta @ B.T])
        z = np.concatenate([y, np.zeros(len(B))])
        gamma = np.concatenate([R, np.ones(len(B))])
    dth, dF = theta - theta.mean(axis=0), F - F.mean(axis=0)
    C_tF, C_FF = dth.T @ dF / (J - 1), dF.T @ dF / (J - 1)
    X = np.linalg.solve(C_FF + alpha * np.diag(gamma), (z - F).T)    # a linear solve, no inverse
    return theta + (C_tF @ X).T


def step_limits(theta, new, floor):
    """Standardized displacements. sd per coefficient is the current ensemble's; where
    sd < floor (collapsed) the floor is used, and the coefficient is listed."""
    sd = theta.std(axis=0, ddof=1)
    collapsed = sd < floor
    den = np.where(collapsed, floor, sd)
    mean_move = np.abs(new.mean(axis=0) - theta.mean(axis=0)) / den
    member_move = np.abs(new - theta) / den
    return mean_move, member_move.max(axis=0), collapsed


def controlled_step(cfg, theta, G, B, admissible):
    """Double alpha from alpha0 until (1) every coefficient's mean move <= max_mean_sd,
    (2) every member's move in every coefficient <= max_member_sd (current ensemble
    sd), and (3) every member is admissible. admissible(theta_row) -> bool. No
    projection is applied: admissibility is reached by damping. Returns
    (new theta, record)."""
    s = cfg["fit"]["step"]
    names = theta_names(cfg)
    floor = float(s["sd_floor"]) * penalty_scale_vector(cfg)
    y, R = np.zeros(G.shape[1]), np.ones(G.shape[1])
    alpha, tried = float(s["alpha0"]), []
    while True:
        new = eki_update(theta, G, y, R, B, alpha)
        mean_move, member_move, collapsed = step_limits(theta, new, floor)
        bad = [j for j, row in enumerate(new) if not admissible(row)]
        ok = (mean_move.max() <= s["max_mean_sd"] and member_move.max() <= s["max_member_sd"] and not bad)
        tried.append({"alpha": alpha, "max_mean_move_sd": float(mean_move.max()),
                      "max_member_move_sd": float(member_move.max()), "inadmissible_members": bad})
        if ok:
            break
        if alpha >= float(s["alpha_limit"]):
            raise StudyError(f"no alpha up to {alpha:g} meets the step limits")
        alpha *= 2.0
    record = {"alpha": alpha, "tried": tried, "limits": {k_: s[k_] for k_ in ("max_mean_sd", "max_member_sd")},
              "mean_move_sd": dict(zip(names, map(float, mean_move))),
              "max_member_move_sd": dict(zip(names, map(float, member_move))),
              "collapsed_coefficients": [n for n, c in zip(names, collapsed) if c],
              "boundary_projections": 0,
              "note": "admissibility enforced by damping (alpha), never by projection or clipping"}
    return new, record


def write_ensemble(cfg, path, iteration, theta, parent):
    """An ensemble file (run_id, iteration, member, theta, speed multipliers) and its
    .json, refusing any inadmissible member (no clipping or projection) and an
    existing file."""
    names = theta_names(cfg)
    theta = np.asarray(theta, float)
    baseline = load_baseline(cfg)
    for row in theta:
        corrected_weather(cfg, baseline, dict(zip(names, row)))   # raises if inadmissible
    k = len(cfg["corrections"]["knots_h"])
    table = pd.DataFrame({"run_id": [member_id(iteration, j + 1) for j in range(len(theta))],
                          "iteration": iteration, "member": np.arange(1, len(theta) + 1),
                          **{n: theta[:, i] for i, n in enumerate(names)},
                          **{f"speed_mult_k{r + 1}": np.exp(theta[:, r]) for r in range(k)},
                          **({"p_eps_value": [p_eps(cfg, row) for row in theta]} if "p" in names else {})})
    if path.exists():
        raise StudyError(f"{path} exists; ensemble files are never edited")
    path.parent.mkdir(parents=True, exist_ok=True)
    write_text(path, table.to_csv(index=False, float_format="%.17g"))
    write_json(path.with_suffix(".json"), {"created_utc": utc_now(), "file": path.name, "sha256": sha256(path),
                                           "config_sha256": cfg.sha256, "parent": parent})
    return path
