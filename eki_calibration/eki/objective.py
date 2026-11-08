"""The fitting objective of the first controlled update (config `fit`).

    Phi(theta) = 1/2 sum_i w_i (d_i / L)^2 + R_amp + R_smooth

d_i is the distance from constraint i's fixed sampling-cell centre to the cumulative
burned area at its deadline u_i (eki.pilot.distance_to_burned); it is 0 when the
cell burned by u_i.

w_i = 1 / (K n_b): every occupied fit block b (K of them, n_b constraints each)
carries the same total weight, and the weights sum to 1. This is an area-balanced
objective, not a claim of independence.

The regularization uses q_g(t) = c_g(t) / scale_g, the correction curve of family g
(log-speed, direction, moisture) divided by its penalty scale:
    R_amp    = amplitude_weight * 1/2 sum_g (1 / T) integral_0^T q_g(t)^2 dt,   T = stop_s
    R_smooth = 1/2 sum_g sum_adjacent (ref_h / dt_h) (q_next - q_prev)^2
The integral is exact for the piecewise-linear curve: (h/3)(q_a^2 + q_a q_b + q_b^2)
per segment, plus the constant part after the last knot. regularization_matrix
returns B with 1/2 |B theta|^2 = R_amp + R_smooth exactly, so the EKI outputs are
[sqrt(w) d / L ; B theta] with a zero target and Gamma = I.
"""
import numpy as np

from .common import StudyError
from .corrections import FAMILIES, theta_names


def block_weights(blocks):
    """w_i = 1 / (K n_b) for a sequence of block labels (one per constraint)."""
    blocks = np.asarray(blocks)
    labels, inverse, counts = np.unique(blocks, return_inverse=True, return_counts=True)
    return 1.0 / (len(labels) * counts[inverse])


def mass_matrix(knots_h, end_h):
    """M with q^T M q = (1/T) integral_0^T q(t)^2 dt for the hat-function curve on
    knots_h (constant outside them), T = end_h. Requires knots_h[0] >= 0."""
    k = np.asarray(knots_h, float)
    if k[0] < 0 or k[-1] > end_h:
        raise StudyError("knots must lie within the run")
    M = np.zeros((len(k), len(k)))
    M[0, 0] += k[0]                                   # constant before the first knot
    for r in range(len(k) - 1):
        h = k[r + 1] - k[r]
        M[r:r + 2, r:r + 2] += h / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])
    M[-1, -1] += end_h - k[-1]                         # constant after the last knot
    return M / end_h


def regularization_matrix(cfg):
    """B (rows x 9): amplitude rows (per family, sqrt(weight) L_M^T q_g with
    L_M L_M^T = M) and then smoothness rows (per family and adjacent pair,
    sqrt(ref_h / dt_h) (q_next - q_prev))."""
    f = cfg["fit"]
    names = theta_names(cfg)
    knots = np.asarray(cfg["corrections"]["knots_h"], float)
    k = len(knots)
    scales = {g: float(f["penalty_scales"][g]) for g in FAMILIES}
    if min(scales.values()) <= 0:
        raise StudyError("penalty scales must be positive")
    Lm = np.linalg.cholesky(mass_matrix(knots, cfg["simulation"]["stop_s"] / 3600.0))
    rows = []
    w = float(f["amplitude_weight"])
    if w > 0:
        for gi, g in enumerate(FAMILIES):
            block = np.zeros((k, len(names)))
            block[:, gi * k:(gi + 1) * k] = np.sqrt(w) * Lm.T / scales[g]
            rows.extend(block)
    ref = float(f["smoothness_reference_h"])
    if ref > 0:
        for gi, g in enumerate(FAMILIES):
            for r in range(k - 1):
                row = np.zeros(len(names))
                c = np.sqrt(ref / (knots[r + 1] - knots[r])) / scales[g]
                row[gi * k + r], row[gi * k + r + 1] = -c, c
                rows.append(row)
    return np.array(rows).reshape(-1, len(names))


def regularization_terms(cfg, theta, n_grid=24001):
    """(R_amp, R_smooth) per member, computed directly from the curves (amplitude by
    fine trapezoid quadrature), independently of regularization_matrix."""
    from .corrections import basis
    f = cfg["fit"]
    theta = np.atleast_2d(np.asarray(theta, float))
    knots = np.asarray(cfg["corrections"]["knots_h"], float)
    k = len(knots)
    end = cfg["simulation"]["stop_s"] / 3600.0
    t = np.linspace(0.0, end, n_grid)
    Bt = basis(t, knots)
    amp = np.zeros(len(theta))
    smooth = np.zeros(len(theta))
    for gi, g in enumerate(FAMILIES):
        q = theta[:, gi * k:(gi + 1) * k] / float(f["penalty_scales"][g])
        amp += 0.5 * np.trapezoid((q @ Bt.T) ** 2, t, axis=1) / end
        smooth += 0.5 * (float(f["smoothness_reference_h"]) / np.diff(knots) * np.diff(q, axis=1) ** 2).sum(axis=1)
    return float(f["amplitude_weight"]) * amp, smooth


def distance_outputs(d_m, w, L_m):
    """sqrt(w_i) d_i / L (members x constraints). Distances must be finite and >= 0."""
    d = np.asarray(d_m, float)
    if not np.isfinite(d).all() or (d < 0).any():
        raise StudyError("distances must be finite and non-negative")
    if not np.isclose(np.sum(w), 1.0):
        raise StudyError("weights must sum to 1")
    return np.sqrt(np.asarray(w, float)) * d / float(L_m)


def objective(cfg, theta, d_m, w):
    """Phi and its parts per member: distance, amplitude, smoothness, total."""
    theta = np.atleast_2d(np.asarray(theta, float))
    F = distance_outputs(np.atleast_2d(d_m), w, cfg["fit"]["length_scale_m"])
    B = regularization_matrix(cfg)
    k = len(cfg["corrections"]["knots_h"])
    n_amp = 3 * k if float(cfg["fit"]["amplitude_weight"]) > 0 else 0
    Bt = theta @ B.T
    parts = {"distance": 0.5 * (F ** 2).sum(axis=1), "amplitude": 0.5 * (Bt[:, :n_amp] ** 2).sum(axis=1),
             "smoothness": 0.5 * (Bt[:, n_amp:] ** 2).sum(axis=1)}
    parts["total"] = parts["distance"] + parts["amplitude"] + parts["smoothness"]
    return parts
