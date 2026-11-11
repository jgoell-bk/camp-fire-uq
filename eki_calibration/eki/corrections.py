"""Coefficients -> correction curves -> corrected hourly weather.

theta = (a_k1, a_k2, a_k3, d_k1, d_k2, d_k3, m_k1, m_k2, m_k3): the values of the
log-speed correction a (dimensionless), the direction correction d (degrees) and
the moisture correction m (percentage points) at the knots. Each curve is piecewise
linear between knots and constant before the first and after the last. At every
weather band:

    ws   = ws0 * exp(a)
    wd   = (wd0 + d) mod 360          d is interpolated unwrapped; only wd is wrapped
    m_k  = m_k0 + w_k * m             k = 1-, 10-, 100-h

Zero coefficients give the baseline exactly. ELMFIRE then interpolates between the
hourly bands itself (wind as u/v components).

v3 adds one constant coefficient when corrections.p_eps is set (Jacob, 2026-11-10):

    p     = ln(P_EPS / reference)     one value for the whole run, appended to theta
    P_EPS = reference * exp(p)        written to the namelist (eki.case)

A member with P_EPS below corrections.p_eps.floor (or >= 1) is not admissible: like a
moisture history below its floor, it is refused, never clipped. Without p_eps in the
config, theta and everything built from it are as in v1 and v2.
"""
import numpy as np

from .common import StudyError

FAMILIES = ("a", "d", "m")
MOISTURE = ("m1", "m10", "m100")


def constants(cfg):
    """The constant (not time-varying) coefficients: ["p"] when corrections.p_eps is
    configured (v3), else none."""
    return ["p"] if cfg["corrections"].get("p_eps") else []


def theta_names(cfg):
    k = len(cfg["corrections"]["knots_h"])
    return [f"{f}_k{r}" for f in FAMILIES for r in range(1, k + 1)] + constants(cfg)


def p_eps(cfg, theta):
    """P_EPS = reference * exp(p) for a theta dict or vector (theta_names order), or
    None when P_EPS is not a coefficient."""
    if not constants(cfg):
        return None
    p = theta["p"] if isinstance(theta, dict) else np.asarray(theta, float)[theta_names(cfg).index("p")]
    return float(cfg["corrections"]["p_eps"]["reference"]) * float(np.exp(p))


def check_p_eps(cfg, theta):
    """Raises if P_EPS is a coefficient and outside [floor, 1)."""
    v = p_eps(cfg, theta)
    if v is None:
        return
    floor = float(cfg["corrections"]["p_eps"]["floor"])
    if not (np.isfinite(v) and floor <= v < 1.0):
        raise StudyError(f"inadmissible P_EPS {v:.6g}: must be in [{floor:g}, 1)")


def basis(hours, knots_h):
    """Hat functions: B[i, r] is the weight of knot r at time hours[i]. Rows sum to 1."""
    t, k = np.asarray(hours, float), np.asarray(knots_h, float)
    if k.ndim != 1 or len(k) < 1 or (np.diff(k) <= 0).any():
        raise StudyError(f"knots_h must be strictly increasing, not {knots_h}")
    B = np.zeros((len(t), len(k)))
    for i, ti in enumerate(t):
        if ti <= k[0]:
            B[i, 0] = 1.0
        elif ti >= k[-1]:
            B[i, -1] = 1.0
        else:
            j = int(np.searchsorted(k, ti, side="right")) - 1
            w = (ti - k[j]) / (k[j + 1] - k[j])
            B[i, j], B[i, j + 1] = 1.0 - w, w
    return B


def band_hours(cfg):
    """Time of each weather band in hours after t0."""
    return np.arange(cfg.nbands) * cfg["simulation"]["dt_weather_s"] / 3600.0


def check_knots(cfg):
    k = cfg["corrections"]["knots_h"]
    stop_h = cfg["simulation"]["stop_s"] / 3600.0
    basis([0.0], k)
    if k[0] < 0 or k[-1] > stop_h:
        raise StudyError(f"knots_h {k} must lie within the run, 0 .. {stop_h} h")


def curves(cfg, theta):
    """a(t), d(t), m(t) at every weather band, from a theta vector or dict."""
    names = theta_names(cfg)
    vec = np.array([theta[n] for n in names], float) if isinstance(theta, dict) else np.asarray(theta, float)
    if vec.shape != (len(names),) or not np.isfinite(vec).all():
        raise StudyError(f"theta must be {len(names)} finite numbers ({', '.join(names)})")
    k = len(cfg["corrections"]["knots_h"])
    B = basis(band_hours(cfg), cfg["corrections"]["knots_h"])
    return {f: B @ vec[i * k:(i + 1) * k] for i, f in enumerate(FAMILIES)}


def corrected_weather(cfg, baseline, theta):
    """The five hourly series written to the rasters, and the correction curves.
    Raises if the result is not admissible (see check_admissible), or if P_EPS is a
    coefficient outside its limits (check_p_eps)."""
    c = curves(cfg, theta)
    check_p_eps(cfg, theta)
    w = cfg["corrections"]["moisture_weights"]
    series = {"ws": baseline.ws_mph.to_numpy(float) * np.exp(c["a"]),
              "wd": np.mod(baseline.wd_deg.to_numpy(float) + c["d"], 360.0)}
    series["wd"] = np.where(series["wd"] >= 360.0, 0.0, series["wd"])
    for k in MOISTURE:
        series[k] = baseline[f"{k}_pct"].to_numpy(float) + float(w[k]) * c["m"]
    check_admissible(cfg, series)
    return series, c


def check_admissible(cfg, series):
    """Every band, not only the knots: finite values, wind speed >= 0, direction in
    [0, 360), and dead moisture >= corrections.moisture_min_pct (ELMFIRE would raise
    lower values to 1 % without telling)."""
    floor = float(cfg["corrections"]["moisture_min_pct"])
    problems = []
    if not all(np.isfinite(v).all() for v in series.values()):
        problems.append("non-finite value")
    if (series["ws"] < 0).any():
        problems.append("negative wind speed")
    if ((series["wd"] < 0) | (series["wd"] >= 360)).any():
        problems.append("wind direction outside [0, 360)")
    for k in MOISTURE:
        low = float(np.min(series[k]))
        if low < floor:
            problems.append(f"{k} down to {low:.3f} %, below {floor} %")
    if problems:
        raise StudyError("inadmissible weather history: " + "; ".join(problems))
