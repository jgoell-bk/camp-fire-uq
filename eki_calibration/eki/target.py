"""One-sided presence constraints: deadlines, the primary constraint per cell, residuals.

A fire_present report says fire was at its place by the end of its window. Its frozen
sampling cell should therefore have arrival T_i <= u_i, the deadline:
    u_i = window end (NIST start + Obs Window) + target.deadline_allowance_s,
or a tighter bound from target.deadline_overrides, where a description supports it.

selection() decides, for every report and from metadata only:
- excluded: no sampling cell, start outside the run, no window, window ending after the
  run, or a residual report;
- primary: among the comparable reports of one fixed cell, the one with the smallest
  deadline (ties to the first in file order);
- not selected: a comparable report in a cell whose primary constraint is another.
A reason is recorded for every report. Model output is never used.

deadline_residuals(): r_i = max(0, (T_i - u_i) / s_i) for a finite arrival. A cell that
is not reached by the end of the run has no residual (NaN; censored = True): it is late
by more than (stop - u_i) / s_i, but that number is never used as its residual. The
update path (fit_arrays) refuses any non-arrival until a treatment is decided.
"""
import numpy as np
import pandas as pd

from .common import BASELINE, IGNITION, NOT_REACHED, REACHED, RUN_FAILED, StudyError


def _prefix(cfg):
    return "r{:g}_".format(float(cfg["observations"]["proxy_radius_m"]))


def deadlines(cfg, obs):
    """(u_i in model seconds, basis) for every report; NaN where there is no window."""
    t = cfg["target"]
    u = obs.window_end_s + float(t["deadline_allowance_s"])
    basis = pd.Series(np.where(u.notna(), "window end" + (f" + {t['deadline_allowance_s']} s"
                                                         if t["deadline_allowance_s"] else ""), ""), index=obs.index)
    for oid, o in (t.get("deadline_overrides") or {}).items():
        i = obs.index[obs.obs_id == oid]
        if len(i) != 1:
            raise StudyError(f"deadline override for unknown {oid}")
        if not float(o["u_model_s"]) <= float(u[i[0]]):
            raise StudyError(f"deadline override for {oid} is not tighter than its window end")
        u[i[0]], basis[i[0]] = float(o["u_model_s"]), f"override: {o['reason']}"
    return u, basis


def selection(cfg, obs):
    """One row per report (obs order): deadline_s, deadline_basis, selection, reason."""
    t, s, p = cfg["target"], cfg["simulation"], _prefix(cfg)
    if t["one_per_cell"] != "smallest_deadline":
        raise StudyError(f"target.one_per_cell must be 'smallest_deadline', not {t['one_per_cell']!r}")
    u, basis = deadlines(cfg, obs)
    reason = pd.Series("", index=obs.index, dtype=object)
    rule = pd.Series("", index=obs.index, dtype=object)
    rules = [(~obs[p + "status"].isin(["direct", "proxy"]), "no_cell",
              "excluded: no sampling cell (" + obs[p + "status"] + ")"),
             (~obs.in_window, "outside_run", "excluded: starts outside the run"),
             (u.isna(), "no_window", "excluded: no observation window"),
             (u > s["stop_s"], "window_after_run", "excluded: window ends after the run"),
             ((obs.event_type == "residual") if t["exclude_residual"] else pd.Series(False, index=obs.index),
              "residual", "excluded: residual report (after the front)")]
    for cond, code, text in rules:
        new = (reason == "") & cond
        reason, rule = reason.where(~new, text), rule.where(~new, code)
    comparable = reason == ""
    c = obs[comparable].assign(u=u[comparable]).sort_values(["u", "feature_index"])
    first = c.drop_duplicates([p + "row", p + "col"])
    winner = c[[p + "row", p + "col"]].merge(first[[p + "row", p + "col", "obs_id", "u"]], how="left",
                                             on=[p + "row", p + "col"]).set_index(c.index)
    n_in_cell = c.groupby([p + "row", p + "col"]).obs_id.transform("size")
    for i in c.index:
        w = winner.loc[i]
        if w.obs_id == obs.obs_id[i]:
            rule[i], reason[i] = "primary", f"primary: smallest deadline among {n_in_cell[i]} comparable report(s) in its cell"
        elif u[i] == w.u:
            rule[i], reason[i] = "same_deadline", f"not selected: same deadline as {w.obs_id}, which comes first in the file"
        else:
            rule[i], reason[i] = "later_deadline", f"not selected: {w.obs_id} in the same cell has an earlier deadline"
    sel = np.where(rule == "primary", "primary",
                   np.where(rule.isin(["same_deadline", "later_deadline"]), "not selected", "excluded"))
    return pd.DataFrame({"obs_id": obs.obs_id, "deadline_s": u, "deadline_basis": basis,
                         "selection": sel, "selection_rule": rule, "selection_reason": reason})


def candidates(cfg, obs, partition, which=None):
    """The primary constraints of one partition (default target.partition), in obs order,
    with deadline_s. partition is a Series aligned with obs."""
    s = selection(cfg, obs)
    keep = (s.selection == "primary").to_numpy() & partition.eq(which or cfg["target"]["partition"]).to_numpy()
    return obs[keep].assign(deadline_s=s.deadline_s[keep].to_numpy())


def deadline_residuals(T, status, u, s):
    """r = max(0, (T - u) / s) where reached; NaN elsewhere. Returns (r, censored),
    censored = eligible cell not reached by the end of the run. Ignition-marker cells
    and failed runs are refused: neither should be among the constraints."""
    T, status, u, s = (np.asarray(a) for a in (T, status, u, s))
    if np.isin(status, [IGNITION, RUN_FAILED]).any():
        raise StudyError("constraint cells include the ignition cell or a failed run")
    if (np.asarray(s, float) <= 0).any():
        raise StudyError("lateness scales must be positive")
    reached = status == REACHED
    r = np.where(reached, np.maximum(0.0, (np.where(reached, T, 0.0) - u) / s), np.nan)
    return r, status == NOT_REACHED


def fit_arrays(z, obs_ids):
    """(member run IDs, theta J x p, T J x N arrival seconds) for obs_ids from an extract
    npz (BASELINE is not a member). Refuses, rather than drops or fills, when any member
    failed or any constraint cell is unreached by any member."""
    index = {o: i for i, o in enumerate(z["obs_id"])}
    cols = np.array([index[o] for o in obs_ids])
    rows = np.array([i for i, r in enumerate(z["run_id"]) if r != BASELINE])
    ids = z["run_id"][rows]
    status, arrival = z["status_primary"][np.ix_(rows, cols)], z["arrival_s_primary"][np.ix_(rows, cols)]
    failed = (status == RUN_FAILED).all(axis=1)
    if failed.any():
        raise StudyError(f"{int(failed.sum())} member runs failed ({list(ids[failed])[:5]}); "
                         "retry them; a member is never dropped from one update")
    unreached = (status != REACHED).any(axis=0)
    if unreached.any():
        raise StudyError(f"{int(unreached.sum())} of {len(cols)} constraint cells are not reached by "
                         "every member; the update needs a non-arrival treatment first")
    if not np.isfinite(arrival).all():
        raise StudyError("non-finite arrival among reached observations")
    return ids, np.asarray(z["theta"], float)[rows], arrival.astype(float)
