"""Which fuel cells can carry an ELMFIRE arrival time, from the fixed landscape only."""
import numpy as np

URBAN = 91   # nonburnable for surface spread, but burns through the building model
FUEL_NAMES = {91: "urban", 92: "snow/ice", 93: "agriculture", 98: "water", 99: "barren",
              256: "NB in the case's fuel table"}


def nonburnable(code):
    """ELMFIRE's ISNONBURNABLE for the Rothermel model (elmfire_init.f90)."""
    code = np.asarray(code)
    return ((code >= 90) & (code <= 100)) | (code == 256) | (code <= 0) | (code > 303)


def eligible(fbfm, bldg_area):
    """A surface-burnable fuel model, or urban fuel 91 with building area > 0 (burned by
    the building model). In the Camp case every fuel-91 cell has building area."""
    fbfm = np.asarray(fbfm)
    return ~nonburnable(fbfm) | ((fbfm == URBAN) & (np.asarray(bldg_area) > 0))


def fuel_class(code, review_fuels):
    code = int(code)
    if not nonburnable(code):
        return "wildland"
    if code == URBAN:
        return "urban"
    if code == 256:
        return "nonburnable_256"
    if code in review_fuels:
        return "nonburnable_review"
    return "nonburnable_other"
