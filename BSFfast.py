import numpy as np
import warnings
from pathlib import Path
from collections import defaultdict

from scipy.interpolate import interp1d, LinearNDInterpolator

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

DATA_DIR = Path(__file__).parent / "BSFfast_Data"

DEFAULT_SCHEME = "cutoff"
ALPHA0_DEFAULT = 0.1   # reference coupling for rescaled tables
M0_DEFAULT = None      # read from table (unique mass)

C_POT_SCALE = 2.0   # q = m * sqrt(C_POT_SCALE / x)

# --- Unitarity warning (Sec. 4.3; fit used for dQCD rescaling) ---
WARN_UNI = True
# Fit: log10(alpha) = A + B log10(1/v) + C log10(rho), with rho = sigma_BSF/sigma_uni
UNI_FIT_A = -0.166
UNI_FIT_B = -0.251
UNI_FIT_C = -0.250

# Avoid spamming in Boltzmann loops: warn once per model per threshold
_warned = set()


# ------------------------------------------------------------
# Physics helpers
# ------------------------------------------------------------

def _potential_scale(m, x):
    """
    Typical momentum / potential scale used for approximate running-coupling results,
    see Eq. (39) and Sec. 3 of the paper.
    """
    return m * np.sqrt(C_POT_SCALE / x)

def _thermal_v_from_x(x):
    """
    Characteristic relative velocity used for the unitarity-warning estimate.
    """
    return np.sqrt(6.0 / x)

def _unitarity_ratio_estimate(alpha, v):
    """
    Estimate r = sigma_BSF / sigma_uni from the Sec. 4.3 fit.
    The fit is written with a negative coefficient in front of log10(rho),
    so the fitted rho corresponds to ~ sigma_uni / sigma_BSF, hence r = 1/rho.
    """
    if alpha <= 0.0 or v <= 0.0:
        return 0.0

    log10_alpha = np.log10(alpha)
    log10_1_over_v = np.log10(1.0 / v)

    # log10(alpha) = A + B log10(1/v) + C log10(rho)  with A,B,C negative (as in paper)
    log10_rho = (log10_alpha - UNI_FIT_A - UNI_FIT_B * log10_1_over_v) / UNI_FIT_C
    rho = 10.0 ** log10_rho
    
    # invert to get the physically relevant ratio r = sigma_BSF / sigma_uni
    return 1.0 / rho


def _maybe_warn_unitarity(model, x, alpha_eff):
    """
    Emit warnings if the (estimated) cross section exceeds 10% or 100% of unitarity.
    Implemented for dQCD-* rescaled results (Sec. 4.3 fit).
    """
    if not WARN_UNI:
        return
    if not model.startswith("dQCD-"):
        return

    v = _thermal_v_from_x(x)
    r = _unitarity_ratio_estimate(alpha_eff, v)

    # warn once per model and level
    if r > 1.0 and (model, "100") not in _warned:
        _warned.add((model, "100"))
        warnings.warn(
            "Warning: cross section estimated to violate unitarity limit:",
            RuntimeWarning
        )
    elif r > 0.1 and (model, "10") not in _warned:
        _warned.add((model, "10"))
        warnings.warn(
            "Warning: cross section estimated to reach >10% of unitarity limit:",
            RuntimeWarning
        )


# ------------------------------------------------------------
# Internal registry (auto-filled at import)
# ------------------------------------------------------------

_registry = defaultdict(dict)


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def _load_csv(path):
    return np.loadtxt(path, delimiter=",")


def _unique(arr):
    return np.unique(arr)


def _build_rescaled_interpolator(data):
    """
    data columns: m0, x, sigma_ref
    """
    m_vals = _unique(data[:, 0])
    if len(m_vals) != 1:
        raise ValueError("Rescaling table must contain exactly one mass")

    m0 = float(m_vals[0])
    x = data[:, 1]
    sig = data[:, 2]

    # interpolate in log-log for stability
    interp = interp1d(
        np.log(x),
        np.log(sig),
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
    )

    def xs(x, m, alpha, alpha0=ALPHA0_DEFAULT):
        xprime = x * (alpha / alpha0) ** 2
        sig_ref = np.exp(interp(np.log(xprime)))
        return sig_ref * (alpha / alpha0) ** 2 * (m0 / m) ** 2

    return xs


def _build_2d_interpolator(data):
    """
    data columns: m, x, sigma
    """
    m = data[:, 0]
    x = data[:, 1]
    sig = data[:, 2]

    # log-log interpolation
    pts = np.column_stack([np.log(m), np.log(x)])
    vals = np.log(sig)

    interp = LinearNDInterpolator(pts, vals, fill_value=np.nan)

    def xs(xq, mq):
        y = interp(np.log(mq), np.log(xq))
        if np.isnan(y):
            # optionally: raise instead of returning nan
            return np.nan
        return float(np.exp(y))

    return xs


def _parse_filename(name):
    """
    Extract model and scheme from filename.
    Examples:
      QCD-SU_cutoff.csv -> model=QCD-SU, scheme=cutoff
      dQCD-S.csv        -> model=dQCD-S, scheme=None
    """
    stem = name.replace(".csv", "")
    if "_" in stem:
        model, scheme = stem.rsplit("_", 1)
        if scheme in ("cutoff", "plateau"):
            return model, scheme
    return stem, None


# ------------------------------------------------------------
# Load all tables at import
# ------------------------------------------------------------

for path in DATA_DIR.glob("*.csv"):
    model, scheme = _parse_filename(path.name)
    data = _load_csv(path)

    masses = _unique(data[:, 0])

    if len(masses) == 1:
        # rescaled table
        xs_func = _build_rescaled_interpolator(data)
        _registry[model][scheme or DEFAULT_SCHEME] = {
            "mode": "rescaled",
            "xs": xs_func,
        }
    else:
        # true 2D grid
        xs_func = _build_2d_interpolator(data)
        _registry[model][scheme or DEFAULT_SCHEME] = {
            "mode": "2D",
            "xs": xs_func,
        }


# ------------------------------------------------------------
# Public interface
# ------------------------------------------------------------

def fastXS(model, x, m, alpha=None):
    """
    Evaluate effective BSF cross section.

    Parameters
    ----------
    model : str
        Model name, e.g. 'dQCD-S', 'QCD-SU'
    x : float
        Temperature parameter, x = m / T
    m : float
        Particle mass m [GeV]
    alpha : None, float, or callable(q), or str ('cutoff'/'plateau')
        - None: use reference alpha (rescaled models) or ignore (2D models)
        - float: constant coupling for rescaled models
        - callable: running running alpha(q) for approximate rescaling
        - str: scheme override (for QCD-* grids), e.g. "plateau"

    Returns
    -------
    float
        Effective thermally averaged BSF cross section [GeV^(-2)]
    """

    if model not in _registry:
        raise ValueError(f"Unknown model '{model}'")

    # scheme via alpha-string hack
    scheme = DEFAULT_SCHEME
    if isinstance(alpha, str):
        scheme = alpha if alpha in ("cutoff", "plateau") else DEFAULT_SCHEME
        alpha = None

    entry = _registry[model].get(scheme)
    if entry is None:
        raise ValueError(f"No scheme '{scheme}' for model '{model}'")

    mode = entry["mode"]
    xs_func = entry["xs"]

    if mode == "2D":
        return float(xs_func(x, m))

    # rescaling branch
    if alpha is None:
        alpha_eff = ALPHA0_DEFAULT
    elif callable(alpha):
        # user provides alpha(q); q choice is internal
        q = _potential_scale(m, x)
        alpha_eff = float(alpha(q))
    else:
        alpha_eff = float(alpha)

    # Unitarity warnings only for rescaled results (Sec. 4.3)
    _maybe_warn_unitarity(model, x, alpha_eff)

    return float(xs_func(x, m, alpha_eff))


