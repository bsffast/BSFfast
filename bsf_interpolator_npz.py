import os
import numpy as np
from scipy.interpolate import LinearNDInterpolator, interp1d

INTERP_DIR = "bsf_interpolation_data/npz"

_interpolators = {}

def load_interpolator(model):
    path = os.path.join(INTERP_DIR, f"{model}.npz")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Interpolator data for '{model}' not found at: {path}")

    npz = np.load(path)
    data = npz["data"]
    interp_type = str(npz["interp_type"])

    if interp_type == "2D":
        log_inputs = np.log(data[:, :2])
        log_sig = np.log(data[:, 2])
        interpolator = LinearNDInterpolator(log_inputs, log_sig)

    elif interp_type == "1D":
        log_x = np.log(data[:, 1])
        log_sig = np.log(data[:, 2])
        interpolator = interp1d(log_x, log_sig, kind='linear', bounds_error=False, fill_value=-np.inf)

    else:
        raise ValueError(f"Unknown interpolation type '{interp_type}' in file {path}")

    _interpolators[model] = (interpolator, interp_type)

def sigma_eff(model, m, x):
    if model not in _interpolators:
        load_interpolator(model)

    interpolator, interp_type = _interpolators[model]

    if interp_type == "2D":
        point = np.log([m, x])
        log_val = interpolator(point)
        return float(np.exp(log_val)) if log_val is not None and np.isfinite(log_val) else 0.0

    elif interp_type == "1D":
        logx = np.log(x)
        log_sig = interpolator(logx)
        if log_sig is None or not np.isfinite(log_sig):
            return 0.0
        sigma_base = np.exp(log_sig)
        return sigma_base * (m / 1.0)**-2

    else:
        raise ValueError(f"Unknown interpolation type '{interp_type}' for model '{model}'")

# --- Beispiel/Test ---
if __name__ == "__main__":
    test_cases = [
        ("stop_asPlateau", 139980.0, 398.11),
        ("QED", 100000.0, 250.0),
        ("QCDnoTrans", 1890000.0, 500.0),
        ("sbottom_noTrans_asCutoff", 300.0, 800000.0)
    ]

    for model, m, x in test_cases:
        result = sigma_eff(model, m, x)
        print(f"sigma_eff({model}, m={m}, x={x}) = {result:.3e} cm^3/s")

