import os
import numpy as np
import pickle
import re
from scipy.interpolate import LinearNDInterpolator, interp1d

def load_m_table(filepath):
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    # Capture Mathematica-like entries including *^ for exponents
    matches = re.findall(r"\{([^{}]+)\}", content)

    data = []
    for match in matches:
        # Replace Mathematica-style exponents: *^ with e
        cleaned = match.replace("*^", "e")

        # Remove backtick precision markers (e.g., 0.123`14.)
        cleaned = re.sub(r"`\d*", "", cleaned)
        cleaned = re.sub(r"\.(?=\s|$)", "", cleaned)  # remove trailing decimal dots

        try:
            numbers = [float(x) for x in cleaned.split(",") if x.strip()]
        except ValueError:
            continue  # skip malformed lines

        if 2 <= len(numbers) <= 3:
            data.append(numbers)

    return np.array(data)



def generate_interpolators():
    input_dir = "BSFfast_DataM"
    output_dir_npz = "BSFfast_DataPy/npz"
    output_dir_pkl = "BSFfast_DataPy/pkl"
    os.makedirs(output_dir_npz, exist_ok=True)
    os.makedirs(output_dir_pkl, exist_ok=True)

    models = {
        # 2D models
        "stop_asPlateau": "sigBSFeff_stop_asPlateau1.RD.m",
        "stop_asCutoff": "sigBSFeff_stop_asCutoff1.RD.m",
        "sbottom_asPlateau": "sigBSFeff_sbottom_asPlateau1.RD.m",
        "sbottom_asCutoff": "sigBSFeff_sbottom_asCutoff1.RD.m",
        "sbottom_noTrans_asPlateau": "sigBSFeff-NoTrans_sbottom_asPlateau1.RD.m",
        "sbottom_noTrans_asCutoff": "sigBSFeff-NoTrans_sbottom_asCutoff1.RD.m",

        # 1D rescalable models (mass fixed to 1, interpolated in x)
        "QCDnoTrans": "xgridSigBSFeff_QCDnoTrans_asConst0.1_m1GeV.m",
        "QED": "xgridSigBSFeff_QED_asConst0.1_m1GeV.m",
        "QEDnoTrans": "xgridSigBSFeff_QEDnoTrans_asConst0.1_m1GeV.m",
    }

    for name, filename in models.items():
        input_path = os.path.join(input_dir, filename)
        npz_path = os.path.join(output_dir_npz, name + ".npz")
        pkl_path = os.path.join(output_dir_pkl, name + ".pkl")

        data = load_m_table(input_path)
        if data.shape[1] == 3 and len(np.unique(data[:, 0])) > 1:
            # 2D case: log-log-log interpolation
            log_inputs = np.log(data[:, :2])
            log_sig = np.log(data[:, 2])
            interpolator = LinearNDInterpolator(log_inputs, log_sig)
            interp_type = "2D"
        else:
            # 1D case: fixed m=1, only x varies
            log_x = np.log(data[:, 1])
            log_sig = np.log(data[:, 2])
            interpolator = interp1d(log_x, log_sig, kind='linear', bounds_error=False, fill_value=-np.inf)
            interp_type = "1D"

        # Save npz
        np.savez_compressed(npz_path, data=data, interp_type=interp_type)

        # Save pkl
        with open(pkl_path, "wb") as f:
            pickle.dump((interpolator, interp_type), f)

        print(f"Generated interpolator for {name} ({interp_type}) with {len(data)} points.")

if __name__ == "__main__":
    generate_interpolators()

