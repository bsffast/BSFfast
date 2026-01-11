# example_use.py

import numpy as np
import warnings

# Make warnings one-line and yellow (for nicer terminal output)
warnings.formatwarning = lambda msg, *args, **kwargs: f"\033[1;33m{msg}\033[0m\n"

from BSFfast import fastXS


def header(title):
    print("\n" + "-" * 60)
    print(title)
    print("-" * 60)


def main():

    # Example parameters
    x = 1.5e5          # m / T
    m = 2.5e6          # GeV

    # ------------------------------------------------------------
    # 1) Rescaled dQCD model with constant alpha
    # ------------------------------------------------------------
    alpha_dQCD = 0.15
    header(f"1) dQCD-S (constant alpha = {alpha_dQCD})")
    xs_dqcd = fastXS("dQCD-S", x, m, alpha_dQCD)
    print(f"<sigma v> = {xs_dqcd:.3e} GeV^(-2)")

    # ------------------------------------------------------------
    # 2) Rescaled dQCD model with approximate running alpha
    # ------------------------------------------------------------
    def alpha_running(q):
        # simple 1-loop QCD running (toy example)
        return 4.0 * np.pi / (7.67 * np.log((q*q + 1.0) / (0.2*0.2)))

    header("2) dQCD-S (approximate running alpha; toy 1-loop running)")
    xs_dqcd_run = fastXS("dQCD-S", x, m, alpha_running)
    print(f"<sigma v> = {xs_dqcd_run:.3e} GeV^(-2)")

    # ------------------------------------------------------------
    # 3) Rescaled dQED model with constant alpha
    # ------------------------------------------------------------
    alpha_dQED = 0.01
    header(f"3) dQED-S (constant alpha = {alpha_dQED})")
    xs_dqed = fastXS("dQED-S", x, m, alpha_dQED)
    print(f"<sigma v> = {xs_dqed:.3e} GeV^(-2)")

    # ------------------------------------------------------------
    # 4) SM QCD model (2D grid, default cutoff scheme)
    # ------------------------------------------------------------
    header("4) QCD-SU (cutoff scheme)")
    xs_qcd_cut = fastXS("QCD-SU", x, m)
    print(f"<sigma v> = {xs_qcd_cut:.3e} GeV^(-2)")

    # ------------------------------------------------------------
    # 5) SM QCD model (2D grid, plateau scheme)
    # ------------------------------------------------------------
    header("5) QCD-SU (plateau scheme)")
    xs_qcd_plat = fastXS("QCD-SU", x, m, "plateau")
    print(f"<sigma v> = {xs_qcd_plat:.3e} GeV^(-2)")

    # ------------------------------------------------------------
    # 5) SM QED model (default coupling alpha=1/128.9)
    # ------------------------------------------------------------
    
    header("6) QED-S (plateau alpha_em)")
    xs_qeds = fastXS("QED-S", x, m)
    print(f"<sigma v> = {xs_qeds:.3e} GeV^(-2)\n")
    

if __name__ == "__main__":
    main()





