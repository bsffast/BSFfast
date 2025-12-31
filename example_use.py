# example_use.py

import numpy as np
from BSFfast import fastXS


def main():

    # Example parameters
    x = 100.0          # m / T
    m = 1000.0         # GeV

    # ------------------------------------------------------------
    # 1) Rescaled model with constant alpha
    # ------------------------------------------------------------
    alpha = 0.2
    xs_dqcd = fastXS("dQCD-S", x, m, alpha)
    print(f"dQCD-S (alpha={alpha}):  <sigma v> = {xs_dqcd:.3e}")

    # ------------------------------------------------------------
    # 2) Rescaled model with approximate running alpha
    # ------------------------------------------------------------
    def alpha_running(q):
        # simple 1-loop QCD running (toy example)
        return 4.0 * np.pi / (7.67 * np.log((q*q + 1.0) / (0.2*0.2)))
        
    xs_dqcd_run = fastXS("dQCD-S", x, m, alpha_running)
    print(f"dQCD-S (running alpha): <sigma v> = {xs_dqcd_run:.3e}")

    # ------------------------------------------------------------
    # 3) QED-like model (alpha fixed internally)
    # ------------------------------------------------------------
    xs_qed = fastXS("dQED-S", x, m)
    print(f"dQED-S (alpha_em):      <sigma v> = {xs_qed:.3e}")

    # ------------------------------------------------------------
    # 4) SM QCD model (2D grid, cutoff scheme)
    # ------------------------------------------------------------
    xs_qcd_cut = fastXS("QCD-SU", x, m)
    print(f"QCD-SU (cutoff):        <sigma v> = {xs_qcd_cut:.3e}")

    # ------------------------------------------------------------
    # 5) SM QCD model (2D grid, plateau scheme)
    # ------------------------------------------------------------
    xs_qcd_plat = fastXS("QCD-SU", x, m, "plateau")
    print(f"QCD-SU (plateau):       <sigma v> = {xs_qcd_plat:.3e}")


if __name__ == "__main__":
    main()
