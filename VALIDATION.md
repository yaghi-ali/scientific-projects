# Validation — 2026-10-04

Environment: Windows, Python 3.12, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2.

- `python -m unittest -v test_projects`: **7 tests passed**. Checks cover telecom round trips and SNR, fibre-mode families, hologram energy/error, turbulence scaling and PSF normalization, exact Grover probabilities, counting state transitions, and input validation.
- `telecom_report.ipynb`: executed end-to-end with nbclient on the generated synthetic trace; saved outputs are from this run. The notebook needs no external measurement file.
- EDFA: reducing dz from 0.001 m to 0.0005 m changes the maximum gain by 0.004041 dB or less; pump and signal remain positive. This is a numerical convergence check, not experimental validation.
- `grover_qiskit.ipynb`: execution attempted but **not validated**: Windows Application Control blocked the Qiskit native extension. Outputs remain empty. The separate NumPy Grover implementation passes analytical checks. Aer sampling uses seed 42.

No hardware, client deployment, laboratory accuracy or full-vector electromagnetic model is validated by these checks. Synthetic datasets and simplified assumptions are documented in the README.
