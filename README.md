# Scientific projects — Ali Yaghi

Reproducible studies in photonics, signal processing, numerical modelling and embedded-vision counting. French academic projects, with an English overview for readers and collaborators.

These are cleaned-up portfolio editions, not claims of commercial deployments. Some modules are reconstructed from earlier notebooks; the scope and limits are explicit below.

## Run locally

Python 3.12 recommended. From this repository:
```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest -v test_projects
python telecom_fdm.py
python fiber_modes.py
python holography_gs.py
python adaptive_optics.py
python grover_search.py
python tracking_counter.py
python edfa_study.py
```

For notebooks, install `requirements-notebooks.txt` and select that Python environment in Jupyter or VS Code.

## Projects / Projets

| Study | Code | Method and scope |
|---|---|---|
| Optical telecommunications | [telecom_fdm.py](telecom_fdm.py), [report](telecom_report.ipynb) | Three OOK subcarriers, I/Q integration, UTF-8 decoding, AWGN and measured BER. Synthetic default; optional local trace is not distributed. |
| Step-index fibre modes | [fiber_modes.py](fiber_modes.py) | Scalar LP boundary equation solved with Bessel functions and bracketing. Infinite cladding, weak guidance; not a full-vector Tidy3D result. |
| Phase-only holography | [holography_gs.py](holography_gs.py) | Gerchberg–Saxton phase retrieval, generated target, energy-normalized Fourier propagation. Rebuilt from the nanotechnology project. |
| Atmospheric/adaptive optics | [adaptive_optics.py](adaptive_optics.py) | Von Kármán phase screen, pupil and normalized PSF. Rebuilt from the optical wavefront notebook. Finite periodic grid without subharmonics underrepresents low frequencies. |
| Grover search | [grover_search.py](grover_search.py), [Qiskit notebook](grover_qiskit.ipynb) | Two-/three-qubit examples, exact state-vector probabilities and circuit sampling. Simulator only. |
| Embedded vision counting | [tracking_counter.py](tracking_counter.py) | Standalone reconstruction of IN/OUT counting on already assigned track IDs. Hysteresis, reversed crossings and track expiry. Not the original STM32 firmware or an object detector. |
| Erbium fibre amplifier | [edfa_study.py](edfa_study.py) | Coupled signal/pump Euler model and HTML plots from a TD. ASE is post-processed, not coupled into saturation. NF follows the exercise expression; do not treat it as a general amplifier noise model. |

### Fibre interpretation

The study uses radius 5 µm, wavelength 1.55 µm and NA 0.26. NA alone does not determine effective indices: **n_core = 1.45 is an explicit added assumption**. V is about 5.27. Count LP families separately from azimuthal and polarization degeneracy. V²/2 is a large-V estimate, not an exact count of LP families. Values are recalculated by the script; old unverified numbers are not reused.

### Embedded-vision provenance

The project context involved STM32N6 detection/tracking and ESP32 supervision. The public module accepts tracked normalized positions and emits JSON events. The original hardware integration, detector weights and recordings are not part of this reconstruction; no field-accuracy result is claimed here. ST software remains available under its own terms from [STMicroelectronics](https://github.com/STMicroelectronics).

## Other projects — existing repositories, no duplicates

- [Beam propagation / FFT-BPM](https://github.com/yaghi-ali/beam-propagation-ml)
- [Automated laser attenuator](https://github.com/yaghi-ali/Automatized-variable-laser-power-attenuator)
- [Signal processing and machine learning](https://github.com/yaghi-ali/ml-neural-network-project)

## Reproducibility and contribution rules

See [VALIDATION.md](VALIDATION.md) for actual checks and [CONTRIBUTING.md](CONTRIBUTING.md) for conventions. Generated data are examples, not laboratory measurements. Physical units, assumptions and limitations must remain documented. Each file retains its applicable attribution; no blanket relicensing of upstream code is implied.

References: [SciPy special functions](https://docs.scipy.org/doc/scipy/reference/special.html), [Brent root solver](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html), [Optiwave fibre-mode guide](https://optiwave.com/optifiber-manuals/optical-fiber-fiber-modes/).

Ali Yaghi · [OptiIA Consulting](https://www.optiia-consulting.fr) · contact@optiia-consulting.fr
