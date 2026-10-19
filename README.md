# VQE electronic-structure simulation of MoS₂ quantum dots

Code and data behind the NanoAfrica 2026 poster (Poster #31, 21 Oct 2026):

> **Variational Quantum Eigensolver Approach to Electronic Structure Simulation of Molybdenum Disulfide Quantum Dots: Towards Quantum-Accelerated Nanomaterial Design**
> C.B. Moyo, A. Ogunlaja, Z.R. Tshentu, Department of Chemistry, Nelson Mandela University, Gqeberha, South Africa

This pilot project of the **NMU Quantum Computing Society** is, to our knowledge, the first application of quantum computing to chemistry at Nelson Mandela University.

## Results at a glance

| | MoS₂ | Mo₃S₆ |
|---|---|---|
| CAS(4e,4o) natural occupations | 1.98 / 1.89 / 0.11 / 0.02 | 1.65 / 1.44 / 0.56 / 0.35 |
| HOMO–LUMO gap, PBE → PBE0 (eV) | 0.11 → 2.26 | 0.58 → 1.87 |
| Hartree–Fock error vs CASCI (mHa) | 34.8 | 123.7 |
| UCCSD error (mHa) / two-qubit gates on Heron | 0.04 / 540 | 1.71 / 621 |
| GateFabric, 2 layers (mHa) / two-qubit gates | 0.08 / 210 | 32.7 / 207 |
| **ADAPT-2 compact circuit** (mHa) / two-qubit gates | **0.70 / 6** | 40.8 / 6 |
| Heron noise-model simulation: raw → TREX → ZNE (mHa) | 23.6 → 7.4 → 0.6 ± 0.7 | 56.9 → 44.4 → 41.6 ± 2.7 |
| **Real QPU, ibm_kingston** (156-qubit IBM Heron r2): raw → TREX + ZNE (mHa) | **123.1 ± 7.3 → 25.8 ± 14.4** | not run |

Errors are total-energy errors relative to CASCI in the same active space. Chemical accuracy is 1.6 mHa (≈ 43 meV ≈ 1 kcal/mol); 1 mHa = 27.2 meV.

**Hardware details.** IBM Quantum job IDs: raw `dbb2shs24mrc73867ne0`, mitigated `dbb2si7i0auc7392rgb0`. 2,048 shots per measurement basis, 8 twirling randomisations, XY4 dynamical decoupling, ZNE at noise factors λ = 1, 3, 5. Result file: `MoS2_cas44_A2_run1_hardware_ibm_kingston.json`.

**What this is not.** The clusters are small stoichiometric models, not full quantum dots. The hardware result is not yet chemically accurate, and no quantum-advantage claim is made: every circuit here is classically simulable.

## Method

1. **Model systems.** (MoS₂)ₙ clusters, n = 1 and 3, optimised at PBE/def2-SVP with the def2 ECP on Mo (PySCF + PyBerny).
2. **Classical references.** PBE, PBE0 and HSE06 gaps; RHF; CASSCF(4e,4o) and (6e,6o); CASCI as the exact reference in each active space.
3. **Qubit Hamiltonian.** CASCI natural orbitals, Jordan–Wigner mapping in PennyLane. (4e,4o) gives 8 qubits, 193 Pauli terms and 44 qubit-wise-commuting measurement groups. Exact diagonalisation reproduces CASCI to better than 10⁻⁹ Ha.
4. **Ansätze.** UCCSD and GateFabric as benchmarks; ADAPT-VQE grown one excitation at a time. Disjoint pair excitations are prepared exactly with one R_y and 3 CNOTs each, so the MoS₂ ADAPT-2 circuit needs only 6 two-qubit gates.
5. **Error mitigation.** TREX-style readout correction, then zero-noise extrapolation by global unitary folding. Parameters are optimised classically, so the QPU only evaluates the energy (minutes of QPU time, not hours).

## Repository layout

| File | Purpose |
|---|---|
| `01_classical.py` | Geometry, DFT, CASSCF/CASCI, correlation diagnostics; writes `*_opt.xyz`, `*_classical.json`, `*_cas44_mo.npy` |
| `01b_finish.py` | Time-boxed finisher used for Mo₃S₆ when CASSCF(6e,6o) stalled |
| `02_vqe_sim.py` | Qubit Hamiltonian, CASCI check, UCCSD / GateFabric VQE |
| `02b_adapt.py` | ADAPT-VQE accuracy versus circuit size |
| `03_noisy_zne.py` | Shot-based simulation with IBM fake-backend noise models, TREX, ZNE, error bars |
| `04_hardware.py` | Real IBM Quantum run (EstimatorV2) |
| `get_results.py`, `check_job.py` | Retrieve jobs and check QPU usage |
| `save_ibm_account.py` | One-time IBM account setup (prompts for your key; nothing is stored in the repo) |
| `gatecount.py` | Two-qubit gate counts on a Heron layout |
| `05_figures.py` | All poster figures → `fig/` |
| `*.json`, `*.xyz` | Results and optimised geometries |
| `PILOT_PLAN.md` | The original 24-hour working plan, kept for transparency |

`MoS2_cas44_A2_hardware_ibm_kingston.json` records a second submission that was cancelled. The hardware numbers above come from the `run1` file.

## Reproduce

```bash
pip install -r requirements.txt
python 01_classical.py MoS2          # also writes MoS2_cas44_mo.npy, needed by step 2
python 02_vqe_sim.py MoS2 4 4
python 02b_adapt.py MoS2 4 4 6
python 03_noisy_zne.py MoS2 4 4 A2 8192 5 FakeMarrakesh
python 05_figures.py
```

Repeat with `Mo3S6` for the larger cluster. The CASSCF orbital files (`*_mo.npy`) are not committed, so run step 1 first.

For real hardware you need your own IBM Quantum account: run `python save_ibm_account.py` once, then `python 04_hardware.py MoS2 4 4 A2 ibm_kingston`. **Never commit your API key**; it is stored outside the repository, in your user folder.

## Cite

See `CITATION.cff`, or cite the NanoAfrica 2026 poster above.

## Contact and acknowledgements

Cyprian Moyo · cyprian.moyo@nithecs.ac.za · [LinkedIn](https://www.linkedin.com/in/cyprian-moyo-compchem-moleculardesign-cpu-gpu-qpu)

Nelson Mandela University · NITheCS · NMU Quantum Computing Society · African Quantum Consortium · IBM Quantum · CHPC.
