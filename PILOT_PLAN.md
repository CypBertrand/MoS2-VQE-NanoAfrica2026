# MoS₂ quantum-dot VQE

NanoAfrica 2026 poster: *Variational Quantum Eigensolver Approach to Electronic Structure Simulation of MoS₂ Quantum Dots*
C.B. Moyo, A. Ogunlaja, Z.R. Tshentu (Nelson Mandela University) · NMU Quantum Computing Society pilot project

Every number on the poster draft was produced by these scripts. 
---

## 1. Methodology 

**Model systems.** Two stoichiometric (MoS₂)ₙ clusters, n = 1 (MoS₂) and n = 3 (Mo₃S₆), in a trigonal-prismatic Mo₃ motif cut from the 2H monolayer. Both were optimised at PBE/def2-SVP with the def2 effective core potential on Mo, density fitting and second-order SCF (PySCF + PyBerny). These are the smallest units of a (MoS₂)ₙ size series. Larger, H-passivated dots are listed under future work.

**Classical references.** At the PBE geometry:
- PBE, PBE0 and HSE06 single points, to show how much the HOMO–LUMO gap depends on the functional
- an unrestricted PBE0 triplet, for the singlet–triplet gap
- RHF
- CASSCF(4e,4o) and CASSCF(6e,6o)
- CASCI in the same active spaces (the *exact* reference that VQE is scored against)

Strong correlation is quantified by natural-orbital occupations and the Mo-4d character of the active orbitals.

**Active space and qubit Hamiltonian.** The CASSCF active orbitals were rotated to CASCI natural orbitals. This folds orbital relaxation into the basis, so the quantum circuit only needs double excitations. One- and two-electron integrals and the frozen-core energy came from PySCF. The Hamiltonian was built in second quantisation, mapped with Jordan–Wigner in PennyLane, and checked by exact diagonalisation against CASCI (agreement better than 10⁻⁹ Ha). (4e,4o) gives 8 qubits, 193 Pauli terms and 44 qubit-wise-commuting measurement groups. (6e,6o) gives 12 qubits.

**Ansätze (noiseless statevector, L-BFGS-B, adjoint gradients).**
- UCCSD, as the accuracy benchmark
- GateFabric, a particle-conserving hardware-efficient ansatz with 1–2 layers
- ADAPT-VQE: excitations are added greedily by energy gradient, giving an accuracy-vs-circuit-size curve

**Hardware-ready circuit.** When ADAPT selects double excitations on disjoint orbital pairs, the state is prepared *exactly* with 1 R_y + 3 CNOT per excitation instead of the textbook ~14 CNOT. MoS₂ (4e,4o) with 2 excitations needs 6 two-qubit gates after transpilation. All parameters are optimised classically, so the QPU only *evaluates* the energy at fixed parameters. That makes each system cost minutes of QPU time, not hours.

**Noise and error mitigation (device-calibrated simulation).** Circuits were transpiled at optimisation level 3 for an Eagle device (ibm_brisbane calibration snapshot) and a Heron device (ibm_marrakesh snapshot). They were then sampled shot by shot in Qiskit Aer with the full device noise model: 8,192 shots per group, 5 independent repeats for error bars. Mitigation was applied in this order:
1. **TREX-style readout correction.** Per-qubit factors (1 − p₀₁ − p₁₀) are calibrated from |0⟩ and |1⟩ preparations.
2. **Zero-noise extrapolation.** Global unitary folding at λ = 1, 3, 5, then linear and Richardson extrapolation to λ = 0.

The Hartree–Fock circuit (X gates only) is measured alongside as a calibration control.

**Hardware.** `04_hardware.py` submits the same two circuits through Qiskit Runtime `EstimatorV2` in a `Batch`. Settings: TREX, ZNE (1, 3, 5), gate and measurement twirling, XY4 dynamical decoupling, plus an unmitigated (resilience 0) run for comparison.

**Figure of merit.** Error relative to CASCI in the same active space. The target is chemical accuracy (1.6 mHa).

## 2. Step-by-step: from zero to finished poster

### Files
| File | Purpose |
|---|---|
| `01_classical.py` | geometry, DFT, CASSCF/CASCI, correlation diagnostics |
| `01b_finish.py` | time-boxed finisher (DFT gaps, triplet, CAS(4e,4o)) when CASSCF(6e,6o) stalls |
| `gatecount.py` | two-qubit gate counts of UCCSD/GateFabric on a Heron layout |
| `02_vqe_sim.py` | qubit Hamiltonian, convention check, UCCSD / GateFabric VQE |
| `02b_adapt.py` | ADAPT-VQE accuracy-vs-size curve |
| `03_noisy_zne.py` | device-noise sampling, TREX, ZNE, error bars |
| `04_hardware.py` | real IBM Quantum run (EstimatorV2, Batch) |
| `05_figures.py` | all poster figures |
