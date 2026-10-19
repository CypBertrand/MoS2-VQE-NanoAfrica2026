"""Stage 4: run the compact VQE circuit on real IBM Quantum hardware (Qiskit Runtime EstimatorV2).
Parameters come from the noiseless ADAPT optimisation (no on-device optimisation -> ~1-3 QPU minutes total).
Usage:
    python 04_hardware.py MoS2 4 4 A2 [backend_name]
First time only:  QiskitRuntimeService.save_account(channel="ibm_quantum_platform", token="<YOUR_TOKEN>", instance="<CRN>")
"""
import numpy as np, json, sys, datetime

from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
import os
from qiskit_ibm_runtime import QiskitRuntimeService, EstimatorV2
name,ne,no,ans=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),sys.argv[4]
nq=2*no; hfstate=[1]*ne+[0]*(nq-ne)   # Hartree-Fock occupation (Jordan-Wigner, spin-orbital order)
terms=json.load(open(f"{name}_cas{ne}{no}_H.json")); e_ref=json.load(open(f"{name}_cas{ne}{no}_vqe.json"))["E_CASCI"]
const=sum(c for c,s in terms if set(s)<={"I"})
op=SparsePauliOp.from_list([(s[::-1],c) for c,s in terms if not set(s)<={"I"}]).simplify()   # PennyLane wire i -> Qiskit qubit i
ad=json.load(open(f"{name}_cas{ne}{no}_adapt.json")); k=int(ans[1:]); st=ad["steps"][k]
sel=[tuple(x["added"]) for x in ad["steps"][1:k+1]]
def compact(sel,thetas):
    qc=QuantumCircuit(nq); used=set()
    for (t,w),th in zip(sel,thetas):
        i,j,a,b=w; used|=set(w); qc.ry(np.pi+th,i); qc.cx(i,j); qc.x(a); qc.x(b); qc.cx(i,a); qc.cx(i,b)
    for q in range(nq):
        if hfstate[q] and q not in used: qc.x(q)
    return qc
service=QiskitRuntimeService()
backend=service.backend(sys.argv[5]) if len(sys.argv)>5 else service.least_busy(operational=True,simulator=False,min_num_qubits=nq)
print("Backend:",backend.name)
pm=generate_preset_pass_manager(optimization_level=3,backend=backend,seed_transpiler=7)
circs={"HF":compact([],[]),"VQE":compact(sel,st["params"])}
pubs=[]
for lab,qc in circs.items():
    isa=pm.run(qc); pubs.append((isa,op.apply_layout(isa.layout))); print(lab,"ECR/CZ:",isa.count_ops())
SHOTS=int(os.environ.get("SHOTS","2048"))          # per measurement basis; raise only on paid plans
WITH_HF=os.environ.get("WITH_HF","0")=="1"          # HF calibration circuit costs extra QPU time
if not WITH_HF:
    pubs=pubs[1:]; circs={"VQE":circs["VQE"]}
def estimator(mode,mitigated):
    est=EstimatorV2(mode=mode); o=est.options
    o.default_shots=SHOTS
    o.max_execution_time=300                        # hard cap: 5 QPU-minutes per job
    if mitigated:
        o.resilience_level=1                        # TREX readout mitigation
        o.resilience.measure_mitigation=True
        o.resilience.zne_mitigation=True
        o.resilience.zne.noise_factors=(1,3,5)
        o.resilience.zne.extrapolator=("linear","polynomial_degree_2")
        o.twirling.enable_gates=True; o.twirling.enable_measure=True
        o.twirling.num_randomizations=8
        o.dynamical_decoupling.enable=True; o.dynamical_decoupling.sequence_type="XY4"
    else:
        o.resilience_level=0
    return est
tag=f"{name}_cas{ne}{no}_{ans}"
j_raw=estimator(backend,False).run(pubs)
j_mit=estimator(backend,True).run(pubs)
ids={"raw":j_raw.job_id(),"mitigated":j_mit.job_id(),"backend":backend.name,"circuits":list(circs)}
json.dump(ids,open(f"{tag}_jobids.json","w"),indent=1)
print("Job IDs saved to",f"{tag}_jobids.json",":",ids)
print("Waiting for results (you can close this window; retrieve later with: python get_results.py",tag+")")
import get_results
get_results.collect(tag)
