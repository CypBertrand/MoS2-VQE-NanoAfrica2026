"""Stage 3: device-calibrated noisy simulation of the hardware experiment.
- circuits: HF reference (X gates only) and compact ADAPT circuit (A1/A2) or GateFabric/UCCSD export
- transpiled for IBM Eagle (FakeBrisbane = ibm_brisbane calibration snapshot), opt level 3
- shot-based measurement of qubit-wise-commuting Pauli groups (incl. readout error)
- readout mitigation: TREX-style per-qubit rescaling (1 - p01 - p10)
- ZNE: global unitary folding, lambda = 1,3,5; linear + Richardson extrapolation
- R independent repeats -> mean +/- std
Usage: python 03_noisy_zne.py MoS2 4 4 A2 [shots] [repeats]"""
import numpy as np, json, sys
import pennylane as qml
from pennylane import numpy as pnp
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime.fake_provider import FakeBrisbane
name,ne,no,ans=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),sys.argv[4]
shots=int(sys.argv[5]) if len(sys.argv)>5 else 8192; R=int(sys.argv[6]) if len(sys.argv)>6 else 5
terms=json.load(open(f"{name}_cas{ne}{no}_H.json")); vqe=json.load(open(f"{name}_cas{ne}{no}_vqe.json"))
nq=2*no; e_ref=vqe["E_CASCI"]; hfstate=qml.qchem.hf_state(ne,nq)
const=sum(c for c,s in terms if set(s)<={"I"}); terms=[(c,s) for c,s in terms if not set(s)<={"I"}]
# qubit-wise commuting grouping (greedy)
groups=[]
for c,s in sorted(terms,key=lambda t:-abs(t[0])):
    for g in groups:
        if all(a=="I" or b=="I" or a==b for a,b in zip(g["basis"],s)):
            g["basis"]=[b if a=="I" else a for a,b in zip(g["basis"],s)]; g["terms"].append((c,s)); break
    else: groups.append({"basis":list(s),"terms":[(c,s)]})
def compact(sel,thetas):
    qc=QuantumCircuit(nq); used=set()
    for (t,w),th in zip(sel,thetas):
        i,j,a,b=w; assert t=="D" and not used&set(w); used|=set(w)
        qc.ry(np.pi+th,i); qc.cx(i,j); qc.x(a); qc.x(b); qc.cx(i,a); qc.cx(i,b)
    for q in range(nq):
        if hfstate[q] and q not in used: qc.x(q)
    return qc
def pl_export(p,key):
    singles,doubles=qml.qchem.excitations(ne,nq)
    def circ(p):
        if key=="UCCSD": qml.AllSinglesDoubles(p,wires=range(nq),hf_state=hfstate,singles=singles,doubles=doubles)
        else:
            shp=qml.GateFabric.shape(n_layers=int(key[-1]),n_wires=nq); qml.GateFabric(pnp.reshape(p,shp),wires=range(nq),init_state=hfstate,include_pi=False)
    q=qml.to_openqasm(qml.tape.make_qscript(circ)(pnp.array(p)),wires=range(nq),measure_all=False,rotations=False,precision=12)
    return QuantumCircuit.from_qasm_str("\n".join(l for l in q.splitlines() if not l.startswith("gphase")))
out={"system":name,"active":f"({ne}e,{no}o)","ansatz":ans,"backend":"FakeBrisbane (ibm_brisbane calibration snapshot)",
     "shots_per_group":shots,"repeats":R,"n_groups":len(groups),"E_CASCI":e_ref}
if ans.startswith("A"):
    ad=json.load(open(f"{name}_cas{ne}{no}_adapt.json")); k=int(ans[1:]); st=ad["steps"][k]
    sel=[tuple(x["added"]) for x in ad["steps"][1:k+1]]
    jobs=[("HF",compact([],[]),ad["steps"][0]["E"]),("VQE",compact(sel,st["params"]),st["E"])]
else:
    key={"UCCSD":"UCCSD","GF1":"GateFabric_L1","GF2":"GateFabric_L2"}[ans]
    jobs=[("HF",compact([],[]),None),("VQE",pl_export(vqe[key]["params"],key),vqe[key]["E"])]
BK=sys.argv[7] if len(sys.argv)>7 else "FakeBrisbane"
import qiskit_ibm_runtime.fake_provider as fp
backend=getattr(fp,BK)(); sim=AerSimulator.from_backend(backend); out["backend"]=BK
props=backend.properties()
for label,qc,e_ideal in jobs:
    tc=transpile(qc,backend,optimization_level=3,seed_transpiler=7)
    phys=tc.layout.final_index_layout()            # logical -> physical qubit
    ro=np.array([1-props.readout_error(p)*1.0 for p in phys])  # approx symmetric: 1-p01-p10 ~ 1-2*avg -> use calibrated below
    # calibrate readout factors (as TREX does) with |0> and |1> preparation circuits
    cal=[]
    for b in (0,1):
        c=QuantumCircuit(backend.num_qubits,nq)
        if b: c.x(phys)
        c.measure(phys,range(nq)); cal.append(c)
    rc=sim.run(transpile(cal,backend,optimization_level=0),shots=20000,seed_simulator=1).result()
    fac=[]
    for q in range(nq):
        p10=sum(v for k,v in rc.get_counts(0).items() if k[::-1][q]=="1")/20000
        p01=sum(v for k,v in rc.get_counts(1).items() if k[::-1][q]=="0")/20000
        fac.append(1-p10-p01)
    base=tc.copy(); base.remove_final_measurements(inplace=True)
    res={lam:[] for lam in (1,3,5)}; res_m={lam:[] for lam in (1,3,5)}
    circs=[]; meta=[]
    for lam in (1,3,5):
        f=base.copy()
        for _ in range((lam-1)//2): f=f.compose(base.inverse()).compose(base)
        for gi,g in enumerate(groups):
            c=QuantumCircuit(backend.num_qubits,nq); c.compose(f,inplace=True); c.barrier()
            for q,bq in enumerate(g["basis"]):
                if bq=="X": c.h(phys[q])
                elif bq=="Y": c.sdg(phys[q]); c.h(phys[q])
            c.measure(phys,range(nq)); circs.append(c); meta.append((lam,gi))
    circs=transpile(circs,backend,optimization_level=0)
    for r in range(R):
        rr=sim.run(circs,shots=shots,seed_simulator=100+r).result()
        acc={lam:const for lam in (1,3,5)}; accm={lam:const for lam in (1,3,5)}
        for ci,(lam,gi) in enumerate(meta):
            cnt=rr.get_counts(ci); tot=sum(cnt.values())
            bits=[(k[::-1],v) for k,v in cnt.items()]
            for c0,s in groups[gi]["terms"]:
                qs=[q for q,p in enumerate(s) if p!="I"]
                ev=sum(v*(-1)**sum(int(k[q]) for q in qs) for k,v in bits)/tot
                acc[lam]+=c0*ev; accm[lam]+=c0*ev/np.prod([fac[q] for q in qs])
        for lam in acc: res[lam].append(acc[lam]); res_m[lam].append(accm[lam])
    lam=np.array([1,3,5.])
    def zne(d):
        Y=np.array([d[l] for l in (1,3,5)]).T  # R x 3
        lin=[np.polyval(np.polyfit(lam,y,1),0) for y in Y]; ric=[np.polyval(np.polyfit(lam,y,2),0) for y in Y]
        return lin,ric
    lin,ric=zne(res_m)
    st=lambda v: [1000*(float(np.mean(v))-e_ref),1000*float(np.std(v,ddof=1))]
    ops=tc.count_ops()
    out[label]={"ECR":int(ops.get("ecr",0)+ops.get("cz",0)),"depth":tc.depth(),"physical_qubits":list(map(int,phys)),"readout_factors":fac,
        "E_ideal":e_ideal,"err_ideal_mHa":None if e_ideal is None else 1000*(e_ideal-e_ref),
        "err_mHa_mean_std":{"raw":st(res[1]),"readout_mit":st(res_m[1]),"lam3_mit":st(res_m[3]),"lam5_mit":st(res_m[5]),
                             "ZNE_linear":st(lin),"ZNE_richardson":st(ric)}}
    print(label,json.dumps({k:v for k,v in out[label].items() if k not in("readout_factors",)},indent=0))
json.dump(out,open(f"{name}_cas{ne}{no}_{ans}_{BK}_noisy.json","w"),indent=1)
