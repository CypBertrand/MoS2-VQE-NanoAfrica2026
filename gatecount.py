import json,sys,numpy as np, pennylane as qml
from pennylane import numpy as pnp
from qiskit import QuantumCircuit, transpile
import qiskit_ibm_runtime.fake_provider as fp
name,ne,no=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]); nq=2*no
v=json.load(open(f"{name}_cas{ne}{no}_vqe.json")); hf=qml.qchem.hf_state(ne,nq); s,d=qml.qchem.excitations(ne,nq)
def exp(circ,p):
    q=qml.to_openqasm(qml.tape.make_qscript(circ)(pnp.array(p)),wires=range(nq),measure_all=False,rotations=False,precision=12)
    return QuantumCircuit.from_qasm_str("\n".join(l for l in q.splitlines() if not l.startswith("gphase")))
out={}
for key in ["UCCSD","GateFabric_L1","GateFabric_L2"]:
    if key=="UCCSD": c=lambda p: qml.AllSinglesDoubles(p,wires=range(nq),hf_state=hf,singles=s,doubles=d)
    else:
        L=int(key[-1]); shp=qml.GateFabric.shape(n_layers=L,n_wires=nq); c=lambda p,shp=shp: qml.GateFabric(pnp.reshape(p,shp),wires=range(nq),init_state=hf,include_pi=False)
    qc=exp(c,v[key]["params"])
    for bk in ["FakeMarrakesh"]:
        t=transpile(qc,getattr(fp,bk)(),optimization_level=3,seed_transpiler=7); o=t.count_ops()
        out[key]={"2q":int(o.get("cz",0)+o.get("ecr",0)),"depth":t.depth(),"err_mHa":v[key]["err_mHa"],"nparam":v[key]["nparam"]}
print(json.dumps(out)); json.dump(out,open(f"{name}_cas{ne}{no}_gatecounts.json","w"))
