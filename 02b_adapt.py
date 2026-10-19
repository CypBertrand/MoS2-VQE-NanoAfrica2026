"""Stage 2b: ADAPT-style compact ansatz — greedily add the excitation with largest energy gradient,
re-optimise, record error vs number of gates (accuracy-vs-depth trade-off for hardware)."""
import numpy as np, json, sys
import pennylane as qml
from pennylane import numpy as pnp
from scipy.optimize import minimize
name,ne,no=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]); kmax=int(sys.argv[4]) if len(sys.argv)>4 else 6
terms=json.load(open(f"{name}_cas{ne}{no}_H.json")); vqe=json.load(open(f"{name}_cas{ne}{no}_vqe.json"))
nq=2*no; e_ref=vqe["E_CASCI"]
H=qml.dot([c for c,_ in terms],[qml.pauli.string_to_pauli_word(s,wire_map={i:i for i in range(nq)}) for _,s in terms])
hf=qml.qchem.hf_state(ne,nq); singles,doubles=qml.qchem.excitations(ne,nq)
pool=[("D",d) for d in doubles]+[("S",s) for s in singles]
dev=qml.device("default.qubit",wires=nq)
def make(sel):
    @qml.qnode(dev,diff_method="adjoint")
    def c(p):
        qml.BasisState(np.array(hf),wires=range(nq))
        for i,(t,w) in enumerate(sel):
            (qml.DoubleExcitation if t=="D" else qml.SingleExcitation)(p[i],wires=w)
        return qml.expval(H)
    return c
sel=[]; params=np.zeros(0); hist=[{"k":0,"E":float(make([])(pnp.array([]))),}]
hist[0]["err_mHa"]=1000*(hist[0]["E"]-e_ref)
for k in range(kmax):
    grads=[]
    for g in pool:
        if g in sel: grads.append(0); continue
        c=make(sel+[g]); p=pnp.array(np.append(params,0.0),requires_grad=True)
        grads.append(abs(float(qml.grad(c)(p)[-1])))
    j=int(np.argmax(grads))
    if grads[j]<1e-5: break
    sel.append(pool[j]); c=make(sel)
    r=minimize(lambda x: float(c(pnp.array(x))),np.append(params,0.0),jac=lambda x: np.array(qml.grad(c)(pnp.array(x,requires_grad=True))),method="L-BFGS-B",options={"ftol":1e-12,"gtol":1e-8})
    params=r.x
    hist.append({"k":k+1,"added":[pool[j][0],list(map(int,pool[j][1]))],"grad":grads[j],"E":float(r.fun),"err_mHa":1000*(float(r.fun)-e_ref),"params":params.tolist()})
    print(hist[-1]["k"],hist[-1]["added"],f'{hist[-1]["err_mHa"]:.3f} mHa')
json.dump({"E_CASCI":e_ref,"hf":list(map(int,hf)),"steps":hist},open(f"{name}_cas{ne}{no}_adapt.json","w"),indent=1)
