"""Stage 2: active-space qubit Hamiltonian (Jordan-Wigner) + noiseless VQE (UCCSD and GateFabric)."""
import numpy as np, json, sys, time
from pyscf import gto, scf, mcscf, lib
import pennylane as qml
from pennylane import numpy as pnp
from scipy.optimize import minimize
lib.num_threads(2)
name, ne, no = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
mol=gto.M(atom=f"{name}_opt.xyz",basis="def2-svp",ecp={"Mo":"def2-svp"},verbose=0)
hf=scf.RHF(mol); hf.max_cycle=200; hf=hf.density_fit().newton(); hf.kernel()
cas=mcscf.CASCI(hf,no,ne); cas.mo_coeff=np.load(f"{name}_cas{ne}{no}_mo.npy")
e_casci=cas.kernel()[0]
# rotate active orbitals to CASCI natural orbitals (folds orbital relaxation into the basis -> shallower circuits)
dm=cas.fcisolver.make_rdm1(cas.ci,no,ne); occ,U=np.linalg.eigh(dm); idx=np.argsort(occ)[::-1]; U=U[:,idx]
C=cas.mo_coeff.copy(); C[:,cas.ncore:cas.ncore+no]=C[:,cas.ncore:cas.ncore+no]@U; cas.mo_coeff=C
e_casci=cas.kernel()[0]
h1,ecore=cas.get_h1eff(); h2=lib.unpack_eri if False else None
from pyscf import ao2mo
h2=ao2mo.restore(1,cas.get_h2eff(),no)          # chemist (pq|rs)
two=np.swapaxes(h2,1,3).transpose(0,2,1,3)     # -> PennyLane convention <pq|rs> (verified below)
two=np.einsum("pqrs->psqr",h2)
H=None
for cand in [np.einsum("pqrs->psqr",h2), np.einsum("pqrs->prsq",h2), np.einsum("pqrs->pqrs",h2), np.einsum("pqrs->prqs",h2)]:
    f=qml.qchem.fermionic_observable(np.array([ecore]),h1,cand)
    Hq=qml.jordan_wigner(f); Hq=qml.simplify(Hq)
    nq=2*no
    ev=np.linalg.eigvalsh(qml.matrix(Hq,wire_order=range(nq)))
    # restrict to correct particle number by checking ground state of N-conserving sector
    Nop=qml.matrix(qml.qchem.particle_number(nq),wire_order=range(nq))
    w,v=np.linalg.eigh(qml.matrix(Hq,wire_order=range(nq)))
    nums=np.real(np.einsum("ik,ij,jk->k",v.conj(),Nop,v))
    e_sector=w[np.isclose(nums,ne)].min()
    if abs(e_sector-e_casci)<1e-6: H=Hq; break
assert H is not None, "integral convention check failed"
res={"system":name,"active":f"({ne}e,{no}o)","qubits":nq,"E_CASCI":e_casci,"E_exactdiag":float(e_sector),
     "E_RHF":hf.e_tot,"n_pauli":len(H.operands) if hasattr(H,'operands') else None}
groups=qml.pauli.group_observables([o for o in H.terms()[1]], grouping_type="qwc")
res["n_qwc_groups"]=len(groups)
hfstate=qml.qchem.hf_state(ne,nq)
res["E_HF_qubit"]=float(np.real(np.diag(qml.matrix(H,wire_order=range(nq)))[int("".join(map(str,hfstate)),2)]))
dev=qml.device("default.qubit",wires=nq)
singles,doubles=qml.qchem.excitations(ne,nq)
def run(circ,nparam,label):
    @qml.qnode(dev,diff_method="adjoint")
    def cost(p): circ(p); return qml.expval(H)
    g=qml.grad(cost); hist=[]
    t=time.time()
    r=minimize(lambda x: float(cost(pnp.array(x))), np.zeros(nparam), jac=lambda x: np.array(g(pnp.array(x,requires_grad=True))),
               method="L-BFGS-B", callback=lambda x: hist.append(float(cost(pnp.array(x)))), options={"ftol":1e-12,"gtol":1e-7,"maxiter":500})
    specs=qml.specs(qml.QNode(lambda p: (circ(p),qml.expval(H))[1],dev))(pnp.array(r.x))["resources"]
    tape=qml.workflow.construct_tape(qml.QNode(lambda p:(circ(p),qml.expval(H))[1],dev))(pnp.array(r.x))
    return {"E":float(r.fun),"err_mHa":1000*(float(r.fun)-e_casci),"nparam":nparam,"iters":len(hist),"hist":[e-e_casci for e in hist],
            "time_s":time.time()-t,"params":r.x.tolist()}
def uccsd(p): qml.AllSinglesDoubles(p,wires=range(nq),hf_state=hfstate,singles=singles,doubles=doubles)
res["UCCSD"]=run(uccsd,len(singles)+len(doubles),"UCCSD")
for L in (1,2):
    shp=qml.GateFabric.shape(n_layers=L,n_wires=nq)
    def gf(p,L=L,shp=shp): qml.GateFabric(pnp.reshape(p,shp),wires=range(nq),init_state=hfstate,include_pi=False)
    res[f"GateFabric_L{L}"]=run(gf,int(np.prod(shp)),"GF")
qml.qchem.molecular_hamiltonian  # noqa
json.dump(res,open(f"{name}_cas{ne}{no}_vqe.json","w"),indent=1)
# save hamiltonian as pauli list for qiskit stage
terms=[(float(np.real(c)),qml.pauli.pauli_word_to_string(o,wire_map={i:i for i in range(nq)})) for c,o in zip(*H.terms())]
json.dump(terms,open(f"{name}_cas{ne}{no}_H.json","w"))
print(json.dumps({k:(v if not isinstance(v,dict) else {kk:vv for kk,vv in v.items() if kk not in('hist','params')}) for k,v in res.items()},indent=1))
