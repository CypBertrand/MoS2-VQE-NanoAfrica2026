"""Stage 1: geometry (PBE/def2-SVP), DFT single points, CASSCF/CASCI references, diagnostics."""
import numpy as np, json, sys, time
from pyscf import gto, dft, scf, mcscf, lib
from pyscf.geomopt.berny_solver import optimize
lib.num_threads(2)

def mos2_1():
    a=np.deg2rad(114/2); r=2.13
    return [("Mo",(0,0,0)),("S",(r*np.sin(a),0,r*np.cos(a))),("S",(-r*np.sin(a),0,r*np.cos(a)))]
def mos2_3():
    at=[]; rMo=1.62; rS=1.90; z=1.40
    for k in range(3):
        t=np.deg2rad(120*k); at.append(("Mo",(rMo*np.cos(t),rMo*np.sin(t),0)))
    for k in range(3):
        t=np.deg2rad(120*k+60)
        for s in (+1,-1): at.append(("S",(rS*np.cos(t),rS*np.sin(t),s*z)))
    return at

systems={"MoS2":mos2_1(),"Mo3S6":mos2_3()}
name=sys.argv[1]
out={"system":name}
mol=gto.M(atom=systems[name],basis="def2-svp",ecp={"Mo":"def2-svp"},charge=0,spin=0,verbose=3)
t0=time.time()
import os
if os.path.exists(f"{name}_opt.xyz"):
    mol_eq=gto.M(atom=f"{name}_opt.xyz",basis="def2-svp",ecp={"Mo":"def2-svp"})
else:
    mf=dft.RKS(mol); mf.xc="pbe"; mf.max_cycle=200; mf.level_shift=0.2; mf.damp=0.3
    mf=mf.density_fit().newton()   # second-order SCF: needed for small-gap Mo clusters
    mol_eq=optimize(mf,maxsteps=40)
out["t_opt_s"]=time.time()-t0
mol=gto.M(atom=mol_eq.atom_coords(unit="Ang").tolist() and [(mol_eq.atom_symbol(i),tuple(mol_eq.atom_coord(i,unit="Ang"))) for i in range(mol_eq.natm)],
          basis="def2-svp",ecp={"Mo":"def2-svp"},spin=0,verbose=3)
os.path.exists(f'{name}_opt.xyz') or open(f"{name}_opt.xyz","w").write(f"{mol.natm}\nPBE/def2-SVP opt\n"+"\n".join(f"{mol.atom_symbol(i)} {c[0]:.5f} {c[1]:.5f} {c[2]:.5f}" for i,c in enumerate(mol.atom_coords(unit='Ang'))))
# DFT functional comparison + HOMO-LUMO gaps
for xc in ["pbe","pbe0","hse06"]:
    m=dft.RKS(mol); m.xc=xc; m.max_cycle=200; m=m.density_fit().newton(); e=m.kernel()
    ho=mol.nelectron//2-1
    out[f"E_{xc}"]=e; out[f"gap_{xc}_eV"]=(m.mo_energy[ho+1]-m.mo_energy[ho])*27.2114
# triplet check (PBE0, unrestricted)
mt=gto.M(atom=mol.atom,basis="def2-svp",ecp={"Mo":"def2-svp"},spin=2,verbose=0)
u=dft.UKS(mt); u.xc="pbe0"; u.max_cycle=200; u=u.density_fit().newton(); out["E_pbe0_triplet"]=u.kernel()
out["ST_gap_pbe0_eV"]=(out["E_pbe0_triplet"]-out["E_pbe0"])*27.2114
# RHF + CASSCF(n,n)
hf=scf.RHF(mol); hf.max_cycle=200; hf=hf.density_fit().newton(); hf.kernel()
out["E_RHF"]=hf.e_tot
for ne,no in [(4,4),(6,6)]:
    mc=mcscf.CASSCF(hf,no,ne); mc.max_cycle_macro=100; mc.kernel()
    # natural occupations in active space
    dm1=mc.fcisolver.make_rdm1(mc.ci,no,ne); occ=np.sort(np.linalg.eigvalsh(dm1))[::-1]
    # CASCI on HF orbitals in same space (for VQE reference consistency we use CASSCF orbitals)
    ci=mcscf.CASCI(hf,no,ne); ci.mo_coeff=mc.mo_coeff; eci=ci.kernel()[0]
    # Mo-d character of active orbitals (Mulliken-like)
    C=mc.mo_coeff[:,mc.ncore:mc.ncore+no]; S=mol.intor("int1e_ovlp")
    lbl=mol.ao_labels(); dmask=np.array(["Mo" in l and " 4d" in l for l in lbl])
    pop=np.einsum("mi,mn,ni->mi",C,S,C); dchar=pop[dmask].sum(0).tolist()
    out[f"CAS{ne}{no}"]={"E_CASSCF":mc.e_tot,"E_CASCI":eci,"nat_occ":occ.tolist(),"Mo4d_char":dchar}
    np.save(f"{name}_cas{ne}{no}_mo.npy",mc.mo_coeff)
json.dump(out,open(f"{name}_classical.json","w"),indent=1)
print(json.dumps(out,indent=1))
