"""Finish the classical stage without CASSCF(6e,6o) (time-boxed): DFT gaps, triplet, RHF, CAS(4e,4o) diagnostics."""
import numpy as np, json, sys
from pyscf import gto, dft, scf, mcscf, lib
lib.num_threads(2); name=sys.argv[1]
mol=gto.M(atom=f"{name}_opt.xyz",basis="def2-svp",ecp={"Mo":"def2-svp"},verbose=3)
out={"system":name,"note":"CASSCF(6e,6o) not converged within time box; (4e,4o) used"}
ho=mol.nelectron//2-1
for xc in ["pbe","pbe0","hse06"]:
    m=dft.RKS(mol); m.xc=xc; m.max_cycle=200; m=m.density_fit().newton(); out[f"E_{xc}"]=m.kernel(); out[f"gap_{xc}_eV"]=(m.mo_energy[ho+1]-m.mo_energy[ho])*27.2114
mt=gto.M(atom=f"{name}_opt.xyz",basis="def2-svp",ecp={"Mo":"def2-svp"},spin=2,verbose=3)
u=dft.UKS(mt); u.xc="pbe0"; u.max_cycle=200; u=u.density_fit().newton(); out["E_pbe0_triplet"]=u.kernel(); out["ST_gap_pbe0_eV"]=(out["E_pbe0_triplet"]-out["E_pbe0"])*27.2114
hf=scf.RHF(mol); hf.max_cycle=200; hf=hf.density_fit().newton(); hf.kernel(); out["E_RHF"]=hf.e_tot
cas=mcscf.CASCI(hf,4,4); cas.mo_coeff=np.load(f"{name}_cas44_mo.npy"); e=cas.kernel()[0]
dm=cas.fcisolver.make_rdm1(cas.ci,4,4); out["CAS44"]={"E_CASCI":e,"nat_occ":sorted(np.linalg.eigvalsh(dm).tolist(),reverse=True)}
json.dump(out,open(f"{name}_classical.json","w"),indent=1); print(json.dumps(out,indent=1))
