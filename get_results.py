"""Retrieve hardware results for jobs submitted by 04_hardware.py (works even after closing the window).
Usage: python get_results.py MoS2_cas44_A2"""
import json, sys, datetime
from qiskit_ibm_runtime import QiskitRuntimeService
def collect(tag):
    ids=json.load(open(f"{tag}_jobids.json"))
    name,cas=tag.split("_")[0],tag.split("_")[1]
    terms=json.load(open(f"{name}_{cas}_H.json")); e_ref=json.load(open(f"{name}_{cas}_vqe.json"))["E_CASCI"]
    const=sum(c for c,s in terms if set(s)<={"I"})
    service=QiskitRuntimeService()
    out={"backend":ids["backend"],"retrieved_utc":datetime.datetime.utcnow().isoformat(),"jobs":ids,"E_CASCI":e_ref}
    for kind in ("raw","mitigated"):
        job=service.job(ids[kind])
        try:
            res=job.result()
        except Exception as e:
            print(kind,"job",ids[kind],"status",job.status(),"->",str(e)[:200]); out[kind]={"status":str(job.status())}; continue
        try: q=job.metrics().get("usage",{}).get("quantum_seconds")
        except Exception: q=None
        out[kind]={"quantum_seconds":q}
        for i,lab in enumerate(ids["circuits"]):
            E=float(res[i].data.evs)+const; sd=float(res[i].data.stds)
            out[kind][lab]={"E":E,"std_mHa":1000*sd,"err_mHa":1000*(E-e_ref)}
            print(f"{kind:9} {lab}: error vs CASCI = {1000*(E-e_ref):8.2f} +/- {1000*sd:.2f} mHa   (QPU s: {q})")
    fn=f"{tag}_hardware_{ids['backend']}.json"; json.dump(out,open(fn,"w"),indent=1); print("Saved",fn)
if __name__=="__main__": collect(sys.argv[1])
