"""Check an IBM Quantum job: status, queue estimate and QPU time used.
Usage:  python check_job.py <JOB_ID> [<JOB_ID> ...]"""
import sys
from qiskit_ibm_runtime import QiskitRuntimeService

service = QiskitRuntimeService()
for jid in sys.argv[1:]:
    job = service.job(jid)
    print(f"\n=== Job {jid} ===")
    print("Backend :", job.backend().name)
    print("Status  :", job.status())
    try:
        est = job.usage_estimation or {}
        if est.get("quantum_seconds") is not None:
            print(f"Estimated QPU time : {est['quantum_seconds']:.1f} s")
    except Exception:
        pass
    try:
        m = job.metrics() or {}
        ts = m.get("timestamps", {})
        for k in ("created", "running", "finished"):
            if ts.get(k): print(f"{k.capitalize():8}: {ts[k]}")
        for k in ("position_in_queue", "estimated_start_time", "estimated_completion_time"):
            if m.get(k) is not None: print(f"{k.replace('_',' ').capitalize()}: {m[k]}")
        u = m.get("usage", {})
        if u.get("quantum_seconds") is not None:
            print(f"QPU time used      : {u['quantum_seconds']} s")
    except Exception as e:
        print("(metrics not available yet:", e, ")")
