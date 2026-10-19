"""One-time setup: save your IBM Quantum account on this computer.
Run:  python save_ibm_account.py
Your API key is stored only in C:\\Users\\<you>\\.qiskit\\qiskit-ibm.json - never share that file."""
from getpass import getpass
from qiskit_ibm_runtime import QiskitRuntimeService

print("Paste your IBM Quantum API key (it will stay hidden), then press Enter:")
token = getpass("API key: ").strip()
crn = input("Paste your instance CRN (starts with crn:v1:...), then press Enter: ").strip()

QiskitRuntimeService.save_account(
    channel="ibm_quantum_platform",
    token=token,
    instance=crn,
    set_as_default=True,
    overwrite=True,
)
print("Saved. Testing the connection...")
service = QiskitRuntimeService()
names = [b.name for b in service.backends()]
print("Success! Quantum processors you can use:", names)
