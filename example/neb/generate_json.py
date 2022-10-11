from typing import List, Dict, Any

from molecule import Molecule
from qcelemental.models import Molecule as qcelmol, AtomicResult
import qcportal as ptl
import json
import numpy as np

client = ptl.PortalClient('http://localhost:7777', verify=False)

sp = {
    'program': 'psi4',
    'driver': 'gradient',
    'method': 'hf',
    'basis': '6-31g'
}

bohr2ang = 0.529177210903

geo_M = Molecule('NEB_ready.xyz')
spec = {
    "keywords": {
        "images": 21,
        "spring_constant": 1.0,
        "energy_weighted": True
    },
    "qc_specification": {
        "basis": "sto-3g",
        "keywords": {},
        "method": "b3lyp",
        "program": "psi4"
    },
    "program": "geometric"
}



mols = []

for i, M in enumerate(geo_M):
    geometry = M.xyzs.flatten() / bohr2ang
    M_el = qcelmol(**{'symbols': M.elem, 'geometry': geometry})
    mols.append(M_el)
#M = geo_M[0]
#geometry = M.xyzs[0].flatten() / bohr2ang
#M_el = qcelmol(**{'symbols': geo_M[0].elem, 'geometry': geometry})
#

#meta, id = client.add_singlepoints([M_el], **sp)


#recs = client.get_singlepoints(id)
#
#print(recs[0])
meta, id = client.add_nebs([mols], **spec)



#grads_all = np.loadtxt('grads.txt').tolist()
#energies_all = np.loadtxt('energies.txt')  # .tolist()
#
# result = {}
# result['molecule'] = qcelmol(**mols[0])
# result['driver'] = 'gradient'
# # result['properties'] = {'return_energy':energies_all[0],
# #                        'return_gradient':grads_all[0]}
# result['return_result'] = grads_all[0]
# result['success'] = True
# result['provenance'] = {
#     "creator": "Psi4",
#     "routine": "psi4.schema_runner.run_qcschema",
#     "version": "1.4a2.dev1136",
# }
# result['model'] = {'method': 'b3lyp', 'basis': 'sto-3g'}
#
#atomic_result = AtomicResult(**result)
#
#schema["molecule"] = mols
#schema["record_type"] = 'neb'
#schema["result"] = result
#schema["specification"] = spec
#
#json_str = json.dumps(schema, indent=4)
#
#fobj = open('test.json', 'w')
#fobj.write(json_str)
