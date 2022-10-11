from qcfractal.config import read_configuration
from qcfractal.db_socket.socket  import SQLAlchemySocket
from qcportal.records.neb import NEBSpecification, NEBKeywords
from qcportal.records.singlepoint import QCSpecification
#from qcportal.records.optimization import OptimizationSpecification
from qcelemental.models import Molecule
from qcfractal import FractalSnowflake
print(FractalSnowflake)
s = FractalSnowflake()
client = s.client()

mol1 = Molecule(symbols=["H", "H"], geometry=[0, 0, 0, 0, 0, 2])
mol2 = Molecule(symbols=["O", "O"], geometry=[0, 0, 0, 0, 0, 2])
mol3 = Molecule(symbols=["O", "H", "H"], geometry=[0, 0, 0, 0, 0, 2, 0, 2, 0])

chian = [mol1,  mol2, mol3]
sp_spec = QCSpecification(
            program="psi4",
            driver="energy",
            method="hf",
            basis="6-31g",
        )

#sp_spec = OptimizationSpecification(
#    program="geometric",
#    keywords={},
#    qc_specification=sp_spec
#)
    
neb_kw = NEBKeywords()

meta, rxn_id1 = client.add_nebs([chian], "neb", sp_spec, neb_kw)
ids = rxn_id1
s.await_results()

neb = client.get_nebs(ids)
