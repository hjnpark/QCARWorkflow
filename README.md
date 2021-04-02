<img src="images/QCARW.png" width = "700">

This package can refine chemical reaction pathways from molecular dynamic (MD) simulation trajectories.[1] It employs QCArchive Infrastructure to achieve efficient data storage and computing resource distribution.[2] Current version can communicate with geomeTRIC[3] and Psi4[4] to refine MD trajectories (xyz files) to provide smoothed pathways. QCARWorkflow consists of two main steps which are `optimization` and `smoothing`. The smoothed pathways then can be used for the Nudged Elastic Band (NEB) method to locate transition state (TS) structures roughly. The NEB method will be implemented in the next version along with the TS optimization calculation step and the Intrinsic Reaction Coordinates (IRC) method. Compatibility with other computational software packages such as Q-Chem and TeraChem will be added as well.    


[1] Wang, L.-P.; McGibbon, R. T.; Pande, V. S.; Martinez, T.J. Automated Discovery and Refinement of Reactive Molecular Dynamics Pathways. *J. Chem. Theory Comput.* **2016**, 12(2), 638–649.[https://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00830](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00830)  
[2] Smith, D. G. A.; Altarawy, D.; Burns, L. A.; Welborn, M.;Naden, L. N.; Ward, L.; Ellis, S.; Pritchard, B. P.; Crawford,T. D. The MOLSSI QCARCHIVE Project: An Open‐source Platform to Compute, Organize, and Share Quantum Chemistry Data. *WIREs Comput. Mol. Sci.* **2020**.[https://doi.org/10.1002/wcms.1491](https://doi.org/10.1002/wcms.1491)   
[3] Wang, L.-P.; Song, C. Geometry Optimization Made Simple with Translation and Rotation Coordinates. *J. Chem. Phys.* **2016**, 144 (21), 214108.[https://doi.org/10.1063/1.4952956](https://doi.org/10.1063/1.4952956)  
[4] Turney, J.M.; Simmonett, A.C.; Parrish, R.M.; Hohenstein, E.G.; Evangelista, F.A.; Fermann, J.T.; Mintz, B.J.; Burns, L.A.; Wilke, J.J.; Abrams, M.L.; Russ, N.J.; Leininger, M.L.; Janssen, C.L.; Seidl, E.T.; Allen, W.D.; Schaefer, H.F.; King, R.A.; Valeev, E.F.; Sherrill, C.D.; Crawford, T.D. Psi4: an open‐source ab initio electronic structure program. *WIREs Comput. Mol. Sci.* **2012** 2: 556-565.[https://doi.org/10.1002/wcms.93](https://doi.org/10.1002/wcms.93)

Author: Heejune Park

Contact Email: heepark@ucdavis.edu

## Quick Set Up
### 1. Creating a Psi4 environment and activating it 


```shell
conda update conda
conda create -n p4env python=3.7.9 psi4 -c psi4/label/dev
conda activate p4env
```
Note that we need the "nightly build" Psi4 channel to make QCARW work.

### 2. Installing QCArchive Infrastructure

In the Psi4 environment, install following compartments.
```shell
conda install qcfractal -c conda-forge
conda install qcportal -c conda-forge
```

### 3. Installing geomeTRIC

geomeTRIC github repository:

[https://github.com/leeping/geomeTRIC](https://github.com/leeping/geomeTRIC)

### 4. Installing QCARWorkflow
Run the command below in ~/QCARWorkflow directory after downloading this repository to a local machine.
```shell
python setup.py install

```
### 5. Other Dependencies
SciPy needs to be installed.
```shell
conda install -c anaconda scipy
```
Done!

## User Guide
### 1. Setting up a QCFractal server 
QCFractal server handles data storage and job submissions.
```shell
qcfractal-server init
nohup qcfractal-server start &
```
The commands above will create a server and run the server. Once the server is started, it will be running until either a user kills it manually or the machine shuts down. Next, we can add users to the server. 
```shell
qcfractal-server user add User1 --password 1234 --permissions admin
``` 
`User1` is the ID and `--password 1234` is the password.
If you don't provide a password, it will automatically generate a long password for the user (I recommend you to assign one).
`--permissions` is user permissions for the `qcfractal-server`. `admin` allows all the possible permissions (read, write, compute, and queue).
More details can be found [here](http://docs.qcarchive.molssi.org/projects/QCFractal/en/stable/server_user.html). 
    
### 2. Running optimization calculations
Once the package is successfully installed and the server is up and running, the python script below will make a dataset named "ds_test" and submit optimization jobs. 
```python
from QCARWorkflow.qcaw import User, Dataset, Workflow

client = User(user = 'User1', password = '1234').server() # This will connect you to the server
ds = Dataset(name = 'ds_test', ds_type = 'OptimizationDataset', client = client).setting('make')

"""
name : Name of the dataset
ds_type : Dataset type. For now "Dataset" and "OptimizationDataset" are supported.
client : client object that we defined above.
setting : 'make' will create a dataset with a given name. 'load' will load the datset. 'delete' will delete the dataset. 'reset' will delete and re-create the dataset.
"""

wf = Workflow(ds = ds, client = client, spec_name = 'spec_test', initial = 'some_example.xyz')

"""
spec_name : Calculation specification name. The specification includes method, basis, program, calculation type, and some other details. Multiple sets of different molecules can be optimized with the same specification. 
initial : This is the xyz file name that contains chemical reactions which will be refined.
"""

wf.optimization(method = 'b3lyp', basis = '6-31g(d)', subsample = 10, compute = True)
print (ds.status(collapse=False))

""" 
subsample : Frame interval of subsampling trajectories. You can set it equal to 1 if you want to optimize all of the frames. 
compute : "True" will automatically submit the jobs after specification and molecules are loaded.  
The very last line will print the current status of the dataset. 
"""
```
Once the jobs are ready, [`qcfractal-manager`](http://docs.qcarchive.molssi.org/projects/QCFractal/en/stable/managers.html) needs to be submitted to computing resources. 
```shell
qcfractal-manager --fractal-uri=https://localhost:7777/  --verify False -u User1 -p 1234
```
You might want to change `localhost` to the machine name where `qcfractal-server start` is running if the computing resources and the local machine aren't same. The portnumber here `7777` is default value from the qcfractal-server. Running `qcfractal-server info` in terminal will provide detailed server information with the portnumber. 

### 3. Smoothing procedure
First, let's check the status of the calculations.
```python
from QCARWorkflow.qcaw import User, Dataset, Workflow
    
client = User(user = 'User1', password = '1234').server() # This will connect you to the server
ds = Dataset(name = 'ds_test', ds_type = 'OptimizationDataset', client = client).setting('load') # Note that the setting is "load" now since we already created the ds_test for optimizations. 
print (ds.status(collapse = False))
```
Running the script above will show job status of the dataset. Once all the jobs are 'COMPLETED', running the script below will detect reactions and smooth them.  
```python
from QCARWorkflow.qcaw import User, Dataset, Workflow
    
client = User(user = 'User1', password = '1234').server() # This will connect you to the server
ds = Dataset(name = 'ds_test', ds_type = 'OptimizationDataset', client = client).setting('load')
wf = Workflow(ds = ds, client = client, spec_name = 'spec_test', initial = 'some_example.xyz').smoothing() # Smoothing function

```
It will create directories /some_example/xx-xx/ and write final result xyz files in there. There are example jupyter notbook and python scripts in /example directory.

## Useful Tips
### 1. Resubmitting "ERROR" status jobs

The current version won't be able to smooth OptimizationDataset results that contain any "ERROR" status in it. You can resubmit jobs to the server with "ERROR" results. The resubmitted calculations will be carried with the same specification.
```python
from QCARWorkflow.qcaw import User, resubmit_all
    
client = User(user = 'User1', password = '1234').server() # This will connect you to the server
client.resubmit_all()

"""
This will resubmit ALL the failed jobs in all of the dataset in the server.
"""
```
If "ERROR" status persists, you can "reset" the datset with a new specification (or add a new specification to the previous dataset with different methods and bases) and run the calculations again until all the jobs are "COMPLETED" successfully. 

If you want to resubmit failed jobs only from a specific dataset, you can use the script below.
 ```python
from QCARWorkflow.qcaw import User, Dataset, Workflow
    
client = User(user = 'User1', password = '1234').server() # This will connect you to the server
ds = Dataset(name = 'ds_test', ds_type = 'OptimizationDataset', client = client).setting('load')

wf = Workflow(ds = ds, client = client, spec_name = 'spec_test', initial = 'some_example.xyz').resubmit() # resubmit function

"""
This will resubmit the failed jobs in the given dataset and specification to the server.
"""
```
### 2. A set of single point energy calculations can be done.

```python
from QCARWorkflow.qcaw import User, Dataset, Workflow
    
client = User(user = 'User1', password = '1234').server() # This will connect you to the server
ds = Dataset(name = 'sp_ds_test', ds_type = 'Dataset', client = client).setting('make')

"""
Note that the ds_type is "Dataset" now.
""" 

wf = Workflow(ds = ds, client = client, spec_name = 'sp_spec_test', initial = 'some_example.xyz').energy(method = 'b3lyp', basis = '6-31g(d)', compute = True)
print(ds.get_records(method = 'b3lyp', program = 'psi4'))
```
The input xyz file here has to have molecule names in the comment lines. The names are used as the molecule names stored in the dataset for the single point energy calculations. More details about the collections (dataset) can be found [here](http://docs.qcarchive.molssi.org/projects/QCPortal/en/stable/collections.html). The next version of QCARWorkflow will provide more convenient ways of navigating calculation result and status throughout the different types of dataset.

### 3. QCArchive Infrastructure
The full documentation of QCArchive Infrastructure can be found [here](http://docs.qcarchive.molssi.org/en/latest/). It conatins all the information regarding different compartments and APIs. 
  






