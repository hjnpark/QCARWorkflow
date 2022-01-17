#!/usr/bin/env python

"""
qcarw.py

"""
import os, sys, subprocess, time, shutil, socket
import qcengine as qcng
import qcelemental as qcel
import qcportal as ptl
import matplotlib.pyplot as plt
import numpy as np
from .params import parse_refine_args
from .molecule import PeriodicTable, Molecule, EqualSpacing, TopEqual, MolEqual
from collections import OrderedDict
from .errors import NoDatasetError, InvalidCommandError, QCFractalError, QCEngineError, OptimizeInputError
from .constant import bohr2ang, au2kcal

def formula(M):
    """
    Provides a molecular formula
        
    Parameters
    ----------
    M : Molecular object

    Return
    ------
    formula : string
    """
    elem_list = M.elem
    name = OrderedDict()
    for elem in elem_list:
        if not elem in name.keys():
            name[elem] = 0
        name[elem] += 1

    formula = ''
    for key, val in name.items():
        formula += key + str(val)
    return formula

def qc_to_geo(qc_M, comment='', b2a=False):
    """
    Convert QCArchive molecule object to geomeTRIC molecule objects
   
    Parameters
    ----------
    qc_M : QCArchive molecule object
        
    comment : string
        Comment for xyz file

    b2a : boolean 
        b2a = True will convery length unit from Bohr to Angstrom

    Return
    ----------
    geo_M : geomeTRIC molecule object
    """
    geo_M = Molecule()
    geo_M.comms = [comment]
    geo_M.elem = list(qc_M.symbols)
    geom = np.array(qc_M.geometry, dtype = np.float64).reshape(-1, 3)
    if b2a == True:
        geom *= bohr2ang 
    geo_M.xyzs = [geom]
    geo_M.build_bonds()
    geo_M.build_topology()
    return geo_M

def resubmit_all(client):
    """
    This function will submit ALL the failed jobs again to server.
    """
    errors = client.query_tasks(status = 'ERROR')
    for i in range(len(errors)):
        client.modify_tasks('restart', errors[i].base_result)
    print ('All the failed calculations were submitted again.')

class User(object):
    """
    This class helps users to connect to the server.
    """
    def __init__(self, user=None, password=None):
        """
        Parameters
        ----------
        user : user ID for QCFractal server
    
        password ; password for QCFractal server
        """
        self.user = user
        self.password = password

    def server(self):
        """
        Setting up or starting a QCFractal server
        
        Return
        ----------
        client : client object 
            With the client object, users can create/access dataset
        """
        info = os.popen('qcfractal-server info').readlines()
        host = socket.gethostname()
        for line in info:
           # if 'host' in line:
           #     host = line.strip().split(' ')[-1]
            if 'port' in line:
                port = int(line.strip().split(' ')[-1])
        try:
            client = ptl.FractalClient('%s:%i'%(host, port), verify = False, username = self.user, password = self.password)
        except:
            client = None
            for i in range(3):
                try:
                    subprocess.run('qcfractal-server start &', shell = True, stdout = subprocess.DEVNULL)
                    time.sleep(3.0)
                    client = ptl.FractalClient('%s:%i'%(host, port), verify = False, username = self.user, password = self.password)
                except:
                    if client != None:
                        break 
                    else:
                        continue
            print("Server is running") 

            if client == None:
                raise QCFractalError("Client could not be claimed properly. Try to restart the qcfractal-server.")
 
        print("Client is ready")
       
        return client           

class Dataset(object):
    """
    This class helps users with handling dataset. 
    """
    def __init__(self, name, client, ds_type='OptimizationDataset'):
        """
        Parameters
        ----------
        name : str
            Name of the dataset
        
        clinet : client object from the User class 

        ds_type : str
            Type of the dataset. Currently 'OptimizationDataset' and 'Dataset' are supported.
        """
        self.name = name
        self.ds_type = ds_type
        self.client = client

    def setting(self, command = None):
        """
        parameters
        ----------
        command : str
            1. 'make' will create a new dataset with a given name.
            2. 'load' will load the dataset with a given name.
            3. 'delete' will delete the dataset with a given name.
            4. 'reset' will delete and re-create the dataset with a give name.

        """
        if command == 'make':
            try:
                if self.ds_type == 'OptimizationDataset':
                    new_ds = ptl.collections.OptimizationDataset(name = self.name, client = self.client)            
                elif self.ds_type == 'Dataset':
                    new_ds = ptl.collections.Dataset(name = self.name, client = self.client)    
                new_ds.save()        
                ds = self.client.get_collection(self.ds_type, name = self.name)
            except:
                ds = self.client.get_collection(self.ds_type, name = self.name)
        elif command == 'load':
            try:
                ds = self.client.get_collection(self.ds_type, name = self.name)
            except:
                raise NoDatasetError("\'%s\' dataset could not be loaded. Try to restart the qcfractal server." % (self.name))
        elif command == 'delete':
            try:
                self.client.delete_collection(self.ds_type, name = self.name)
                ds = None
            except:
                raise NoDatasetError("\'%s\' dataset can't be deleted since there is no %s named %s." % (self.name, self.ds_type, self.name))
        elif command == 'reset':
            try:
                self.client.delete_collection(self.ds_type, name = self.name)
                if self.ds_type == 'OptimizationDataset':
                    ds = ptl.collections.OptimizationDataset(name = self.name, client = self.client)
                elif self.ds_type == 'Dataset':
                    ds = ptl.collections.Dataset(name = self.name, client = self.client) 
            except:
                ds = ptl.collections.OptimizationDataset(name = self.name, client = self.client)
        else:
           raise InvalidCommandError("Please provide a valid command for the dataset.")  
        return ds
    
       

class Workflow(object):
    """
    This class perfroms the automated reaction refinement workflow using QCArchive Infrastructure
    """
    def __init__(self, client=False, ds=None, spec_name=None):
        """
        Parameters
        ----------
        ds : Dataset object
            Dataset object created from the Dataset class

        spec_name : string
            QCSpecification name            

        client : client object
        """
        
        self.client = client
        self.ds = ds
        self.spec_name = spec_name

    def resubmit(self): 
        """
        This function detects ERROR calculation results in a given dataset with a specficiation and submit them again. 
        """
        opts = self.ds.df[self.spec_name].tolist()
        num = 0
        for i in range (len(opts)):
            if opts[i].error != None:
                num += 1
                self.client.modify_tasks('restart', opts[i].id)
        print ("%i failed jobs in \'%s\' dataset with \'%s\' specfication have been submitted again." %(num, self.ds.name, self.spec_name))

    def equal(self, m1, m2):
        return TopEqual(m1, m2) if True else MolEqual(m1, m2)

   # def energy(self, method=None, basis=None, compute=False):
   #     """
   #     This function converts xyz file to QCAI molecule objects and creates single point energy calculation jobs.          

   #     Parameters
   #     -----------
   #     method, basis : string
   #         Electron structure method and basis sets

   #     compute : boolean
   #         'compute = False' will only save the molecules in the given dataset and specfication. 'compute = True' will submit the jobs to server.

   #     spec_overwrite : boolean
   #         'spec_overwrite = True' will overwrite the spec if the same name spec exists

   #     Return
   #     ----------
   #     ds_sp : dataset object
   #         print(ds_opt.status(collapse = False)) will show the dataset's status 

   #     """  
   #     ds_sp = self.ds 
   #     M = self.M 
   #     key = ptl.models.KeywordSet(values = {'maxiter':1000,
   #                                          'e_convergence': 1e-6,
   #                                          'guess' : 'sad',
   #                                          'scf_type' : 'df'})        
   #     ds_sp.add_keywords(self.spec_name, 'psi4', key)

   #     spec = {
   #             'program' : 'psi4',
   #             'method' : method,
   #             'basis' : basis,
   #             'keywords': self.spec_name,
   #             'tag' : None}
   #     frames = list(range(len(M)))

   #     for frm in frames:               
   #         mol = qcel.models.Molecule(**{'symbols': M[frm].elem, 'geometry': np.array(M[frm].xyzs)/0.529177210, 'molecular_charge' : self.charge, 'molecular_multiplicity' : self.mult}) 
   #         #if len(M[frm].comms[0]) == 0 or len(M[frm].comms[0]) > 10:
   #         #    raise RuntimeError('Please provide a short name, less than 10 letters, of the molecule in comment line in the input xyz file (between number of atoms and coordinates).')
   #         try:
   #             ds_sp.add_entry('%s_%s' %(formula(M[frm])+str(frm)), mol)
   #         except:
   #             pass 
   #     ds_sp.save()
   #     if compute:
   #         ds_sp.compute(**spec)
   #         print ("Single point energy calculations in %s have been submitted. Run the QCFractal manager to carry the calculations." %(ds_sp.name))
   #     return ds_sp 
        

    def dsoptimize(self, initial, charge=0, mult=1, method='b3lyp', basis='6-31g(d)', subsample=10, maxiter=500, coordsys='tric', compute=False, spec_overwrite=True): 
        """
        This function converts the xyz file to QCAI molecule objects and sets up optimization jobs.          

        Parameters
        -----------
        initial : string
            MD trajectory (xyz file name) that needs to be refined
    
        charge, mult : int
            Molecular charge and multiplicity
            
        method, basis : string
            Electron structure method and basis sets

        subsample : int
            Frame interval for subsampling trajectories

        compute : boolean
            'compute = False' will only save the molecules in the given dataset and specfication. 'compute = True' will submit the jobs to server.

        maxiter : int
            maximum scf iteration number

        spec_overwrite : boolean
            'spec_overwrite = True' will overwrite the spec if the same name spec exists

        Return
        ----------
        ds_opt : dataset object
            print(ds_opt.status(collapse = False)) will show the dataset's status 

        mass : float
            Molecular mass
        """  

        M = Molecule(initial) 
        ds_opt = self.ds
        key = [ptl.models.KeywordSet(values = {'maxiter':maxiter})]        
        key_id = self.client.add_keywords(key)[0]
        optimize = {
            'name' : self.spec_name,
            'optimization_spec' : {'program': 'geometric', 'keywords': {'coordsys': coordsys}},
            'qc_spec' : {
                    'driver': 'gradient', 
                    'method': method, 
                    'basis': basis,
                    'keywords': key_id,
                    'program': 'psi4'
                     }
                       }
        try:
            ds_opt.add_specification(**optimize, overwrite = spec_overwrite)
            print("Specification %s was added into %s" %(self.spec_name, ds_opt))
        except:
            print("Specification %s is either already added or it has key values that are not allowed" %self.spec_name)
            pass

        mass = sum([PeriodicTable.get(M.elem[i], 0.0) for i in range(M.na)])
        frames = list(range(0, len(M), subsample))
        if (len(M)-1) not in frames:
            frames.append(len(M)-1)
        for frm in frames:               
            mol = qcel.models.Molecule(**{'symbols': M[frm].elem, 'geometry': np.array(M[frm].xyzs)/bohr2ang,  'molecular_charge' : charge, 'molecular_multiplicity' : mult}) 
            try:
                ds_opt.add_entry('%s_%i' %(initial.split('.')[0], frm), mol, save = False)
            except:
                pass 
        ds_opt.save()
        if compute:
            ds_opt.compute(self.spec_name)
            print ("Calculations in \'%s\' with \'%s\' specification have been submitted. Run the QCFractal manager to carry the calculations." %(ds_opt.name, self.spec_name))
        return ds_opt, mass

    
    def smoothing(self, initial):
        """
        Once the optimization is done, smoothing function will detect reactions and smooth them for the NEB calculation. 

        Parameters
        -----------
        initial : string
            MD trajectory (xyz file name) that was refined with the dataset
 
        """
        M = Molecule(initial)
        opt = self.ds.status(self.spec_name, collapse = False)
        opts = self.ds.df[self.spec_name].tolist()
        stats = [opt.status for opt in opts] 
        mol_names = opt.index.tolist() 
        mol_name = '_'.join(str(elem) for elem in mol_names[0].split('_')[:-1]) 

        OptMols = OrderedDict()
        Iter = OrderedDict()
        err = 0
        for name, calc in zip(mol_names, stats):
            frm = int(name.split('_')[-1])
            if calc.value == 'ERROR' or calc.value == 'INCOMPLETE':
                Iter[frm] = 'ERROR'
                err += 1
                continue
            record = self.ds.get_record(name, self.spec_name)
            init_M = record.get_initial_molecule()
            geo_M = qc_to_geo(init_M, b2a = True) 
        
            input_check = np.allclose(geo_M.xyzs[0], M[frm].xyzs[0])
            if not input_check:
                print('Please double check your MD trajectory (xyz file) and dataset name.')
            qcel_M = record.get_final_molecule()
            OptMols[frm] = qc_to_geo(qcel_M, b2a = True)  
            Iter[frm] = len(self.ds.get_record(name, self.spec_name).trajectory)

        if err > 0 : 
            print ("WARNING: %i ERROR or/and INCOMPLETE results detected." % err) 

        if os.path.exists('%s_iterations.txt' %mol_name):
            os.remove('%s_iterations.txt' %mol_name)
            
        with open ('%s_iterations.txt' %mol_name, 'w') as fn:
            fn.write('Optimization Iteration Numbers of %s\n' %mol_name)
            for frame, ite in Iter.items():
                fn.write('%i : %s\n' %(frame, ite))

        print ("Detecting reactions now.")

        path_initial = [] 
        path_final = []
        for fi, fj in zip(list(OptMols.keys())[:-1], list(OptMols.keys())[1:]): 
            if not self.equal(OptMols[fi], OptMols[fj]): 
                path_initial.append(fi)
                path_final.append(fj)

        MolPairs = []
        FramePairs = []
            
        for fi in path_initial:
            for fj in path_final:
                if fj > fi and (not self.equal(OptMols[fi], OptMols[fj])):
                    if (fj - fi) > 1000: continue
                    NewPair = True
                    for imp, (m1, m2) in enumerate(MolPairs):
                        if self.equal(OptMols[fi], m1) and self.equal(OptMols[fj], m2):
                            FramePairs[imp].append((fi, fj))
                            NewPair = False
                            break
                        elif self.equal(OptMols[fi], m2) and self.equal(OptMols[fj], m1):
                            FramePairs[imp].append((fi, fj))
                            NewPair = False
                            break
                    if NewPair:
                        MolPairs.append((OptMols[fi], OptMols[fj]))
                        FramePairs.append([(fi, fj)])

        if len(MolPairs) != len(FramePairs) or len (MolPairs) == 0:
            raise RuntimeError ("No reactions are detected or the Number of detected pairs of reacting molecules and frames don't match.")            

        geo_mol_Traj = None

        path = './%s/' %(mol_name)
        if os.path.exists(path):
            shutil.rmtree(path) 
        os.mkdir(path)
        neb_inputs = {}

        for i in range(len(MolPairs)): 
            (a,b) = FramePairs[i][np.argmin([(jb-ja) for (ja, jb) in FramePairs[i]])]
            qc_mol_Traj1 = self.ds.get_record(mol_name + '_' + str(a), self.spec_name).get_molecular_trajectory()
            qc_mol_Traj2 = self.ds.get_record(mol_name + '_' + str(b), self.spec_name).get_molecular_trajectory()

            for j in range(len(qc_mol_Traj1)-1):
                if geo_mol_Traj == None:
                    geo_mol_Traj = qc_to_geo(qc_mol_Traj1[-1], b2a = True)
                geo_mol_Traj += qc_to_geo(qc_mol_Traj1[::-1][j+1], b2a = True) 
            geo_mol_Traj += M[a:b]
            for k in range(len(qc_mol_Traj2)):
                geo_mol_Traj += qc_to_geo(qc_mol_Traj2[k], b2a = True)    
            
            fnum =  str(a) + '-' + str(b)
            #fname = str(mol_name +'_'+ fnum)
            path = './%s/' %(mol_name)
            
            NEB_path = path + fnum
            os.mkdir(NEB_path)
            geo_mol_Traj.write(os.path.join(path + fnum,'connected.xyz'))
            equal = EqualSpacing(geo_mol_Traj, dx = 0.05) 
            equal.write(os.path.join(path + fnum, 'spaced.xyz'))
            geo_mol_Traj = None 
            command ='Nebterpolate.py --morse 1e-2 --repulsive --allpairs --anchor 2 %s/spaced.xyz %s/NEB_ready.xyz &> %s/interpolate.log' %(NEB_path, NEB_path, NEB_path)
            log = open('%s/interpolate.log' %NEB_path, 'a')
            #err = open('%s/interpolate_%s.log' %(NEB_path, fname), 'a')
            subprocess.Popen(command, shell = True, stdout = log, stderr = log)
            neb_inputs[i] = NEB_path + '/NEB_ready.xyz'
        print("Smoothing Procedure is running on the local machine. NEB ready xyz files will be generated once the smoothing procedure is done.")
        return neb_inputs
        

    def neb(self, initial, charge=0, mult=1, method='b3lyp', basis='6-31+g(d)', images=21, coordsys='cart', ew = False, nebk = 1, avgg=0.025, maxg=0.05, guessk=0.01, guessw=0.5, tmpdir=None):
        """
        This function will run NEB calculations to locate rough transition state structures.          

        Parameters
        -----------
        initial : str
            Name of the neb ready xyz file.
        
        charge, mult : int
            Molecular charge and multiplicity

        method, basis : string
            Electron structure method and basis sets

        images : integer
            Number of images for the neb chain. 

        coordsys : str
            'cart': Cartesian Coordinates
            'prim': Primitive (a.k.a. redundant) Coordinates
            'tric': Translation-Rotational Coordinates
            'dlc' : Delocalized Internal Coordinates
            'hdlc': Hybrid Delocalized Internal Coordinates

        ew : boolean
            True will perfrom the energy weighted NEB calculation

        nebk : float
            Spring constant in kcal/mol/Ang^2

        avgg, maxg : float
            Average RMS-gradient and max RMS-gradient for convergence. Unit of eV/Ang

        guessk, guessw : float
            Guess Hessian and guess weight for chain coordinates.
    
        tmpdir : str
            Temporary directory for NEB results.

        Return
        ----------
        it will write each iterated chains in .tmp directory and final transition state geometries.

        """  
        band = Molecule(initial) 
        inp = '.'.join(initial.split('.')[:-1])
        
        qcel_mol = qcel.models.Molecule(**{'symbols': band.elem, 'geometry': np.array(band[0].xyzs)/bohr2ang,  'molecular_charge' : charge, 'molecular_multiplicity' : mult}) 

        neb_input = {
            'keywords' : {
                'program' : 'psi4',
                'images': images,
                'avgg' : avgg,
                'maxg' : maxg,
                'nebk' : nebk, 
                'maxcyc' : 200,
                'guessk' : guessk,
                'guessw' : guessw,
                'coords': initial,
                'coordsys': coordsys,
                'engine': 'qcengine',
                'client': self.client, 
                'input': inp
                        },

                
            'input_specification':{
                'driver': 'gradient',
                'model' : {
                    'method': method,
                    'basis': basis
                    },
                                    },
            'initial_molecule':qcel_mol
            }


        if ew:
            neb_input['keywords']['ew'] = 'yes'
        neb_procedure = qcng.compute_procedure(neb_input, 'geometric')
        return neb_procedure

    def optimize(self, initial, charge=0, mult=1, method='b3lyp', basis='6-31+g(d)', coordsys='tric', maxiter=500, ts=False):
        """
        This function will run a single optimization procedure        

        Parameters
        -----------
        initial : str or geomeTRIC Molecule object
            Name of the xyz file containing the initial structure.

        charge, mult : int
            Molecular charge and multiplicity

        method, basis : string
            Electron structure method and basis sets

        maxiter : int
            Maximum iteration number for scf calculations.
        
        ts : boolean
            True will perform a transition state structure optimization

        Return
        ----------
        M_qc : QCAI Molecule object of the optimized structure
        M_geo : geomeTRIC Molecule object of the optimized structure
        energy : energy of the optimized TS structure

        """  
        if ts:
            method = 'ts-' + method
        if isinstance(initial, str):
            init = Molecule(initial) 
        elif isinstance(initial, Molecule):
            init = initial
        else:
            raise OptimizeInputError('Please provide either xyz file name or a geomeTRIC Molecule object.')

        qcel_mol = qcel.models.Molecule(**{'symbols': init.elem, 'geometry': np.array(init.xyzs)/bohr2ang,  'molecular_charge' : charge, 'molecular_multiplicity' : mult}) 
        if not self.client:
            """
            Optimization calculation will be carried locally
            """
            opt_input = {
                'keywords' : {
                    'program' : 'psi4',
                    'coordsys' : coordsys
                    },
                'input_specification':{
                    'driver': 'gradient',
                    'model' : {
                        'method': method,
                        'basis': basis
                        }
                    
                    },
                'initial_molecule':qcel_mol
                }


            opt_result = qcng.compute_procedure(opt_input, 'geometric') #OptimizationResult
            if type(opt_result).__name__ == "FailedOperation":
                raise QCEngineError(opt_result.error)
            else:
                energy = opt_result.energies[-1]
                M_qc = opt_result.final_molecule
                M_geo = qc_to_geo(M_qc, comment='Energy : %.7f Hartree' %energy, b2a=True)
        else:
            """
            Optimization will be carried with QCFractal server
            """
            import time
            opt_qcschema = {
                 'keywords' : {'coordsys' : coordsys},
                 'qc_spec': {
                    'driver': 'gradient',
                    'method': method,
                    'basis': basis,
                    'program': 'psi4',
                    }
                }


            r = self.client.add_procedure('optimization', 'geometric', opt_qcschema, [qcel_mol])
            proc_id = r.ids
            loop = 0 
            while True:
                proc = self.client.query_procedures(id=proc_id)[0] #OptimizationRecord
                status = proc.status.split('.')[-1].upper().strip()
                if status == 'INCOMPLETE':
                    time.sleep(50)
                    loop += 1
                elif status == 'ERROR':
                    print('Error detected')
                    res = self.client.modify_tasks('restart',proc.id)
                    print(res.n_updated,"ERROR status optimization resubmitted")
                    loop += 1
                elif status == 'COMPLETE': 
                    energy = proc.get_final_energy()
                    M_qc = proc.get_final_molecule()
                    M_geo = qc_to_geo(M_qc, comment= 'Energy : %.7f Hartree' %energy, b2a=True)     
                    print("QCAI Optimization is done.")
                    break
                if loop > 1000:
                    raise QCFractalError('Optimization failed.')
 
        return M_qc, M_geo, energy

def main():
    args_dict = parse_refine_args(sys.argv[1:])

    input_name = args_dict['input'].split('.')
    temp_name = '.'.join(input_name[:-1]) 
    spec_name = temp_name + '_spec'

    if 'dataset' in args_dict:
        pass
    else:
        args_dict['dataset'] = temp_name + '_ds'

    dataset     = args_dict.get('dataset')
    user        = args_dict.get('user'      , None)
    password    = args_dict.get('passowrd'  , None)
    initial     = args_dict.get('input')
    charge      = args_dict.get('charge'    , 0)
    mult        = args_dict.get('mult'      , 1)
    subsample   = args_dict.get('subsample' , 10)
    maxiter     = args_dict.get('maxiter'   , 500)
    optmethod   = args_dict.get('optmethod' , 'b3lyp')
    optbasis    = args_dict.get('optbasis'  , '6-31g(d)')
    optcrdsys   = args_dict.get('optcrdsys' , 'tric')
    tsmethod    = args_dict.get('tsmethod'  , 'b3lyp')
    tsbasis     = args_dict.get('tsbasis'   , '6-31+g(d)')
    images      = args_dict.get('images'    , 21)
    nebmethod   = args_dict.get('nebmethod' , 'b3lyp')
    nebbasis    = args_dict.get('nebbasis'  , '6-31+g(d)')
    nebcrdsys   = args_dict.get('nebcrdsys' , 'cart')
    nebk        = args_dict.get('nebk'      , 1)
    avgg        = args_dict.get('avgg'      , 0.025)
    maxg        = args_dict.get('maxg'      , 0.05)
    ew          = args_dict.get('ew'        , False)

    client = User(user, password).server()
    
    ds = Dataset(dataset, client).setting('make')
    wf = Workflow(ds=ds, client=client, spec_name=spec_name)
    ds, Mmass =wf.dsoptimize(initial=initial, charge=charge, mult=mult, method=optmethod, basis=optbasis, subsample=subsample, coordsys=optcrdsys, compute = True)
    
    cycle = 0
    while True: 
        ds = Dataset(dataset, client).setting('load')
        if cycle % 10 == 0 :
            print('Dataset Status')
            print(ds.status(collapse=False))

        ds.status(spec_name)
        opts = ds.df[spec_name].tolist() #OptimizationRecord in a list
        num_opt = len(opts)
        comp = 0
        for opt in opts:
            stat = opt.status.upper().split('.')[-1]
            if stat == 'ERROR':
                client.modify_tasks('restart', opt.id)
            if stat == 'COMPLETE':
                comp += 1
        
        if comp == num_opt or comp/num_opt> 0.8:
            print('Optimization step is completed.')
            break
        print("%i/%i calculations are completed" %(comp,num_opt)) 
        wait = (num_opt-comp)*int(Mmass*0.5)
        print("Molecular mass = %.2f / Waiting %i seconds" %(Mmass, wait))
        time.sleep(wait) 

        if cycle > 1000:
            raise RuntimeError('Stuck in a while loop checking OptimizationDataset status.')

        if comp == 0 and cycle > 10 :
            raise QCFractalError('Jobs are not recognized by QCFractal server. Try to restart the server and manager')
        cycle += 1

    smoothed = wf.smoothing(initial=initial)
    time.sleep(30)
    neb_num = len(smoothed)
    print("Number of reactions detected: %i" %neb_num)
    smoothed_list = list(smoothed.values()) 
    guess_ts_list = []
    tmp_list = [] 
    for neb in smoothed_list:
        inp = '.'.join(neb.split('.')[:-1]) 
        print("Generating %s" %neb) 
        smoothing_cycle = 0
        while True:
            if os.path.exists(neb):
                break 
            elif smoothing_cycle > 1000:
                f = open(inp + '.error','w')
                f.write('Smoothing procedure error. It probably just copied the input file as result. Check the .log file.')
                break
            else:
                smoothing_cycle += 1
                time.sleep(Mmass*0.5)
        print("NEB method will be used to refine %s initial chain" %neb)
        wf.neb(initial=neb, charge=charge, mult=mult, method=nebmethod, basis=nebbasis, images=images, coordsys=nebcrdsys, ew=ew, nebk=nebk, avgg=avgg, maxg=maxg) 
        print("%s initial chain is processed." %neb) 
        ts_name = inp + '.tsClimb.xyz'
        tmp_name = inp + '.tmp'
        guess_ts_list.append(ts_name)
        tmp_list.append(tmp_name)
    
    for i, ts in enumerate(guess_ts_list):
        if os.path.exists(ts): #Sometimes the smoothing function won't be able to smooth the pass. If there are paths that were not smoothed, it will just skip them.
            pass
        else:
            continue 
        inp = '~/' + '/'.join(('.'.join(tmp_list[i].split('.')[:-1]) + '.xyz').split('/')[-3:-1])
        print(inp)
        tmp_dir = '/'.join(tmp_list[i].split('/')[:-1])
        M_qc_ts, M_geo_ts, E_ts= wf.optimize(ts, charge=charge, mult=mult, method=tsmethod, basis=tsbasis, ts=True)
        M_geo_ts.write(tmp_dir +'/ts.xyz')
        M = Molecule(tmp_list[i]+'/chain_0000.xyz')
        reac = M[0]
        prod = M[-1]
        M_qc_reac, M_geo_reac, E_reac = wf.optimize(reac, charge=charge, mult=mult, method=tsmethod, basis=tsbasis)
        M_geo_reac.write(tmp_dir +'/reactant.xyz')
        M_qc_prod, M_geo_prod, E_prod = wf.optimize(prod, charge=charge, mult=mult, method=tsmethod, basis=tsbasis)
        M_geo_prod.write(tmp_dir + '/product.xyz')
        x = ['Reactant', 'Transition', 'Product']
        y = [0, (E_ts - E_reac)*au2kcal, (E_prod-E_reac)*au2kcal]
        fig, ax = plt.subplots()
        ax.set_title('Electronic Energy Differences of %s refinement result' %inp, size = 12)
        ax.plot(x, y, marker = '_', markersize = 50, linestyle='dotted')
        ax.set_ylabel('Energy (kcal/mol)', size = 12)
        ax.text(0.5, y[1], 'Ea=' + str(np.round(y[1], 1)) + 'kcal/mol', size = 12)
        ax.text(1.5, y[-1], 'Ea=' + str(np.round(y[-1], 1)) + 'kcal/mol', size = 12)
        fig.savefig(tmp_dir+'/result.png') 
    print("All the detected reactions were optimized.")
         
 
if __name__ == '__main__':
    main()














