#!/usr/bin/env python

"""
qcaw.py

"""
import os, sys, subprocess, time, shutil, json
import qcengine as qcng
import qcelemental as qcel
import qcportal as ptl
from .molecule import Molecule, EqualSpacing, TopEqual, MolEqual
from collections import OrderedDict
from .errors import NoDatasetError, InvalidCommandError, QCFractalError, QCEngineError
from .constant import bohr2ang
import numpy as np

#def load_xyz(initial, subsample=None):
#    '''
#    This function will accept an xyz file and return geomeTRIC molecule ohjects with subsamples.    
#
#    Parameters
#    ----------
#    initial : string
#        xyz file name 
#            
#    Return
#    ---------
#    frames : molecule object(s)
#        An array of molecule objects
#    '''
#    if subsample != None:
#        subsample = int(subsample)
#    M = Molecule(initial, topframe = ) #add frame numbers here
#    frames = M[::subsample]
#    return frames

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


def qc_to_geo(qc_M, comment="", b2a=False):
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

class User:
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
        for line in info:
            if 'host' in line:
                host = line.strip().split(' ')[-1]
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

class Dataset:
    """
    This class helps users with handling dataset. 
    """
    def __init__(self, name, ds_type, client):
        """
        Parameters
        ----------
        name : str
            Name of the dataset
        
        ds_type : str
            Type of the dataset. Currently 'OptimizationDataset' and 'Dataset' are supported.

        clinet : client object from the User class 
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
                print ("\'%s\' dataset has been created." %self.name)
            except:
                ds = self.client.get_collection(self.ds_type, name = self.name)
        elif command == 'load':
            try:
                ds = self.client.get_collection(self.ds_type, name = self.name)
                print("\'%s\' dataset has been loaded." %self.name)
            except:
                raise NoDatasetError("\'%s\' dataset can't be loaded since there is no %s named %s." % (self.name, self.ds_type, self.name))
        elif command == 'delete':
            try:
                self.client.delete_collection(self.ds_type, name = self.name)
                ds = None
                print ("\'%s\' dataset has been deleted." %self.name)
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
    
       

class Workflow:
    """
    This class perfroms the automated reaction refinement workflow using QCArchive Infrastructure
    """
    def __init__(self, client=False, ds=None, spec_name=None, initial=None, charge=None, mult=None):
        """
        Parameters
        ----------
        ds : Dataset object
            Dataset object created from the Dataset class

        spec_name : string
            QCSpecification name            

        client : client object

        initial : string
            xyz file name 

        charge : int
            molecular charge.

        mult : int
            molecular multiplicity.
        """
        
        self.initial = initial
        self.ds = ds
        self.client = client
        self.spec_name = spec_name
        self.M = Molecule(initial)
        self.charge = charge 
        self.mult = mult

    def resubmit(self): 
        """
        This function detects ERROR calculation results in a given dataset with a specficiation and submit them again. 
        """
        #opt = self.ds.status(self.spec_name, collapse = False)
        opts = self.ds.df[self.spec_name].tolist()
        num = 0
        for i in range (len(opts)):
            if opts[i].error != None:
                num += 1
                self.client.modify_tasks('restart', opts[i].id)
        print ("%i failed jobs in \'%s\' dataset with \'%s\' specfication have been submitted again." %(num, self.ds.name, self.spec_name))

    def equal(self, m1, m2):
        return TopEqual(m1, m2) if True else MolEqual(m1, m2)

    def energy(self, method=None, basis=None, compute=False):
        """
        This function converts xyz file to QCAI molecule objects and creates single point energy calculation jobs.          

        Parameters
        -----------
        method, basis : string
            Electron structure method and basis sets

        compute : boolean
            'compute = False' will only save the molecules in the given dataset and specfication. 'compute = True' will submit the jobs to server.

        spec_overwrite : boolean
            'spec_overwrite = True' will overwrite the spec if the same name spec exists

        Return
        ----------
        ds_sp : dataset object
            print(ds_opt.status(collapse = False)) will show the dataset's status 

        """  
        if method == None:
            method = 'b3lyp'
        if basis == None:
            basis = '6-31g(d)' 
        ds_sp = self.ds 
        M = self.M 
        key = ptl.models.KeywordSet(values = {'maxiter':1000,
                                             'e_convergence': 1e-6,
                                             'guess' : 'sad',
                                             'scf_type' : 'df'})        
        ds_sp.add_keywords(self.spec_name, 'psi4', key)

        spec = {
                'program' : 'psi4',
                'method' : method,
                'basis' : basis,
                'keywords': self.spec_name,
                'tag' : None}
        frames = list(range(len(M)))

        for frm in frames:               
            mol = qcel.models.Molecule(**{'symbols': M[frm].elem, 'geometry': np.array(M[frm].xyzs)/0.529177210, 'molecular_charge' : self.charge, 'molecular_multiplicity' : self.mult}) 
            #if len(M[frm].comms[0]) == 0 or len(M[frm].comms[0]) > 10:
            #    raise RuntimeError('Please provide a short name, less than 10 letters, of the molecule in comment line in the input xyz file (between number of atoms and coordinates).')
            try:
                ds_sp.add_entry('%s_%s' %(formula(M[frm])+str(frm)), mol)
            except:
                pass 
        ds_sp.save()
        if compute:
            ds_sp.compute(**spec)
            print ("Single point energy calculations in %s have been submitted. Run the QCFractal manager to carry the calculations." %(ds_sp.name))
        return ds_sp 
        

    def optimization(self, method='b3lyp', basis='6-31g(d)', subsample=10, maxiter=500, compute=False, spec_overwrite=False): 
        """
        This function converts the xyz file to QCAI molecule objects and sets up optimization jobs.          

        Parameters
        -----------
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

        """  
        ds_opt = self.ds
        key = [ptl.models.KeywordSet(values = {'maxiter':maxiter})]        
        key_id = self.client.add_keywords(key)[0]
        optimize = {
            'name' : self.spec_name,
            'optimization_spec' : {'program': 'geometric', 'keywords': None},
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
        except:
            pass

        M = self.M 
        frames = list(range(0, len(M), subsample))
        if (len(M)-1) not in frames:
            frames.append(len(M)-1)
        for frm in frames:               
            mol = qcel.models.Molecule(**{'symbols': M[frm].elem, 'geometry': np.array(M[frm].xyzs)/bohr2ang,  'molecular_charge' : self.charge, 'molecular_multiplicity' : self.mult}) 
            try:
                ds_opt.add_entry('%s_%i' %(self.initial.split('.')[0], frm), mol, save = False)
            except:
                pass 
        ds_opt.save()
        if compute:
            ds_opt.compute(self.spec_name)
            print ("Calculations in \'%s\' with \'%s\' specification have been submitted. Run the QCFractal manager to carry the calculations." %(ds_opt.name, self.spec_name))
        return ds_opt

    
    def smoothing(self):
        '''
        Once the optimization is done, smoothing function will detect reactions and smooth them for the NEB calculation. 
        '''
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
            qcel_M = self.ds.get_record(name, self.spec_name).get_final_molecule()
            OptMols[frm] = qc_to_geo(qcel_M, b2a = True)  
            Iter[frm] = len(self.ds.get_record(name, self.spec_name).trajectory)

        if err > 0 : 
            print ("WARNING: %i ERROR or/and INCOMPLETE results detected." % err) 

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
        self.neb_inputs = {}

        for i in range(len(MolPairs)): 
            (a,b) = FramePairs[i][np.argmin([(jb-ja) for (ja, jb) in FramePairs[i]])]
            qc_mol_Traj1 = self.ds.get_record(mol_name + '_' + str(a), self.spec_name).get_molecular_trajectory()
            qc_mol_Traj2 = self.ds.get_record(mol_name + '_' + str(b), self.spec_name).get_molecular_trajectory()

            for j in range(len(qc_mol_Traj1)-1):
                if geo_mol_Traj == None:
                    geo_mol_Traj = qc_to_geo(qc_mol_Traj1[-1], b2a = True)
                geo_mol_Traj += qc_to_geo(qc_mol_Traj1[::-1][j+1], b2a = True) 
            geo_mol_Traj += self.M[a:b]
            for k in range(len(qc_mol_Traj2)):
                geo_mol_Traj += qc_to_geo(qc_mol_Traj2[k], b2a = True)    
            
            fnum =  str(a) + '-' + str(b)
            fname = str(mol_name +'_'+ fnum)
            path = './%s/' %(mol_name)
            
            NEB_path = path + fnum
            os.mkdir(NEB_path)
            geo_mol_Traj.write(os.path.join(path + fnum,'connected_%s.xyz' %fname))
            equal = EqualSpacing(geo_mol_Traj, dx = 0.05) 
            equal.write(os.path.join(path + fnum, 'spaced_%s.xyz' %fname))
            geo_mol_Traj = None 
            command ='Nebterpolate.py --morse 1e-2 --repulsive --allpairs --anchor 2 %s/spaced_%s.xyz %s/NEB_ready_%s.xyz &> %s/interpolate_%s.log' %(NEB_path, fname, NEB_path, fname, NEB_path, fname)
            log = open('%s/interpolate_%s.log' %(NEB_path, fname), 'a')
            #err = open('%s/interpolate_%s.log' %(NEB_path, fname), 'a')
            subprocess.Popen(command, shell = True, stdout = log, stderr = log)
            self.neb_inputs[i] = NEB_path + '/NEB_ready_%s.xyz' %fname
        print("Smoothing Procedure is running on the local machine. NEB ready xyz files will be generated once the smoothing procedure is done.")
        

    def neb(self, neb, method='b3lyp', basis='6-31+g(d)', images=21, coordsys='cart', ew = False, nebk = 1, avgg=0.025, maxg=0.05, guessk=0.01, guessw=0.5):
        """
        This function will run NEB calculations to locate rough transition state structures.          

        Parameters
        -----------
        neb : str
            Name of the neb ready xyz file.

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

        nebk : float
            Spring constant in kcal/mol/Ang^2

        avgg, maxg : float
            Average RMS-gradient and max RMS-gradient for convergence. Unit of eV/Ang

        guessk, guessw : float
            Guess Hessian and guess weight for chain coordinates.

        Return
        ----------
        it will write each iterated chains in .tmp directory and final transition state geometries.

        """  
        band = Molecule(neb) 

        qcel_mol = qcel.models.Molecule(**{'symbols': band.elem, 'geometry': np.array(band[0].xyzs)/bohr2ang,  'molecular_charge' : self.charge, 'molecular_multiplicity' : self.mult}) 

        neb_input = {
            'keywords' : {
                'program' : 'psi4'
                },
            'input_specification':{
                'driver': 'gradient',
                'model' : {
                    'method': method,
                    'basis': basis
                    },
                'keywords': {
                    'images': images,
                    'avgg' : avgg,
                    'maxg' : maxg,
                    'nebk' : nebk, 
                    'maxcyc' : 200,
                    'guessk' : guessk,
                    'guessw' : guessw,
                    'coords': neb,
                    'coordsys': coordsys,
                    'engine': 'qcengine',
                    'client': self.client, 
                    'input': 'chains'
                    }
                },
            'initial_molecule':qcel_mol
            }

        if ew:
            neb_input['input_specification']['keywords']['ew'] = None

        neb_procedure = qcng.compute_procedure(neb_input, 'geometric')


    def tsopt(self, tsxyz, method='b3lyp', basis='6-31g(d)', maxiter=500):
        """
        This function will run transition           

        Parameters
        -----------
        tsxyz : str
            Name of the xyz file containing rough guess ts structure.

        method, basis : string
            Electron structure method and basis sets

        maxiter : int
            Maximum iteration number for scf calculations.

        Return
        ----------
        M : QCARW Molecule object of the optimized TS structure.
        energy : energy of the optimized TS structure

        """  
        ts = Molecule(tsxyz) 

        qcel_mol = qcel.models.Molecule(**{'symbols': ts.elem, 'geometry': np.array(ts.xyzs)/bohr2ang,  'molecular_charge' : self.charge, 'molecular_multiplicity' : self.mult}) 
        if not self.client:
            """
            TS optimization will be carried locally
            """
            tsopt_input = {
                'keywords' : {
                    'program' : 'psi4'
                    },
                'input_specification':{
                    'driver': 'gradient',
                    'model' : {
                        'method': 'TS-'+ method,
                        'basis': basis
                        }
                    
                    },
                'initial_molecule':qcel_mol
                }


            tsopt_result = qcng.compute_procedure(tsopt_input, 'geometric') #OptimizationResult
            if type(tsopt_result).__name__ == "FailedOperation":
                raise QCEngineError(tsopt_result.error)
            else:
                M = qc_to_geo(tsopt_result.final_molecule, b2a=True)
                energy = tsopt_result.energies[-1]
        else:
            """
            TS optimization will be carried with QCFractal server
            """
            import time
            #Adding transition in qc_spec will cause the internal server error.
            #Adding transistion in the 1st line or in qc_spec keywords blocks interaction between server and workers.  
            #Adding 'transition': 'yes' in qc_spec does not work either..
            tsopt_qcschema = {
                 'keywords' : None,
                 'qc_spec': {
                    'driver': 'gradient',
                    'method': "TS-" + method,
                    'basis': basis,
                    'program': 'psi4'
                    }
                }


            r = self.client.add_procedure('optimization', 'geometric', tsopt_qcschema, [qcel_mol])
            proc_id = r.ids
            loop = 0 
            while True:
                proc = self.client.query_procedures(id=proc_id)[0] #OptimizationRecord
                status = proc.status.split('.')[-1].upper().strip()
                print(status)
                if status == 'INCOMPLETE':
                    time.sleep(50)
                    loop += 1
                elif status == 'ERROR':
                    print('Error detected')
                    res = self.client.modify_tasks('restart',proc.id)
                    print(res.n_updated,"ERROR status optimization resubmitted")
                    loop += 1
                elif status == 'COMPLETE':
                    M = qc_to_geo(proc.get_final_molecule(), b2a=True)     
                    energy = proc.get_final_energy()
                    print("QCAI Optimization is done.")
                    break
 
        return M, energy
        
       # proc_id = tsopt_procedure.ids
       # loop = 0 
       # while True:
       #     proc = self.client.query_procedures(id=proc_id)[0] #OptimizationRecord
       #     status = proc.status.split('.')[-1].upper().strip()
       #     print(status)
       #     if status == "INCOMPLETE":
       #         if loop == 0:
       #             self.client.modify_tasks("restart",proc.id)
       #             loop += 1
       #         else:
       #             time.sleep(50)
       #             loop += 1
       #     elif status == "ERROR":
       #         print("Error detected")
       #         res = self.engine.client.modify_tasks("restart",proc.id)
       #         print(res.n_updated,"ERROR status optimization resubmitted")
       #         loop += 1
       #     elif status == "COMPLETE":
       #         qc_mol = proc.get_final_molecule()     
       #         energy = proc.get_final_energy()
       #         print("QCAI Optimization is done.")
       #         break
       # M = qc_to_geo(qc_mol, energy, b2a=True)
       # M.write("ts_opt.xyz")
        #return energy, M
        
        
        














