#!/usr/bin/env python

"""
qcaw.py

"""
import os, sys, subprocess, time
import qcengine as qcng
import qcelemental as qcel
import qcportal as ptl
from QCARWorkflow.molecule import Molecule, EqualSpacing, TopEqual, MolEqual
from collections import OrderedDict
import numpy as np

#def load_xyz(initial, subsample=None):
#    """
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
#    """
#    if subsample != None:
#        subsample = int(subsample)
#    M = Molecule(initial, topframe = ) #add frame numbers here
#    frames = M[::subsample]
#    return frames



def qc_to_geo(qc_M, b2a = False):
    """
    Convert QCArchive molecule object to geomeTRIC molecule objects
   
    Parameters
    ----------
    qc_M : QCArchive molecule object

    b2a : boolean 
        b2a = True will convery length unit from Bohr to Angstrom

    Return
    ----------
    geo_M : geomeTRIC molecule object
    """
    geo_M = Molecule()
    geo_M.comms = ["placeholder"]
    geo_M.elem = list(qc_M.symbols)
    geom = np.array(qc_M.geometry, dtype=np.float64).reshape(-1, 3)
    if b2a == True:
        geom *= 0.529177210 
    geo_M.xyzs = [geom]
    geo_M.build_bonds()
    geo_M.build_topology()
    return geo_M



def resubmit_all(client):
        """
        This function will submit ALL the failed jobs again to server.
        """
        errors = client.query_tasks(status="ERROR")
        for i in range(len(errors)):
            client.modify_tasks("restart", errors[i].base_result)
        print ("All the failed calculations in the server were submitted again.")


class User:
    def __init__(self, user, password):
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
        
        Parameters
        ----------
        start : boolean
            start = True will start the server. As long as the running server does not get killed, start = False will work.

        Return
        ----------
        client : client object 
            With the client object, a user can create/access datasets
        """
        info = os.popen("qcfractal-server info").readlines()
        for line in info:
            if "host" in line:
                host = line.strip().split(" ")[-1]
            if "port" in line:
                port = int(line.strip().split(" ")[-1])
        try:
            client = ptl.FractalClient("%s:%i"%(host, port), verify=False, username = self.user, password= self.password)
        except:
            for i in range(0, 3):
                try:
                    subprocess.run("qcfractal-server start &", shell = True, stdout=subprocess.DEVNULL)
                    time.sleep(3.0)
                    client = ptl.FractalClient("%s:%i"%(host, port), verify=False, username = self.user, password= self.password)
                    cnt_error = None
                except Exception as cnt_error:
                    pass
                if cnt_error == None:
                    break 
            print("Server is running") 
        print("Client is ready")
       
        return client           

class Dataset:
    def __init__(self, name, ds_type, client):
        """
        Parameters
        ----------
        name : str
            Name of the dataset
        
        ds_type : str
            Type of the dataset. Currently "OptimizationDataset" and "Dataset" are supported.

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
            1. "make" will create a new dataset with a given name.
            2. "load" will load the dataset with a given name.
            3. "delete" will delete the dataset with a given name.
            4. "reset" will delete and re-create the dataset with a give name.

        """
        if command == "make":
            try:
                if self.ds_type == "OptimizationDataset":
                    new_ds = ptl.collections.OptimizationDataset(name = self.name, client = self.client)            
                elif self.ds_type == "Dataset":
                    new_ds = ptl.collections.Dataset(name = self.name, client = self.client)    
                new_ds.save()        
                ds = self.client.get_collection(self.ds_type, name = self.name)
                print ("\'%s\' dataset has been created." %self.name)
            except:
                ds = self.client.get_collection(self.ds_type, name = self.name)
        elif command == "load":
            try:
                ds = self.client.get_collection(self.ds_type, name = self.name)
                print("\'%s\' dataset has been loaded." %self.name)
            except:
                raise RuntimeError("\'%s\' dataset can't be loaded since there is no %s named %s." % (self.name, self.ds_type, self.name))
        elif command == "delete":
            try:
                self.client.delete_collection(self.ds_type, name = self.name)
                ds = None
                print ("\'%s\' dataset has been deleted." %self.name)
            except:
                raise RuntimeError("\'%s\' dataset can't be deleted since there is no %s named %s." % (self.name, self.ds_type, self.name))
        elif command == "reset":
            try:
                self.client.delete_collection(self.ds_type, name = self.name)
                if self.ds_type == "OptimizationDataset":
                    ds = ptl.collections.OptimizationDataset(name = self.name, client = self.client)
                elif self.ds_type == "Dataset":
                    ds = ptl.collections.Dataset(name = self.name, client = self.client) 
            except:
                ds = ptl.collections.OptimizationDataset(name = self.name, client = self.client)
        else:
           raise RuntimeError("Please provide a valid command for the dataset.")  
        return ds
    
       

class Workflow:
    """
    This class handles the automated workflow using QCArchive Infrastructure
    """
    def __init__(self, ds, spec_name, client, initial):
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

        """
        
        self.initial = initial
        self.ds = ds
        self.client = client
        self.spec_name = spec_name
        self.M = Molecule(initial)
 
    def energy(self, method=None, basis=None, compute = False):
        """
        This function converts xyz file to QCAI molecule objects and creates single point energy calculation jobs.          

        Parameters
        -----------
        method, basis : string
            Electron structure method and basis sets

        compute : boolean
            "compute = False" will only save the molecules in the given dataset and specfication. "compute = True" will submit the jobs to server.

        spec_overwrite : boolean
            "spec_overwrite = True" will overwrite the spec if the same name spec exists

        Return
        ----------
        ds_sp : dataset object
            print(ds_opt.status(collapse = False)) will show the dataset's status 

        """  
        if method == None:
            method = "b3lyp"
        if basis == None:
            basis = "6-31g(d)" 
        ds_sp = self.ds 
        M = self.M 
        key = ptl.models.KeywordSet(values={'maxiter':1000,
                                             'e_convergence': 1e-6,
                                             'guess' : 'sad',
                                             'scf_type' : 'df'})        
        ds_sp.add_keywords(self.spec_name, "psi4",key)

        spec = {
                "program" : "psi4",
                "method" : method,
                "basis" : basis,
                "keywords": self.spec_name,
                "tag" : None}
        frames = list(range(len(M)))
        for frm in frames:               
            mol = qcel.models.Molecule(**{"symbols": M[frm].elem, "geometry": M[frm].xyzs}) 
            if len(M[frm].comms[0]) == 0 or len(M[frm].comms[0]) > 10:
                raise RuntimeError("Please provide a short name, less than 10 letters, of the molecule in comment line in the input xyz file (between number of atoms and coordinates).")
            try:
                ds_sp.add_entry("%s" %(M[frm].comms[0]), mol)
            except:
                pass 
        ds_sp.save()
        if compute:
            ds_sp.compute(**spec)
            print ("Single point energy calculations in %s have been submitted. Run the QCFractal manager to carry the calculations." %(ds_sp.name))
        return ds_sp 
        

    def optimization(self, method=None, basis=None, subsample = None, compute = False, maxiter = None, spec_overwrite = False): 
        """
        This function converts the xyz file to QCAI molecule objects and sets up optimization jobs.          

        Parameters
        -----------
        method, basis : string
            Electron structure method and basis sets

        subsample : int
            Frame interval for subsampling trajectories

        compute : boolean
            "compute = False" will only save the molecules in the given dataset and specfication. "compute = True" will submit the jobs to server.

        maxiter : int
            maximum scf iteration number

        spec_overwrite : boolean
            "spec_overwrite = True" will overwrite the spec if the same name spec exists

        Return
        ----------
        ds_opt : dataset object
            print(ds_opt.status(collapse = False)) will show the dataset's status 

        """  
        if method == None:
            method = "b3lyp"
        if basis == None:
            basis = "6-31g(d)"   
        if subsample == None:
            subsample = 10
        if maxiter == None:
            maxiter = 500
        ds_opt = self.ds
        key = [ptl.models.KeywordSet(values={'maxiter':maxiter})]        
        key_id = self.client.add_keywords(key)[0]
        optimize = {
            "name" : self.spec_name,
            "optimization_spec" : {"program": "geometric", "keywords": None},
            "qc_spec" : {
                    "driver": "gradient", 
                    "method": method, 
                    "basis": basis,
                    "keywords": key_id,
                    "program": "psi4"
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
            mol = qcel.models.Molecule(**{"symbols": M[frm].elem, "geometry": M[frm].xyzs}) 
            try:
                ds_opt.add_entry("%s_%i" %(self.initial.split(".")[0], frm), mol, save = False)
            except:
                pass 
        ds_opt.save()
        if compute:
            ds_opt.compute(self.spec_name)
            print ("Calculations in \'%s\' with \'%s\' specification have been submitted. Run the QCFractal manager to carry the calculations." %(ds_opt.name, self.spec_name))
        return ds_opt

    def resubmit(self): 
        """
        This function detects ERROR calculation results in a given dataset with a specficiation and submit them again. 
        """
        opt = ds.status(self.spec_name, collapse = False)
        opts = ds.df[self.spec_name].tolist()
        num = 0
        for i in range (len(opts)):
            if opts[i].error != None:
                num += 1
                self.client.modify_task("restart", opts[i].id)
        print ("%i failed jobs in %s dataset with %s specfication have been submitted again." %(num, self.ds.name, self.spec_name))

    def Equal(self, m1, m2):
        return TopEqual(m1, m2) if True else MolEqual(m1, m2)

    def smoothing(self):
        """
        Once the optimization is done, smoothing function will detect reactions and smooth them for NEB calculation.
        """
        opt = self.ds.status(self.spec_name, collapse = False)
        opts = self.ds.df[self.spec_name].tolist()
        stats = [opt.status for opt in opts] 
        mol_names = opt.index.tolist() 
        mol_name = '_'.join(str(elem) for elem in mol_names[0].split("_")[:-1])
       # for stat in stats:
       #     if stat.value == "ERROR" or stat.value == "INCOMPLETE":
       #         raise RuntimeError ("ERROR or INCOMPLETE result detected. Make sure all the optimizations in \'%s\' dataset with \'%s\' specification are completed." %(self.ds.name, self.spec_name))
        
        OptMols = OrderedDict()
        for name, calc in zip(mol_names, stats):
            frm = int(name.split("_")[-1])
            err = 0
            if calc.value == "ERROR" or calc.value == "INCOMPLETE":
                err += 1
                if err == 0:
                    print ("WARNING: ERROR or INCOMPLETE result detected.") 
                continue
            qcel_M = self.ds.get_record(name, self.spec_name).get_final_molecule()
            OptMols[frm] = qc_to_geo(qcel_M, b2a = True)  

        path_initial = [] 
        path_final = []
        for fi, fj in zip(list(OptMols.keys())[:-1], list(OptMols.keys())[1:]): 
            if not self.Equal(OptMols[fi], OptMols[fj]): 
                path_initial.append(fi)
                path_final.append(fj)

        MolPairs = []
        FramePairs = []
            
        for fi in path_initial:
            for fj in path_final:
                if fj > fi and (not self.Equal(OptMols[fi], OptMols[fj])):
                    if (fj - fi) > 1000: continue
                    NewPair = True
                    for imp, (m1, m2) in enumerate(MolPairs):
                        if self.Equal(OptMols[fi], m1) and self.Equal(OptMols[fj], m2):
                            FramePairs[imp].append((fi, fj))
                            NewPair = False
                            break
                        elif self.Equal(OptMols[fi], m2) and self.Equal(OptMols[fj], m1):
                            FramePairs[imp].append((fi, fj))
                            NewPair = False
                            break
                    if NewPair:
                        MolPairs.append((OptMols[fi], OptMols[fj]))
                        FramePairs.append([(fi, fj)])

        if len(MolPairs) != len(FramePairs) or len (MolPairs) == 0:
            raise RuntimeError ("No reactions are detected or the Number of detected pairs of reacting molecules and frames don't match.")            

        geo_mol_Traj = None
        for i in range(len(MolPairs)): 
            (a,b) = FramePairs[i][np.argmin([(jb-ja) for (ja, jb) in FramePairs[i]])]
            qc_mol_Traj1 = self.ds.get_record(mol_name + "_" + str(a), self.spec_name).get_molecular_trajectory()
            qc_mol_Traj2 = self.ds.get_record(mol_name + "_" + str(b), self.spec_name).get_molecular_trajectory()

            for j in range(len(qc_mol_Traj1)-1):
                if geo_mol_Traj == None:
                    geo_mol_Traj = qc_to_geo(qc_mol_Traj1[-1])
                geo_mol_Traj += qc_to_geo(qc_mol_Traj1[::-1][j+1]) 
            geo_mol_Traj += self.M[a:b]
            for k in range(len(qc_mol_Traj2)):
                geo_mol_Traj += qc_to_geo(qc_mol_Traj2[k])    
            
            fnum =  str(a) + "-" + str(b)
            fname = str(mol_name +"_"+ fnum)
            path = "./%s/" %(mol_name)
            
            NEB_path = path + fnum
            if not os.path.exists(path):
                os.mkdir(path)
            os.mkdir(NEB_path)
            geo_mol_Traj.write(os.path.join(path + fnum,"connected_%s.xyz" %fname))
            equal = EqualSpacing(geo_mol_Traj, dx = 0.05) 
            equal.write(os.path.join(path + fnum, "spaced_%s.xyz" %fname))
            geo_mol_Traj= None 
            command ="Nebterpolate.py --morse 1e-2 --repulsive --allpairs --anchor 2 %s/spaced_%s.xyz %s/NEB_ready_%s.xyz &> %s/interpolate_%s.log" %(NEB_path, fname, NEB_path, fname, NEB_path, fname)
            log = open('%s/interpolate_%s.log' %(NEB_path, fname), 'a')
            #err = open('%s/interpolate_%s.log' %(NEB_path, fname), 'a')
            subprocess.Popen(command, shell=True, stdout = log, stderr = log)
        print("Smoothing Procedure is running on the local machine. NEB ready xyz files will be generated once the smoothing procedure is done.")
        












