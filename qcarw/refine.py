#!/usr/bin/env python

"""
refine.py

"""
import os, sys, subprocess, time, shutil, socket
import qcengine as qcng
import qcelemental as qcel
import qcportal as ptl
import matplotlib.pyplot as plt
import numpy as np
from .params import parse_refine_args
from .molecule import PeriodicTable, Molecule, EqualSpacing, TopEqual, MolEqual
from collections import OrderedDict, Counter
from .errors import NoDatasetError, InvalidCommandError, QCFractalError, QCEngineError, OptimizeInputError
from .constant import bohr2ang, au2kcal

def find_groups(sl1, sl2):
    """ 
    [This function is copied from https://github.com/leeping/nanoreactor/src/rxndb.py]
    Given two lists of atom lists, find the groups of sets in each
    list that only contain each others' elements (i.e. if somehow we
    have two parallel reactions in one.)
    
    Parameters
    ----------
    sl1, sl2 : list of lists
        List of lists of atom indices corresponding to molecules.

    Returns
    -------
    list of lists
        Atom indices corresponding to separately reacting groups of molecules.
    """
    # Convert to sets
    sl1c = [set(s) for s in sl1]
    sl2c = [set(s) for s in sl2]
    # Iterate this while loop until we find a single set of atom groups
    while set([tuple(sorted(list(s))) for s in sl1c]) != set([tuple(sorted(list(s))) for s in sl2c]):
        # Double loop over molecule atom indices
        for s1 in sl1c:
            for s2 in sl2c:
                # For any pair of molecules that have any overlapping atoms,
                # add all atoms in each atom set to the other.
                if len(s1.intersection(s2)) > 0:
                    s1.update(s2)
                    s2.update(s1)
    result = sorted([list(t) for t in list(set([tuple(sorted(list(s))) for s in sl1c]))])
    return result

def find_reacting_groups(m1, m2):
    """
    [This function is copied from https://github.com/leeping/nanoreactor/src/rxndb.py]
    Given two Molecule objects, determine the groups of atoms that
    reacted with each other (i.e. formed different molecules.)  This will
    remove spectator atoms (ones that didn't react at all) and separate 
    concurrent reactions occuring in different places.

    Parameters
    ----------
    m1, m2 : Molecule
        Length-1 Molecule objects corresponding to reactant and product
        frames.  For the sake of future electronic structure calculations,
        these objects must have qm_mulliken_charges and qm_mulliken_spins.

    Returns
    -------
    extracts: list of 3-tuple of list, int, int
        Each 3-tuple is a group of atoms that reacted, and their associated
        charge / multiplicity.
    """
    if not isinstance(m1, Molecule) or len(m1) != 1:
        raise RuntimeError("Please only pass length-1 Molecule objects")
    if not isinstance(m2, Molecule) or len(m2) != 1:
        raise RuntimeError("Please only pass length-1 Molecule objects")


    # Get a list of atom indices belonging to each molecule.
    m1_mol_atoms = [g.L() for g in m1.molecules]
    m2_mol_atoms = [g.L() for g in m2.molecules]
    # Count the number of atoms in spectator molecules that don't
    # react at all, and store their molecular formulas (for public
    # shaming).
    n_spectator_atoms = 0
    spectator_formulas = []
    strrxns = []
    # The results: extract groups of atoms to extract corresponding to
    # individual reaction pathways, and the net charge / multiplicity
    # belonging to 
    extract_groups = []
    extract_charges = []
    extract_mults = []
    # Separate atoms into groups of separately reacting molecules.
    # This is also effective at finding spectator molecules that don't react at all.
    do_extract = False
    for atom_group in find_groups(m1_mol_atoms, m2_mol_atoms):
        m1g = m1.atom_select(atom_group)
        m2g = m2.atom_select(atom_group)
        spectator_atoms = []
        print("atom group: %s" % str(atom_group))
        print("m1 molecules: %s" % str([[atom_group[i] for i in g.L()] for g in m1g.molecules]))
        print("m2 molecules: %s" % str([[atom_group[i] for i in g.L()] for g in m2g.molecules]))
        for g1 in m1g.molecules:
            print("atoms in molecule: %s" % str(g1.L()))
            for g2 in m2g.molecules:
                # Graphs are usually compared by comparing elements and
                # topology, but a spectator molecule also has the same
                # atom numbers.
                if g1 == g2 and g1.L() == g2.L():
                    spectator_atoms += g1.L()
        # Since we already separated the atoms into groups of separately reacting ones,
        # any atom group with spectator atoms is expected to be a single spectator molecule.
        if len(spectator_atoms) > 0:
            print("spectator atoms: %s" % str([atom_group[i] for i in spectator_atoms]))


        if len(spectator_atoms) == m1g.na:
            if len(m1g.molecules) != 1:
                raise RuntimeError("I expected an atom group with all spectators to be a single molecule")
            n_spectator_atoms += len(spectator_atoms)
            spectator_formulas.append(m1g.molecules[0].ef())
            continue
        elif len(spectator_atoms) > 0:
            raise RuntimeError("I expected an atom group with any spectators to be a single molecule")
        else:
            strrxn = ' + '.join(['%s%s' % (str(j) if j>1 else '', i) for i, j in list(Counter([m.ef() for m in m1g.molecules]).items())])
            strrxn += ' -> '
            strrxn += ' + '.join(['%s%s' % (str(j) if j>1 else '', i) for i, j in list(Counter([m.ef() for m in m2g.molecules]).items())])
            strrxns.append(strrxn)
 
        # Now we have a group of reacting atoms that we can extract from the
        # pathway, but we should perform some sanity checks first.
        mjoin = m1g + m2g
        # A bit of code copied from extract_pop.  Verify that the reacting
        # atoms have consistent charge and spin in the two passed Molecule
        # objects.  If not consistent, then we cannot extract spectator atoms.
        Chgs = np.array([sum(i) for i in mjoin.qm_mulliken_charges])
        SpnZs = np.array([sum(i) for i in mjoin.qm_mulliken_spins])
        chg, chgpass = extract_int(Chgs, 0.3, 1.0, label="charge")
        spn, spnpass = extract_int(abs(SpnZs), 0.3, 1.0, label="spin-z")
        nproton = sum([Elements.index(i) for i in m1g.elem])
        nelectron = nproton + chg
        # If the sanity checks fail, then do not extract the spectator atoms
        # and simply return a list of all the atoms at the end.
        do_extract = True
        if ((nelectron-spn)//2)*2 != (nelectron-spn):
            print("\x1b[91mThe number of electrons (%i; charge %i) is inconsistent with the spin-z (%i)\x1b[0m" % (nelectron, chg, spn))
            do_extract = False
            break
        if (not chgpass or not spnpass):
            print("\x1b[91mCannot determine a consistent set of spins/charges after extracting spectators\x1b[0m")
            do_extract = False
            break
        extract_groups.append(np.array(atom_group))
        extract_charges.append(chg)
        extract_mults.append(abs(spn)+1)
    if do_extract:
        message = "Initial Reaction : " + ' ; '.join(strrxns)
        if n_spectator_atoms > 0:
            # I know it's supposed to be spelled 'spectator', but it's fun to say 'speculator' :)
            message += " ; Speculators (removed) : \x1b[91m%s\x1b[0m" % (' + '.join(['%s%s' % (str(j) if j>1 else '', i) for i, j in list(Counter(spectator_formulas).items())]))
        print(message)
        return list(zip(extract_groups, extract_charges, extract_mults))
    else:
        print("Unable to split reaction pathway into groups")
        return list(zip([np.arange(m1.na)], [m1.charge], [m1.mult]))

def equal(m1, m2):
    """
    To see whether the two Molecule objects have the same geometry

    Parameters
    ----------
    m1, m2: Molecule Object

    Return
    ----------
    True or False
    """
    return TopEqual(m1, m2) if True else MolEqual(m1, m2)

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
    if b2a:
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
    
        password : password for QCFractal server
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
        client = None
        #info = os.popen('qcfractal-server info').readlines()
        #for line in info:
        #    if 'port' in line:
        #        port = int(line.strip().split(' ')[-1])
        #host = socket.gethostname()
        client = ptl.PortalClient('http://%s:%i'%('localhost', 7777), username = self.user, password = self.password, verify = False)
        #try:
        #    client = ptl.PortalClient('%s:%i'%('localhost', port), username = self.user, password = self.password, verify = False)
        #except:
        #    client = ptl.PortalClient('%s:%i'%(host, port), username = self.user, password = self.password, verify = False)

          #for i in range(3):
          #    try:
          #        subprocess.run('qcfractal-server start &', shell = True, stdout = subprocess.DEVNULL)
          #        time.sleep(3.0)
          #        client = ptl.PortalClient('%s:%i'%(host, port), username = self.user, password = self.password, verify = False)
          #    except:
          #        if client != None:
          #            break 
          #        else:
          #            continue
        if client is None:
            raise QCFractalError("Client could not be claimed properly. Make sure to run the server \'qcfractal-server start\'.")
 
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
                    new_ds = ptl.add_dataset.OptimizationDataset(name = self.name, client = self.client)            
                elif self.ds_type == 'Dataset':
                    new_ds = ptl.add_collections.Dataset(name = self.name, client = self.client)    
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
    def __init__(self, initial, charge=0, mult=1, client=False, ds=None, spec_name=None):
        """
        Parameters
        ----------
        initial : string
            MD trajectory (xyz file name) that needs to be refined

        charge, mult: int
            Molecular charge and multiplicity
        ds : Dataset object
            Dataset object created from the Dataset class

        spec_name : string
            QCSpecification name            

        client : client object
        """
        self.initial = Molecule(initial)
        self.charge = charge
        self.mult = mult
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

    def dsoptimize(self, method='b3lyp', basis='6-31g(d)', subsample=10, maxiter=100, coordsys='tric'): 
        """
        This function converts the xyz file to QCAI molecule objects and sets up optimization jobs.          

        Parameters
        -----------
        method, basis : string
            Electron structure method and basis sets

        subsample : int
            Frame interval for subsampling trajectories

        maxiter : int
            maximum scf iteration number

        Return
        ----------
        ds_opt : dataset object
            print(ds_opt.status(collapse = False)) will show the dataset's status 

        mass : float
            Molecular mass
        """  

        M = self.initial
        ds_opt = self.ds
        key = [ptl.models.KeywordSet(values = {'maxiter':maxiter, 'properties':'mulliken_charges'})]        
        key_id = self.client.add_keywords(key)[0]
        optimize = {
            'name' : self.spec_name,
            'optimization_spec' : {'program': 'geometric', 'keywords': {'coordsys': coordsys}},
            'qc_spec' : {
                    'program': 'psi4',
                    'driver': 'gradient', 
                    'method': method, 
                    'basis': basis,
                    'keywords': key_id
                     }
                       }
        try:
            ds_opt.add_specification(**optimize, overwrite = True)
            print("Specification %s was added into %s" %(self.spec_name, ds_opt))
        except:
            print("Specification %s is either already added or it has key values that are not allowed" %self.spec_name)
            pass

        mass = sum([PeriodicTable.get(M.elem[i], 0.0) for i in range(M.na)])
        frames = list(range(0, len(M), subsample))
        if (len(M)-1) not in frames:
            frames.append(len(M)-1)
        for frm in frames:               
            print('1')
            mol = qcel.models.Molecule(**{'symbols': M[frm].elem, 'geometry': np.array(M[frm].xyzs)/bohr2ang,  'molecular_charge' : self.charge, 'molecular_multiplicity' : self.mult}) 
            print('2')
            try:
    
                print('3')
                ds_opt.add_entries('%s_%i' %(self.initial.split('.')[0], frm), mol, save = False)
            except:
                print('4')
                pass 

        print('Optimization dataset saved')
        ds_opt.save()
        ds_opt.compute(self.spec_name)
        print ("Calculations in \'%s\' with \'%s\' specification have been submitted. Run the QCFractal manager to carry the calculations." %(ds_opt.name, self.spec_name))
        return ds_opt, mass
    
    def smoothing(self):
        """
        Once the optimization is done, smoothing function will detect reactions and smooth them for the NEB calculation. 
        Parameters
        -----------
        """
        M = self.initial
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
            print("extraaa", record.keywords)
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
            if not equal(OptMols[fi], OptMols[fj]): 
                path_initial.append(fi)
                path_final.append(fj)

        MolPairs = []
        FramePairs = []
            
        for fi in path_initial:
            for fj in path_final:
                if fj > fi and (not equal(OptMols[fi], OptMols[fj])):
                    if (fj - fi) > 1000: continue
                    NewPair = True
                    for i, (m1, m2) in enumerate(MolPairs):
                        if equal(OptMols[fi], m1) and equal(OptMols[fj], m2):
                            FramePairs[i].append((fi, fj))
                            NewPair = False
                            break
                        elif equal(OptMols[fi], m2) and equal(OptMols[fj], m1):
                            FramePairs[i].append((fi, fj))
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
        frames = []
        for i in range(len(MolPairs)): 
            (a,b) = FramePairs[i][np.argmin([(jb-ja) for (ja, jb) in FramePairs[i]])]
            qc_mol_Traj1 = self.ds.get_record(mol_name + '_' + str(a), self.spec_name).get_molecular_trajectory()
            qc_mol_Traj2 = self.ds.get_record(mol_name + '_' + str(b), self.spec_name).get_molecular_trajectory()

            for j in range(len(qc_mol_Traj1)-1):
                if geo_mol_Traj == None:
                    geo_mol_Traj = qc_to_geo(qc_mol_Traj1[-1], b2a = True).without('qm_mulliken_charges', 'qm_mulliken_spins')
                geo_mol_Traj += qc_to_geo(qc_mol_Traj1[::-1][j+1], b2a = True).without('qm_mulliken_charges', 'qm_mulliken_spins') 
            geo_mol_Traj += M[a:b]
            for k in range(len(qc_mol_Traj2)):
                geo_mol_Traj += qc_to_geo(qc_mol_Traj2[k], b2a = True).without('qm_mulliken_charges', 'qm_mulliken_spins')    
            
            fnum =  str(a) + '-' + str(b)
            frames.append([a,b])
            #fname = str(mol_name +'_'+ fnum)
            path = './%s/' %(mol_name)
            
            NEB_path = path + fnum
            os.mkdir(NEB_path)
            
            reacting_groups = find_reacting_groups(OptMols[a][-1], OptMols[b][-1]) 
            for rgrp, (ratoms, rcharge, rmult) in enumerate(reacting_groups):
                Joined = geo_mol_Traj.atom_select(ratoms)
                Spaced = EqualSpacing(geo_mol_Traj, dx = 0.05)
                pathname = path + fnum
                Joined.write(os.path.join(pathname, 'connected.xyz'))
                Spaced.write(os.path.join(pathname, 'spaced.xyz'))
        

            geo_mol_Traj = None 
            command ='Nebterpolate.py --morse 1e-2 --repulsive --allpairs --anchor 2 %s/spaced.xyz %s/NEB_ready.xyz &> %s/interpolate.log' %(NEB_path, NEB_path, NEB_path)
            log = open('%s/interpolate.log' %NEB_path, 'a')
            #err = open('%s/interpolate_%s.log' %(NEB_path, fname), 'a')
            subprocess.Popen(command, shell = True, stdout = log, stderr = log)
            neb_inputs[i] = NEB_path + '/NEB_ready.xyz'
        print("Smoothing Procedure is running on the local machine. NEB ready xyz files will be generated once the smoothing procedure is done.")
        return neb_inputs, frames
        

    def neb(self, method='b3lyp', basis='6-31+g(d,p)', images=21, coordsys='cart', ew = False, nebk = 1, avgg=0.025, maxg=0.05, guessk=0.01, guessw=0.5, tmpdir=None):
        """
        This function will run NEB calculations to locate rough transition state structures.          

        Parameters
        -----------
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
        band = self.initial 
        #inp = '.'.join(initial.split('.')[:-1])
        
        qcel_mol = qcel.models.Molecule(**{'symbols': band.elem, 'geometry': np.array(band[0].xyzs)/bohr2ang,  'molecular_charge' : self.charge, 'molecular_multiplicity' : self.mult}) 

        neb_input = {
            'keywords' : {
                'program' : 'psi4',
                'neb' : True,
                'images': images,
                'avgg' : avgg,
                'maxg' : maxg,
                'nebk' : nebk, 
                'maxcyc' : 100,
                'guessk' : guessk,
                'guessw' : guessw,
                'coords': initial,
                'coordsys': coordsys,
                'engine': 'qcengine',
                'client': self.client, 
                'prefix': tmpdir
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
            neb_input['keywords']['nebew'] = 'yes'
        neb_procedure = qcng.compute_procedure(neb_input, 'geometric')
        return neb_procedure

    def optimize(self, initial, charge=0, mult=1, method='b3lyp', basis='6-31+g(d,p)', coordsys='tric', maxiter=500, ts=False):
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

        coordsys : str
            'cart': Cartesian Coordinates
            'prim': Primitive (a.k.a. redundant) Coordinates
            'tric': Translation-Rotational Coordinates
            'dlc' : Delocalized Internal Coordinates
            'hdlc': Hybrid Delocalized Internal Coordinates


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
            M = Molecule(initial) 
        elif isinstance(initial, Molecule):
            M = initial
        else:
            raise OptimizeInputError('Please provide either xyz file name or a geomeTRIC Molecule object.')

        qcel_mol = qcel.models.Molecule(**{'symbols': M.elem, 'geometry': np.array(M.xyzs)/bohr2ang,  'molecular_charge' : charge, 'molecular_multiplicity' : mult}) 
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
                M_geo.qm_energies = [energy]
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
                    M_geo.qm_energies = [energy]
                    print("QCAI Optimization is done.")
                    break
                if loop > 1000:
                    raise QCFractalError('Optimization failed.')
 
        return M_qc, M_geo, energy

    def irc(self, initial, charge=0, mult=1, method='b3lyp', basis='6-31+g(d,p)', coordsys='cart', trust=0.1, tmpdir=None):
        """
        This function will perfrom the IRC method

        Parameters
        -----------
        initial : str or geomeTRIC Molecule object
            Name of the xyz file containing the initial transition structure.

        charge, mult : int
            Molecular charge and multiplicity

        method, basis : string
            Electron structure method and basis sets

        coordsys : str
            'cart': Cartesian Coordinates
            'prim': Primitive (a.k.a. redundant) Coordinates
            'tric': Translation-Rotational Coordinates
            'dlc' : Delocalized Internal Coordinates
            'hdlc': Hybrid Delocalized Internal Coordinates
    
        trust : float
            Trust radius for the IRC iterations

        tmpdir : str
            Temporary diretory for the IRC results

        Return
        ----------
        irc_progress : Molecule object consist of the IRC results <- needs to be modified

        """  
        if isinstance(initial, str):
            M = Molecule(initial) 
        elif isinstance(initial, Molecule):
            M = initial
        else:
            raise OptimizeInputError('Please provide either xyz file name or a geomeTRIC Molecule object.')
        
        qcel_mol = qcel.models.Molecule(**{'symbols': M.elem, 'geometry': np.array(M.xyzs)/bohr2ang,  'molecular_charge' : charge, 'molecular_multiplicity' : mult}) 

        irc_input = {
            'keywords' : {
                'program' : 'psi4',
                'coordsys': coordsys,
                'engine': 'qcengine',
                'client': self.client, 
                'irc': True,
                'trust': trust,
                'prefix': tmpdir
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


        irc_procedure = qcng.compute_procedure(irc_input, 'geometric')
        return irc_procedure

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
    maxiter     = args_dict.get('scf_iter'   , 100)
    optmethod   = args_dict.get('optmethod' , 'b3lyp')
    optbasis    = args_dict.get('optbasis'  , '6-31g(d)')
    optcrdsys   = args_dict.get('optcrdsys' , 'tric')
    tsmethod    = args_dict.get('tsmethod'  , 'b3lyp')
    tsbasis     = args_dict.get('tsbasis'   , '6-31+g(d,p)')
    images      = args_dict.get('images'    , 21)
    nebmethod   = args_dict.get('nebmethod' , 'b3lyp')
    nebbasis    = args_dict.get('nebbasis'  , '6-31+g(d,p)')
    nebcrdsys   = args_dict.get('nebcrdsys' , 'cart')
    nebk        = args_dict.get('nebk'      , 1)
    avgg        = args_dict.get('avgg'      , 0.025)
    maxg        = args_dict.get('maxg'      , 0.05)
    ew          = args_dict.get('ew'        , False)
    #irccrdsys   = args_dict.get('irccrdsys' , 'cart')
    trust       = args_dict.get('ircstep'   , 0.1)
    analyze     = args_dict.get('analyze'   , False)

    client = User(user, password).server()
    
    ds = Dataset(dataset, client).setting('make')
    wf = Workflow(initial=initial, charge=charge, mult=mult, ds=ds, client=client, spec_name=spec_name)
    ds, Mmass =wf.dsoptimize(method=optmethod, maxiter=maxiter,basis=optbasis, subsample=subsample, coordsys=optcrdsys)
    
    cycle = 0
    while True: 
        
        error = 0

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
                error += 1
            if stat == 'COMPLETE':
                comp += 1
        
        
        print("%i/%i calculations are completed" %(comp,num_opt)) 
        wait = (num_opt-comp)*int(Mmass*0.5)

        if comp == num_opt or comp/num_opt> 0.80:
            print('Optimization step is completed.')
            break

        print("Molecular mass = %.2f / Waiting %i seconds" %(Mmass, wait))
        time.sleep(wait) 

        if cycle > 100:
            print("Iteration went over 100. Rest of the procedure will be performed based on the successfully optimized geometries.")
            break

        if comp == 0 and cycle > 10 :
            raise QCFractalError('Jobs are not recognized by QCFractal server. Try to restart the server and manager')
        
        if error/num_opt > 0.2 and cycle > 30:
            print("There are optimizations that could not be converged. Rest of the procedure will be performed based on the successfully optimized geometries.")
            break

        cycle += 1

    smoothed, frames = wf.smoothing()
    print('frames', frames)
    time.sleep(10)
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
                f.write('Smoothing procedure error. It probably just copied the input file to generate result. Check the .log file.')
                break
            else:
                smoothing_cycle += 1
                time.sleep(Mmass*0.5)
        print("NEB method will be used to refine %s initial chain" %neb)
        wf.neb(initial=neb, charge=charge, mult=mult, method=nebmethod, basis=nebbasis, images=images, coordsys=nebcrdsys, ew=ew, nebk=nebk, avgg=avgg, maxg=maxg, tmpdir=inp) 
        print("%s initial chain is processed." %neb) 
        ts_name = inp + '.tsClimb.xyz'
        tmp_name = inp + '.tmp'
        guess_ts_list.append(ts_name)
        tmp_list.append(tmp_name)

    M_info = OrderedDict() #Molecule info
    for i, ts in enumerate(guess_ts_list):
        if not os.path.exists(ts): #Sometimes the smoothing function won't be able to smooth a given rxn path. If there are paths that were not smoothed, it will just skip them.
            continue
        inp = '~/' + '/'.join(('.'.join(tmp_list[i].split('.')[:-1]) + '.xyz').split('/')[-3:-1])
        tmp_dir = '/'.join(tmp_list[i].split('/')[:-1])

        # Optimizing TS structure, geomeTRIC Molecule object has qm_energies attribute. 
        M_qc_ts, M_geo_ts, E_ts= wf.optimize(ts, charge=charge, mult=mult, method=tsmethod, basis=tsbasis, ts=True)
        M_geo_ts.write(os.path.join(tmp_dir,'ts.xyz'))

        wf.irc(initial=M_geo_ts, charge=charge, mult=mult, method=tsmethod, basis=tsbasis, coordsys='cart', trust=trust, tmpdir=tmp_dir) 

        M = Molecule(os.path.join(tmp_dir, 'IRC_%.2f.xyz'%trust))
        reac = M[0]
        prod = M[-1]
        M_qc_reac, M_geo_reac, E_reac = wf.optimize(reac, charge=charge, mult=mult, method=tsmethod, basis=tsbasis)
        M_geo_reac.write(os.path.join(tmp_dir, 'reactant.xyz'))
        M_qc_prod, M_geo_prod, E_prod = wf.optimize(prod, charge=charge, mult=mult, method=tsmethod, basis=tsbasis)
        M_geo_prod.write(os.path.join(tmp_dir, 'product.xyz'))
        
        frame = str(frames[i][0]) + '-' + str(frames[i][-1])
       
        #Saving molecule objects 
        M_info[frame]=[M_geo_reac, M_geo_ts, M_geo_prod]
        
    print("All the detected reactions were optimized.")

    if analyze:
    
        from .connect import connect_rxns 

        print ("Unit reactions will be connected to generate uniqe pathways.")

        cwd = os.getcwd()
        unique_rxns = connect_rxns(M_info)
        for i, (k, v) in enumerate(unique_rxns.items()):
            v.write(os.path.join(cwd,"reaction_%i" %i))

  #  start_frm = min(frames[0]) 
  #  end_frm = max(frames[-1])
  #  print(start_frm, end_frm)
  #  print(M_info)
  #  connected = connect_rxns(M_info)
  #  print(connected)
  #  for k, v in connected.items():
  #      v.write(os.path.join(cwd,'%s.xyz' %k))
    
   
 
if __name__ == '__main__':
    main()














