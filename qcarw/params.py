"""
params.py
"""
import argparse

class ArgumentParserWithFile(argparse.ArgumentParser):
    """
    This class is from geomeTRIC/geometric/parmas.py
    ------------------------------------------------
    It will read arguments between the following two lines from a text file.
    $options
    $end
    """
    def __init__(self, *args, **kwargs):
        self.in_options = False
        super(ArgumentParserWithFile, self).__init__(*args, **kwargs)
    
    def convert_arg_line_to_args(self,line):
        line = line.split("#")[0].strip() # Don't lower-case because may be case sensitive
        if '$options' in line:
            self.in_options = True
        elif '$' in line:
            self.in_options = False
        elif self.in_options:
            s = line.split()
            s[0] = '--'+s[0]
            return s
        return []

def str2bool(v):
    """Allows command line options such as 'Yes' and 'True' to be converted into Booleans."""
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'yeah', 'ok', 'on', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'nope', 'off', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def parse_refine_args(*args):

    """ 
    Read user input from the command line interface. 
    Designed to be called by qcarw.main() passing in sys.argv[1:] 
    
    Avoid setting default values for variables here. The default values of certain variables 
    depends on the values of other variables. 
    """

    parser = ArgumentParserWithFile(add_help=False, formatter_class=argparse.RawTextHelpFormatter, fromfile_prefix_chars='@')

    grp_qcf = parser.add_argument_group('qcfserver', 'general information for qcfractal server')
    grp_qcf.add_argument('--user', type=str, help='User ID for QCFractal server\n')
    grp_qcf.add_argument('--password', type=str, help='QCFractal password\n')
    grp_qcf.add_argument('--dataset', type=str, help='OptimizationDataset name, default = str(input filename + _ds) \n')
    
    grp_optparam = parser.add_argument_group('optparam', 'Optimization parameters for the initial MD trajectory')
    grp_optparam.add_argument('input', type=str, help='REQUIRED, xyz file (MD trajectory) name that needs to be refined\n')
    grp_optparam.add_argument('--charge', type=int, help='Molecular charge, default = 0\n')
    grp_optparam.add_argument('--mult', type=int, help='Spin multiplicity, default = 1\n')
    grp_optparam.add_argument('--subsample', type=int, help='Frame interval for subsampling trajectories, default = 10\n')
    grp_optparam.add_argument('--maxiter', type=int, help='Maximum number of optimization steps, default = 100\n')
    grp_optparam.add_argument('--optcrdsys',  type=str, help='Coordinate system for geometry optimization:\n' 
                          '"cart" = Cartesian coordinate system\n'
                          '"tric" for Translation-Rotation Internal Coordinates (default)\n'
                          '"prim" = Primitive (a.k.a redundant internal coordinates)\n'
                          '"dlc" = Delocalized Internal Coordinates\n'
                          '"hdlc" = Hybrid Delocalized Internal Coordinates\n'
                          '"tric-p" for primitive Translation-Rotation Internal Coordinates (no delocalization)\n')
    grp_optparam.add_argument('--optmethod', type=str, help='Electronic structure method, default = b3lyp\n') 
    grp_optparam.add_argument('--optbasis', type=str, help='Basis set, default 6-31g(d)')
    grp_optparam.add_argument('--tsmethod', type=str, help='Electronic structure method for transition state structure optimization, default = b3lyp\n')
    grp_optparam.add_argument('--tsbasis', type=str, help='Basis set for transition state structure optimization, default = 6-31+g(d,p)')   
 
    grp_nebparam = parser.add_argument_group('nebparam', 'The Nudged Elastic Band method parameters')
    grp_nebparam.add_argument('--images', type=int, help='Number of images im the NEB chain, default = 21\n')
    grp_nebparam.add_argument('--nebmethod', type=str, help='Electronic structure method, default = b3lyp\n') 
    grp_nebparam.add_argument('--nebbasis', type=str, help='Basis set, default = 6-31+g(d,p)\n')
    grp_nebparam.add_argument('--nebcrdsys', type=str, help='Coordinate system for NEB:\n' 
                          '"cart" = Cartesian coordinate system (default)\n'
                          '"tric" for Translation-Rotation Internal Coordinates\n'
                          '"prim" = Primitive (a.k.a redundant internal coordinates)\n '
                          '"dlc" = Delocalized Internal Coordinates\n'
                          '"hdlc" = Hybrid Delocalized Internal Coordinates\n'
                          '"tric-p" for primitive Translation-Rotation Internal Coordinates (no delocalization)\n')
    grp_nebparam.add_argument('--nebk', type=float, help='Spring constant in kcal/mol/Ang^2, defualt = 1\n')
    grp_nebparam.add_argument('--avgg', type=float, help='Average RMS-gradient for the convergence in ev/Ang, default = 0.025\n')
    grp_nebparam.add_argument('--maxg', type=float, help='Maximum RMS-gradient for the convergence in ev/Ang, default = 0.05\n')
    grp_nebparam.add_argument('--ew', type=str2bool, help='Provide "yes" to perform the energy weighted NEB\n')

    grp_ircparam = parser.add_argument_group('ircparam', 'The Intrinsic Reaction Coordinate method parameters')
    #grp_ircparam.add_argument('--irccrdsys', type=str, help='Coordinate sustem for IRC. Only Cartesian is supported for now\n')
    grp_ircparam.add_argument('--ircstep', type=float, help='A step size that the IRC method will take, default = 0.2\n')

    grp_help = parser.add_argument_group('help', 'Get help')
    grp_help.add_argument('-h', '--help', action='help', help='Show this help message and exit')

    # Keep all arguments whose values are not None, so that the setting of default values
    args_dict = {}
    for k, v in vars(parser.parse_args(*args)).items():
        if v is not None:
            args_dict[k] = v

    return args_dict


def parse_dsopt_args(*args):

    """ 
    Read user input from the command line interface. 
    Designed to be called by dsopt.main() passing in sys.argv[1:] 
    
    Avoid setting default values for variables here. The default values of certain variables 
    depends on the values of other variables. 
    """

    parser = ArgumentParserWithFile(add_help=False, formatter_class=argparse.RawTextHelpFormatter, fromfile_prefix_chars='@')

    grp_qcf = parser.add_argument_group('qcfserver', 'general information for qcfractal server')
    grp_qcf.add_argument('--user', type=str, help='User ID for QCFractal server\n')
    grp_qcf.add_argument('--password', type=str, help='QCFractal password\n')
    grp_qcf.add_argument('--dataset', type=str, help='OptimizationDataset name, default = str(input filename + _ds) \n')
    
    grp_optparam = parser.add_argument_group('optparam', 'Optimization parameters for the initial MD trajectory')
    grp_optparam.add_argument('input', type=str, help='REQUIRED, xyz file (MD trajectory) name that needs to be refined\n')
    grp_optparam.add_argument('--charge', type=int, help='Molecular charge, default = 0\n')
    grp_optparam.add_argument('--mult', type=int, help='Spin multiplicity, default = 1\n')
    grp_optparam.add_argument('--subsample', type=int, help='Frame interval for subsampling trajectories, default = 10\n')
    grp_optparam.add_argument('--maxiter', type=int, help='Maximum number of optimization steps, default = 500\n')
    grp_optparam.add_argument('--coordsys',  type=str, help='Coordinate system for geometry optimization:\n' 
                          '"cart" = Cartesian coordinate system\n'
                          '"tric" for Translation-Rotation Internal Coordinates (default)\n'
                          '"prim" = Primitive (a.k.a redundant internal coordinates)\n'
                          '"dlc" = Delocalized Internal Coordinates\n'
                          '"hdlc" = Hybrid Delocalized Internal Coordinates\n'
                          '"tric-p" for primitive Translation-Rotation Internal Coordinates (no delocalization)\n')
    grp_optparam.add_argument('--method', type=str, help='Electronic structure method, default = b3lyp\n') 
    grp_optparam.add_argument('--basis', type=str, help='Basis set, default 6-31g(d)')

    grp_help = parser.add_argument_group('help', 'Get help')
    grp_help.add_argument('-h', '--help', action='help', help='Show this help message and exit')

    # Keep all arguments whose values are not None, so that the setting of default values
    args_dict = {}
    for k, v in vars(parser.parse_args(*args)).items():
        if v is not None:
            args_dict[k] = v

    return args_dict

def parse_smooth_args(*args):

    """ 
    Read user input from the command line interface. 
    Designed to be called by smooth.main() passing in sys.argv[1:] 
    
    Avoid setting default values for variables here. The default values of certain variables 
    depends on the values of other variables. 
    """

    parser = ArgumentParserWithFile(add_help=False, formatter_class=argparse.RawTextHelpFormatter, fromfile_prefix_chars='@')
        
    
    grp_qcf = parser.add_argument_group('qcfserver', 'general information for qcfractal server')
    grp_qcf.add_argument('input', type=str, help='REQUIRED, xyz file (MD trajectory) that was optimized\n')
    grp_qcf.add_argument('--user', type=str, help='User ID for QCFractal server\n')
    grp_qcf.add_argument('--password', type=str, help='QCFractal password\n')
    grp_qcf.add_argument('--dataset', type=str, help='OptimizationDataset name that you want to get smoothed pathways \n')
    

    grp_help = parser.add_argument_group('help', 'Get help')
    grp_help.add_argument('-h', '--help', action='help', help='Show this help message and exit')

    # Keep all arguments whose values are not None, so that the setting of default values
    args_dict = {}
    for k, v in vars(parser.parse_args(*args)).items():
        if v is not None:
            args_dict[k] = v

    return args_dict

def parse_neb_args(*args):

    """ 
    Read user input from the command line interface. 
    Designed to be called by neb.main() passing in sys.argv[1:] 
    
    Avoid setting default values for variables here. The default values of certain variables 
    depends on the values of other variables. 
    """

    parser = ArgumentParserWithFile(add_help=False, formatter_class=argparse.RawTextHelpFormatter, fromfile_prefix_chars='@')

    grp_qcf = parser.add_argument_group('qcfserver', 'general information for qcfractal server, they are optional for neb.')
    grp_qcf.add_argument('--user', type=str, help='User ID for QCFractal server\n')
    grp_qcf.add_argument('--password', type=str, help='QCFractal password\n')
 
    grp_nebparam = parser.add_argument_group('nebparam', 'The Nudged Elastic Band method parameters')
    grp_nebparam.add_argument('input', type=str, help='REQURIED, NEB ready xyz file name that can be used as an input')
    grp_nebparam.add_argument('--charge', type=int, help='Molecular charge, default = 0\n')
    grp_nebparam.add_argument('--mult', type=int, help='Spin multiplicity, default = 1\n')
    grp_nebparam.add_argument('--images', type=int, help='Number of images im the NEB chain, default = 21\n')
    grp_nebparam.add_argument('--method', type=str, help='Electronic structure method, default = b3lyp\n') 
    grp_nebparam.add_argument('--basis', type=str, help='Basis set, default = 6-31+g(d, p)\n')
    grp_nebparam.add_argument('--coordsys', type=str, help='Coordinate system for neb:\n' 
                          '"cart" = Cartesian coordinate system (default)\n'
                          '"tric" for Translation-Rotation Internal Coordinates\n'
                          '"prim" = Primitive (a.k.a redundant internal coordinates)\n'
                          '"dlc" = Delocalized Internal Coordinates\n'
                          '"hdlc" = Hybrid Delocalized Internal Coordinates\n'
                          '"tric-p" for primitive Translation-Rotation Internal Coordinates (no delocalization)\n')
    grp_nebparam.add_argument('--nebk', type=float, help='Spring constant in kcal/mol/Ang^2, defualt = 1\n')
    grp_nebparam.add_argument('--avgg', type=float, help='Average RMS-gradient for the convergence in ev/Ang, default = 0.025\n')
    grp_nebparam.add_argument('--maxg', type=float, help='Maximum RMS-gradient for the convergence in ev/Ang, default = 0.05\n')
    grp_nebparam.add_argument('--ew', type=str2bool, help='provide "yes" to perform the energy weighted NEB\n')

    grp_help = parser.add_argument_group('help', 'Get help')
    grp_help.add_argument('-h', '--help', action='help', help='Show this help message and exit')

    # Keep all arguments whose values are not None, so that the setting of default values
    args_dict = {}
    for k, v in vars(parser.parse_args(*args)).items():
        if v is not None:
            args_dict[k] = v

    return args_dict

def parse_opt_args(*args):

    """ 
    Read user input from the command line interface. 
    Designed to be called by opt.main() passing in sys.argv[1:] 
    
    Avoid setting default values for variables here. The default values of certain variables 
    depends on the values of other variables. 
    """

    parser = ArgumentParserWithFile(add_help=False, formatter_class=argparse.RawTextHelpFormatter, fromfile_prefix_chars='@')

    grp_qcf = parser.add_argument_group('qcfserver', 'general information for qcfractal server, they are optional for single optimizations.')
    grp_qcf.add_argument('--user', type=str, help='User ID for QCFractal server\n')
    grp_qcf.add_argument('--password', type=str, help='QCFractal password\n')

    grp_optparam = parser.add_argument_group('optparam', 'Optimization parameters for the initial MD trajectory')
    grp_optparam.add_argument('input', type=str, help='REQUIRED, xyz file (MD trajectory) name that needs to be refined\n')
    grp_optparam.add_argument('--charge', type=int, help='Molecular charge, default = 0\n')
    grp_optparam.add_argument('--mult', type=int, help='Spin multiplicity, default = 1\n')
    grp_optparam.add_argument('--subsample', type=int, help='Frame interval for subsampling trajectories, default = 10\n')
    grp_optparam.add_argument('--maxiter', type=int, help='Maximum number of optimization steps, default = 500\n')
    grp_optparam.add_argument('--coordsys',  type=str, help='Coordinate system for geometry optimization:\n' 
                          '"cart" = Cartesian coordinate system\n'
                          '"tric" for Translation-Rotation Internal Coordinates (default)\n'
                          '"prim" = Primitive (a.k.a redundant internal coordinates)\n'
                          '"dlc" = Delocalized Internal Coordinates\n'
                          '"hdlc" = Hybrid Delocalized Internal Coordinates\n'
                          '"tric-p" for primitive Translation-Rotation Internal Coordinates (no delocalization)\n')
    grp_optparam.add_argument('--method', type=str, help='Electronic structure method, default = b3lyp\n') 
    grp_optparam.add_argument('--basis', type=str, help='Basis set, default 6-31g(d)\n')
    grp_optparam.add_argument('--transition', type=str2bool, help='provide "yes" to perform a TS optimization\n')
 

    grp_help = parser.add_argument_group('help', 'Get help')
    grp_help.add_argument('-h', '--help', action='help', help='Show this help message and exit')

    # Keep all arguments whose values are not None, so that the setting of default values
    args_dict = {}
    for k, v in vars(parser.parse_args(*args)).items():
        if v is not None:
            args_dict[k] = v

    return args_dict

def parse_irc_args(*args):

    """ 
    Read user input from the command line interface. 
    Designed to be called by neb.main() passing in sys.argv[1:] 
    
    Avoid setting default values for variables here. The default values of certain variables 
    depends on the values of other variables. 
    """

    parser = ArgumentParserWithFile(add_help=False, formatter_class=argparse.RawTextHelpFormatter, fromfile_prefix_chars='@')

    grp_qcf = parser.add_argument_group('qcfserver', 'general information for qcfractal server, they are optional for neb.')
    grp_qcf.add_argument('--user', type=str, help='User ID for QCFractal server\n')
    grp_qcf.add_argument('--password', type=str, help='QCFractal password\n')
 
    grp_ircparam = parser.add_argument_group('ircparam', 'The Intrinsic Reaction Coordinate method parameters')
    grp_ircparam.add_argument('input', type=str, help='REQURIED, IRC ready xyz file (optimized TS structure) name that can be used as an input')
    grp_ircparam.add_argument('--charge', type=int, help='Molecular charge, default = 0\n')
    grp_ircparam.add_argument('--mult', type=int, help='Spin multiplicity, default = 1\n')
    grp_ircparam.add_argument('--method', type=str, help='Electronic structure method, default = b3lyp\n') 
    grp_ircparam.add_argument('--basis', type=str, help='Basis set, default = 6-31+g(d,p)\n')
   # grp_nebparam.add_argument('--coordsys', type=str, help='Coordinate system for irc:\n' 
   #                       '"cart" = Cartesian coordinate system (default)\n'
   #                       '"tric" for Translation-Rotation Internal Coordinates\n'
   #                       '"prim" = Primitive (a.k.a redundant internal coordinates)\n'
   #                       '"dlc" = Delocalized Internal Coordinates\n'
   #                       '"hdlc" = Hybrid Delocalized Internal Coordinates\n'
   #                       '"tric-p" for primitive Translation-Rotation Internal Coordinates (no delocalization)\n')
    grp_ircparam.add_argument('--trust', type=float, help='Trust step size, default = 0.2\n')

    grp_help = parser.add_argument_group('help', 'Get help')
    grp_help.add_argument('-h', '--help', action='help', help='Show this help message and exit')

    # Keep all arguments whose values are not None, so that the setting of default values
    args_dict = {}
    for k, v in vars(parser.parse_args(*args)).items():
        if v is not None:
            args_dict[k] = v

    return args_dict














