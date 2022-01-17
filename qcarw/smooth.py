"""
smooth.py
"""
from .params import parse_smooth_args
from .molecule import Molecule
from .refine import *


def main():
    args_dict = parse_smooth_args(sys.argv[1:])    

    input_name = args_dict['input'].split('.')
    temp_name = '.'.join(input_name[:-1]) 
    spec_name = temp_name + '_spec'


    if 'dataset' in args_dict:
        pass
    else:
        args_dict['dataset'] = temp_name + '_ds'
    
    filename    = args_dict.get('input')
    dataset     = args_dict.get('dataset') 
    user        = args_dict.get('user'      , None)
    password    = args_dict.get('passowrd'  , None)

    client = User(user, password).server()
    
    ds = Dataset(dataset, client).setting('load')
    wf = Workflow(client=client, ds=ds, spec_name=spec_name)
    smoothed = wf.smoothing(filename)    

    print("Smoothing procedure is done. %i number of reactions were detected." %len(smoothed))

if __name__ == '__main__':
    main()
