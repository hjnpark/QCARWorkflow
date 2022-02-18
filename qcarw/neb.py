"""
neb.py
"""
from .params import parse_neb_args
from .refine import *

def main():
    args_dict = parse_neb_args(sys.argv[1:])    
    
    initial     = args_dict.get('input')
    user        = args_dict.get('user'      , None)
    password    = args_dict.get('passowrd'  , None)
    method      = args_dict.get('method'    , 'b3lyp')
    basis       = args_dict.get('basis'     , '6-31+g(d,p)')
    images      = args_dict.get('images'    , '21')
    coordsys    = args_dict.get('coordsys'  , 'cart')
    nebk        = args_dict.get('nebk'      , 1)
    avgg        = args_dict.get('avgg'      , 0.025)
    maxg        = args_dict.get('maxg'      , 0.05)
    charge      = args_dict.get('charge'    , 0)
    mult        = args_dict.get('mult'      , 1)
    ew          = args_dict.get('ew'        , False)


    client = User(user, password).server()
    
    wf = Workflow(client=client, ds=None, spec_name=None)
    wf.neb(initial=initial, charge=charge, mult=mult, method=method, basis=basis, images=images, coordsys=coordsys, ew=ew, nebk=nebk, avgg=avgg, maxg=maxg)

    print("NEB calculation is completed.")

if __name__ == '__main__':
    main()
