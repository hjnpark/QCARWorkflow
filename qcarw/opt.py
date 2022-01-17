"""
opt.py
"""
from .params import parse_opt_args
from .refine import *


def main():
    args_dict = parse_opt_args(sys.argv[1:])    
    
    initial     = args_dict.get('input')
    user        = args_dict.get('user'      , None)
    password    = args_dict.get('passowrd'  , None)
    method      = args_dict.get('method'    , 'b3lyp')
    basis       = args_dict.get('basis'     , '6-31+g(d)')
    coordsys    = args_dict.get('coordsys'  , 'tric')
    maxiter     = args_dict.get('maxiter'   , 500)
    charge      = args_dict.get('charge'    , 0)
    mult        = args_dict.get('mult'      , 1)
    ts          = args_dict.get('transition' , False)


    client = User(user, password).server()
    
    wf = Workflow(client=client, ds=None, spec_name=None)
    qcel, M, energy=wf.optimize(initial=initial, charge=charge, mult=mult, method=method, basis=basis, coordsys=coordsys, maxiter=maxiter, ts=ts)
    if ts:
        M.write('ts-optimized.xyz')
        print("Transition state structure was optimized. ts-optimized.xyz was created.")
    else:
        M.write('optimized.xyz')
        print("Optimization is completed. optimized.xyz was created.")

if __name__ == '__main__':
    main()
