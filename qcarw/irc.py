"""
irc.py
"""
from .params import parse_irc_args
from .refine import *

def main():
    args_dict = parse_irc_args(sys.argv[1:])    
    
    initial     = args_dict.get('input')
    user        = args_dict.get('user'      , None)
    password    = args_dict.get('passowrd'  , None)
    charge      = args_dict.get('mult'      , 0)
    mult        = args_dict.get('charge'    , 1)
    method      = args_dict.get('method'    , 'b3lyp')
    basis       = args_dict.get('basis'     , '6-31+g(d)')
    coordsys    = args_dict.get('coordsys'  , 'cart')
    trust       = args_dict.get('trust'     , 0.2)


    client = User(user, password).server()
    
    wf = Workflow(client=client, ds=None, spec_name=None)
    process=wf.irc(initial=initial, charge=charge, mult=mult, method=method, basis=basis, coordsys=coordsys, trust=trust)
    print("IRC step is completed.")

if __name__ == '__main__':
    main()
