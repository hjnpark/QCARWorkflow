"""
dsopt.py
"""
from .params import parse_dsopt_args
from .refine import *

def main():
    args_dict = parse_dsopt_args(sys.argv[1:])    
    
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
    maxiter     = args_dict.get('maxiter'   , 100)
    optmethod   = args_dict.get('method'    , 'b3lyp')
    optbasis    = args_dict.get('basis'     , '6-31g(d)')
    coordsys    = args_dict.get('coordsys'  , 'tric')

    client = User(user, password).server()
    
    ds = Dataset(dataset, client).setting('make')
    wf = Workflow(ds=ds, client=client, spec_name=spec_name)
    ds, Mmass =wf.dsoptimize(initial=initial, charge=charge, mult=mult, method=optmethod, basis=optbasis, subsample=subsample, compute = True)
    
    cycle = 0
    while True: 
        ds = Dataset(dataset,client).setting('load') 
        if cycle % 5 == 0 :
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

        wait = (num_opt-comp)*int(Mmass*0.5)
        print("%i/%i calculations are completed" %(comp,num_opt))
        print("Molecular mass = %.2f / Waiting %i seconds" %(Mmass, wait))
        time.sleep(wait) 

        if cycle > 1000:
            raise RuntimeError('Stuck in a while loop checking OptimizationDataset status.')

        if comp == 0 and cycle > 5 :
            raise QCFractalError('Jobs are not recognized by QCFractal server. Try to restart the server and manager')
        cycle += 1


if __name__ == '__main__':
    main()
