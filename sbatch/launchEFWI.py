import fwiJob
import argparse
import json
import parallelRun
import numpy as np
def runMultiEpsExtFWI(parser):
    """Initialize a FWI launch with multiple epsilons"""   
                  
    args = parser.parse_args()
    
    parfile = args.parfile
    eps = args.epsilon.split(':')
    epsilon = np.arange(float(eps[0]),float(eps[-1]),float(eps[1])).tolist()
    opts = {}
    options = args.sbatch.split(',')
    for opt in options:
        name = opt.split('=')[0]
        opts[name] = opt.split('=')[-1]
    if args.image: opts["singularity"] = args.image

    with open(parfile,'r') as f:
        dic = json.load(f)
    fwiPoll = []
    ind = parfile.find('.json')
    
    save_as = dic['non_linear_solver']['prefix']    
    for eps in epsilon:
        jobtag = args.prefix + "_%s" % eps
        dic["epsilon"] = eps
        dic['non_linear_solver']['prefix'] = save_as + '_eps_%s' % eps 
        filename = parfile[:ind] + jobtag + parfile[ind:]
        with open(filename,'w') as fout:
            json.dump(dic,fout,indent=4)
        fwiPoll.append(fwiJob.fwiJob(args.fwi,filename,tag=jobtag,singularity_env=args.env))
    parSlurm = parallelRun.slurmParallelRun(fwiPoll,"sbatch/",**opts)
    print("Starting multi-epsilon eFWI with eps = %s"%args.epsilon)
    parSlurm.runJobs(sleepTime=args.sleep,maxJobsRunning=args.maxjobs)

def runMultiSubVpExtFWI(parser):
    """Initialize a FWI launch with multiple subsampling"""   
                  
    args = parser.parse_args()
    
    parfile = args.parfile
    subsample = args.subsample.split(',')
    
    opts = {}
    options = args.sbatch.split(',')
    for opt in options:
        name = opt.split('=')[0]
        opts[name] = opt.split('=')[-1]
    if args.image: opts["singularity"] = args.image

    with open(parfile,'r') as f:
        dic = json.load(f)
    fwiPoll = []
    ind = parfile.find('.json')
    lin_save_as = dic['linear_solver']['prefix']
    save_as = dic['non_linear_solver']['prefix']
    for sub in subsample:
        s = sub.split(':')
        type = dic['pre']['type']
        jobtag = "_vp_%s_%s-%s" % (type,s[0],s[1])
        dic["pre"]['ns'] = [int(s[0]),int(s[1])]
        dic['non_linear_solver']['prefix'] =  save_as + '_%s_%s-%s' % (type,int(s[0]),int(s[1]))
        dic['linear_solver']['prefix'] = lin_save_as + '_%s_%s-%s'% (type,int(s[0]),int(s[1])) 
        filename = parfile[:ind] + jobtag + parfile[ind:]
        with open(filename,'w') as fout:
            json.dump(dic,fout,indent=4)
        fwiPoll.append(fwiJob.fwiJob(args.fwi,filename,tag=jobtag,singularity_env=args.env))
    parSlurm = parallelRun.slurmParallelRun(fwiPoll,"sbatch/",**opts)
    print("Starting multi-subsample FWI with eps = %s"%args.epsilon)
    parSlurm.runJobs(sleepTime=args.sleep,maxJobsRunning=args.maxjobs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fwi", type=str,
                help="fwi script to run")
    parser.add_argument("--parfile", type=str,
                    help="input base parfile for inversion")
    parser.add_argument("--subsample", type=str, default=None,
                    help="ns1:ns2,ns1:ns2,...")     
    parser.add_argument("--epsilon", type=str, default=None,
                    help="eps_min:d_eps:eps_max")  
    parser.add_argument("--sbatch", type=str,
                    help="nodes=...,queue=...,etc...")      
    parser.add_argument("--image", type=str,
                    help="Singularity image")            
    parser.add_argument("--env", type=str,
                    help="File describing env variables for Singularity image")   
    parser.add_argument("--sleep", type=float,
                    help="Max sleep time between jobs")  
    parser.add_argument("--maxjobs", type=float,
                    help="Max num of jobs launched at the same time")  
    parser.add_argument("--prefix", type=str,
                        help="Prefix to assign to a job")  
    
    args = parser.parse_args()
    if args.epsilon: runMultiEpsExtFWI(parser)
    if args.subsample: runMultiSubVpExtFWI(parser)