import subprocess
import os

class slurmCluster:
    """Class for creating job scripts for SLURM"""
    def __init__(self, fileout, **kw):
        """Initialize the job creator class and default most parameters"""
        self.singularityCmd = ""
        if "singularity" in kw:
            self.singularityCmd = "singularity exec %s" % kw["singularity"]
        self.initializeDefaultParams()
        self.createJobScript(fileout,**kw)
      
    def createJobScript(self, fileout, **kw):
        """We can override any parameters then write out file"""
        self.overrideDefaultParams(kw)
        self.writeFile(fileout)
        self.script = fileout
    
    def initializeDefaultParams(self):
        """Setup some defaults"""
        self.params = {}
        self.params["name"] = os.environ["USER"] + "_" + "awesomeness" #name of the job in the queue
        self.params["nodes"] = 1
        self.params["cores"] = 1
        self.params["queue"] = "twohour"
        self.params["stdout"] = "/data/cees/%s/stdout"%os.environ["USER"]
        self.params["stderr"] = "/data/cees/%s/stderr"%os.environ["USER"]

    def overrideDefaultParams(self, lst):
        """Override the defaults"""
        for k,v in lst.items():
            self.params[k] = v
          
    def writeFile(self, fileout):
        """Write job file"""
        if "commands" not in self.params:
            raise Exception("Must specify commands in either initialization or createJobScript")
        cmdfile = fileout + '.sh'
        try: 
            f = open(fileout, "w")
            cmd = open(cmdfile,'w')
        except:
            raise Exception("Could not open "%fileout)
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --job-name=%s\n"%self.params["name"])
        f.write("#SBATCH --partition=%s\n"%self.params["queue"])
        f.write("#SBATCH --nodes=%d\n"%int(self.params["nodes"]))
        f.write("#SBATCH --cpus-per-task=%d\n"%int(self.params["cores"]))
        f.write("#SBATCH --error=%s\n"%self.params["stderr"])
        f.write("#SBATCH --output=%s\n"%self.params["stdout"])
        if "ntasks" in self.params:
            f.write("#SBATCH --ntasks=%d\n"%int(self.params["ntasks"]))
        if "mail" in self.params:
            f.write("#SBATCH --mail-user=%s\n"%self.params["mail"])
        if "time" in self.params:
            f.write("#SBATCH --time=%s\n"%self.params["time"])
        f.write("cd $SLURM_SUBMIT_DIR\n")
        if type(self.params["commands"]) is list:
            for c in self.params["commands"]:
                cmd.write("%s\n"%c)
        elif type(self.params["commands"]) is str:
            cmd.write("%s\n"%self.params["commands"])
        else:
            raise Exception("Commands must be a list or string")
        f.write("%s /bin/bash %s\n" % (self.singularityCmd, cmdfile))
        f.close()
        cmd.close()

    def submitJob(self):
        """Submit job and return job id number"""
        out = subprocess.check_output(["sbatch", self.script])
        out = out.decode("utf-8").split("\n")[0].split(" ")
        self.jobId = out[-1]

    def returnJobStatus(self):
        """Return dictionaryjob id-> [error, running, funished, queued]"""
        status = ""
        for state in ["COMPLETED","RUNNING","FAILED"]:
            lines=subprocess.check_output([
                "sacct","-j","%s" % self.jobId,"-o","State"
                ]).decode("utf-8").split("\n")
            lines.pop(0)
            lines.pop(0)
            if lines[0] == state: 
                status = state
                break
        return status