import subprocess

class myJob:
    """A generalized class for describing a parallel job"""
    def __init__(self,tag,singularity_env=False):
        """Initialize a job"""
        self.tag = tag
        self.cmd = []
        if singularity_env:
            self.cmd.append("source %s\n" % singularity_env)
    
    def getTag(self):
        return self.tag

    def preJob(self):
        """How to prepare a job to be run"""
        #By default we don't need to do anything
    
    def checkJobFinishedCorrectly(self):
        """A routine to see if a job finished correctly"""
        return True
        
    def returnJobCommand(self):
        """Return a string containing how to run the job, stdout, stderr"""
        
    def postJob(self):
        """Stuff I need to do after a job has run"""

class fwiJob(myJob):
    """A class for describing a parallel FWI job"""
    def __init__(self, script, parfile, tag, singularity_env=False):
        """Initialize a job"""
        super().__init__(tag,singularity_env)
        self.script = script
        self.parfile = parfile
        self.preJob()
    
    def preJob(self):
        """How to prepare a job to be run"""
        #By default we don't need to do anything
        subprocess.check_output(['mkdir','-p','Inv/'])
        subprocess.check_output(['mkdir','-p','Log/'])
    
    def checkJobFinishedCorrectly(self):
        """A routine to see if a job finished correctly"""
        return True
        
    def returnJobCommand(self):
        """Return a string containing how to run the job, stdout, stderr"""
        stdo = "Log/%s.out" % self.getTag()
        stde = "Log/%s.err" % self.getTag()
        self.cmd.append("python %s parfile=%s\n" % (self.script,self.parfile))
        return self.cmd, stdo, stde

    def postJob(self):
        """Stuff I need to do after a job has run"""

