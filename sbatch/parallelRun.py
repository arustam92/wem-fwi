import slurm
import subprocess
import time
class runParallelJob:
    """A generalized base class for running a series of jobs in parallel"""
    def __init__(self, jobs):
        """Initialization of the base class
        jobs is the list of myJob objects"""
        self.jobPoll = {}
        for job in jobs:
            self.jobPoll[job.getTag()] = job
        self.jobsRunning = {}
        self.failed = {}

    def startJob(self, tag):
        """How to start a job"""
        #Needs to be overwritten
        raise Exception("Need to override how to startJob")
  
    def getJobsFinished(self):
        """Check to see what jobs are running """ 
        #Force this to be overwritten
        raise Exception("getJobsFinished must be overwritten")
    
    def getJobsFailed(self):
        """Check to see what jobs have failed """ 
        #Force this to be overwritten
        raise Exception("getJobsFailed must be overwritten")
    
    def addToFailed(self,tag):
        #check to see if the job failed before
        if self.failed.count(tag) == 0:
            self.failed[tag] = 0;
        #update the count of failed job
        self.failed[tag] += 1
        #if the job has failed more than twice give up on it
        if self.failed[tag] > 2:
            print ("Giving up on %s"%tag)
        #try to run to job again
        else:
            self.jobPoll[tag] = self.jobsRunning[tag]

    def updateJobs(self): 
        """Check to see if the jobs finished correctly"""
        finished = self.getJobsFinished()
        for tag in finished:
            if not self.jobsRunning[tag].checkJobFinishedCorrectly():
                self.addToFailed(tag)
            del self.jobsRunning[tag]
      
    def allJobsFinished(self):
        """What to do when all the jobs are finished"""

    def runJobs(self,sleepTime=1,maxJobsRunning=1):
        """Run a series of parallel jobs"""
        while len(self.jobPoll) > 0 or len(self.jobsRunning) >0:
            self.updateJobs()
            if len(self.jobsRunning) < maxJobsRunning and len(self.jobPoll) > 0:
                tag, job = self.jobPoll.popitem()
                print("Starting job ", tag)
                self.jobsRunning[tag] = job
                self.startJob(tag)
            time.sleep(sleepTime)
        self.allJobsFinished()

class slurmParallelRun(runParallelJob):
    """Running series of jobs using SLURM scheduler"""
    def __init__(self, jobs, prefix_sbatch,**kw):
        super().__init__(jobs)
        self.slurmPoll = {}
        subprocess.check_output(['mkdir','-p',prefix_sbatch])
        for job in jobs:
            tag = job.getTag()
            script = prefix_sbatch + tag + ".job"
            cmd, stdo, stde = job.returnJobCommand()
            self.slurmPoll[tag] = slurm.slurmCluster(script,**kw,name=tag,stdout=stdo,stderr=stde,commands=cmd)
    
    def startJob(self, tag):
        """How to start a job"""
        self.slurmPoll[tag].submitJob()

    def getJobsFinished(self):
        """Return jobs that are finished """ 
        finished = []
        for tag, job in self.jobsRunning.items():
            status = self.slurmPoll[tag].returnJobStatus()
            if status == "COMPLETED": 
                finished.append(tag) 
        return finished
    
    def getJobsFailed(self):
        """Return jobs that are failed """ 
        failed = []
        for tag, job in self.jobsRunning:
            status = self.slurmPoll[tag].returnJobStatus()
            if status == "FAILED": 
                failed.append(job) 
        return failed
    
