import sys
import shutil
import os
import time
import subprocess
import numpy as np
from tempfile import mkstemp
from binary import *
import string


def replace(file, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(file)
    for line in old_file:
        new_file.write(line.replace(pattern, subst))
    #close temp file
    new_file.close()
    os.close(fh)
    old_file.close()
    #Remove original file
    os.remove(file)
    #Move new file
    shutil.move(abs_path, file)


if __name__ == '__main__':
    Njobs = int(sys.argv[1]) # number of jobs to submit
    SetupScriptPathsFile = sys.argv[2]
    w_single_track = sys.argv[3]
    execfile(SetupScriptPathsFile,globals())

    num_lines = sum(1 for line in open(JobListFile))
    f = open(JobListFile,"r")
    for i in range(0,min(Njobs, num_lines)):
        lines = f.readlines()
        (M2, Mbh, P) =  lines[i].split()
        #(M2, Mbh, P, alpha, beta, system) =  lines[0].split()

        grid_label = string.split(ExecutableDir,"/executables/")[1]
        NewSubmissionScriptFile = mesa_runs_path+'/job_logs/Job_'+grid_label+'_'+str(M2)+'_'+str(Mbh)+'_'+str(P)+'.slurm'
        shutil.copy(SubmissionScriptFile,NewSubmissionScriptFile)	
        if machine=='baobab':
            os.chdir("../job_logs/")
            replace(NewSubmissionScriptFile, "#SBATCH --job-name=DummyJobName", "#SBATCH --job-name="+grid_label+'_'+str(M2)+'_'+str(Mbh)+'_'+str(P))
            replace(NewSubmissionScriptFile, "#SBATCH --error=DummyJobName", "#SBATCH --error="+grid_label+'_'+str(M2)+'_'+str(Mbh)+'_'+str(P))
            replace(NewSubmissionScriptFile, "#SBATCH --output=DummyJobName", "#SBATCH --output="+grid_label+'_'+str(M2)+'_'+str(Mbh)+'_'+str(P))
            replace(NewSubmissionScriptFile, 'python StartJobRLO.py', 'python '+mesa_runs_path+'/scripts/StartJobRLO.py '+str(M2)+' '+str(Mbh)+' '+str(P)+' '+SetupScriptPathsFile+' '+str(w_single_track))#+' '+str(alpha)+' '+str(beta)+' '+system)
            subprocess.Popen('sbatch Job_'+grid_label+'_'+str(M2)+'_'+str(Mbh)+'_'+str(P)+'.slurm', shell=True).wait()
            os.chdir("../scripts/")
            time.sleep(0.1)
        if machine=='local':
            os.chdir("../job_logs/")
            replace(NewSubmissionScriptFile, "#SBATCH --job-name=DummyJobName", "#SBATCH --job-name="+grid_label+'_'+str(M2)+'_'+str(Mbh)+'_'+str(P))
            replace(NewSubmissionScriptFile, "#SBATCH --error=DummyJobName", "#SBATCH --error="+grid_label+'_'+str(M2)+'_'+str(Mbh)+'_'+str(P))
            replace(NewSubmissionScriptFile, "#SBATCH --output=DummyJobName", "#SBATCH --output="+grid_label+'_'+str(M2)+'_'+str(Mbh)+'_'+str(P))
            replace(NewSubmissionScriptFile, 'python StartJobRLO.py', 'python '+mesa_runs_path+'/scripts/StartJobRLO.py '+str(M2)+' '+str(Mbh)+' '+str(P)+' '+SetupScriptPathsFile+' '+str(w_single_track))#+' '+str(alpha)+' '+str(beta)+' '+system )
            subprocess.Popen('chmod u+x ./Job_'+grid_label+'_'+str(M2)+'_'+str(Mbh)+'_'+str(P)+'.slurm', shell=True).wait()
            subprocess.Popen('nohup ./Job_'+grid_label+'_'+str(M2)+'_'+str(Mbh)+'_'+str(P)+'.slurm &> '+str(M2)+'_'+str(Mbh)+'_'+str(P)+'.out &', shell=True).wait()
            os.chdir("../scripts/")
            time.sleep(0.1)
    f.close()
