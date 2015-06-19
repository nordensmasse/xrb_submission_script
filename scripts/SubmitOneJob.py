import sys
import shutil
import os
import time
import subprocess
import fileinput
import numpy as np
import time
from tempfile import mkstemp
from binary import *


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
    M2 = float(sys.argv[1])
    Mbh = float(sys.argv[2])
    P = float(sys.argv[3])

    shutil.copy('example.slurm','Job_'+str(M2)+'_'+str(Mbh)+'_'+str(P)+'.slurm')	
    replace('Job_'+str(M2)+'_'+str(Mbh)+'_'+str(P)+'.slurm', "#SBATCH --job-name=DummyJobName", "#SBATCH --job-name="+str(M2)+'_'+str(Mbh)+'_'+str(P))
    replace('Job_'+str(M2)+'_'+str(Mbh)+'_'+str(P)+'.slurm', "#SBATCH --error=DummyJobName", "#SBATCH --error="+str(M2)+'_'+str(Mbh)+'_'+str(P))
    replace('Job_'+str(M2)+'_'+str(Mbh)+'_'+str(P)+'.slurm', "#SBATCH --output=DummyJobName", "#SBATCH --output="+str(M2)+'_'+str(Mbh)+'_'+str(P))
    replace('Job_'+str(M2)+'_'+str(Mbh)+'_'+str(P)+'.slurm', "python StartJobRLO.py", "python StartJobRLO.py "+str(M2)+' '+str(Mbh)+' '+str(P))
#    subprocess.Popen('sbatch Job_'+str(M2)+'_'+str(Mbh)+'_'+str(P)+'.slurm', shell=True).wait()
    subprocess.Popen('chmod u+x ./Job_'+str(M2)+'_'+str(Mbh)+'_'+str(P)+'.slurm', shell=True).wait()
    subprocess.Popen('nohup ./Job_'+str(M2)+'_'+str(Mbh)+'_'+str(P)+'.slurm &> '+str(M2)+'_'+str(Mbh)+'_'+str(P)+'.out &', shell=True).wait()

                     
