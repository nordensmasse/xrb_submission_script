import sys
import shutil
import os
import time
import subprocess
from tempfile import mkstemp
import numpy as np
from binary import *
import fiducialMTseq

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


def StartingModel(M2, Mbh, P):
	SingleTrackFile = SingleGridsDir+'/'+str(M2)+'/LOGS/history.data'
	if not os.path.isfile(SingleTrackFile):
		message = 'No data exist for companion star of '+str(M2)+' solar masses'
		raise NameError(message)

	model_number, star_age, star_mass, log10_R = np.loadtxt(SingleTrackFile,skiprows=6,usecols=(0,1,2,37),unpack=True)
	StarRadius = 10.**log10_R*u.Rsun
	
	RLORadius = roche_lobe(star_mass*u.Msun,Mbh*u.Msun)*period_to_separation(P*u.day,star_mass*u.Msun,Mbh*u.Msun)
	print M2, Mbh, P, StarRadius, RLORadius
	if (RLORadius > StarRadius).all():
		#message = 'Selected star will not fill each Roche Lobe within 13.7Gyr'
		#raise NameError(message)
		StartingModel = -1
	elif RLORadius[0] < StarRadius[0]:
		# message = 'Selected star\`s ZAMS radius is greater than its Roche lobe'
		# raise NameError(message)
		print RLORadius[0], StarRadius[0], log10_R[0]
		StartingModel = -2
	else:
		#We search for the first model at which the star's radius becomes larger than 95% of the Roche Lobe radius
		idx =  np.argmin(np.sign(0.98*RLORadius-StarRadius))
		#We actually pick the model before that as a starting model (checking also if the model before exists)
		StartingModel = model_number[max([idx-1,10])] - model_number[max([idx-1,10])]%10
	return StartingModel 



def SetupOneRLOJob(M2, Mbh, P, w_single_track):
	if w_single_track > 0:
		iModel=int(StartingModel(M2, Mbh, P))
		print iModel
		if iModel < 0:
			if iModel == -1:
				print 'With selected parameter (M2='+str(M2)+', Mbh='+str(Mbh)+', P='+str(P)+'), the star will not fill each Roche Lobe within 13.7Gyr'
			if iModel == -2:
				print 'With selected parameter (M2='+str(M2)+', Mbh='+str(Mbh)+', P='+str(P)+'), the star\`s ZAMS radius is greater than its Roche lobe'
        		if iModel == -3:
            			print 'With selected parameter (M2='+str(M2)+', Mbh='+str(Mbh)+', P='+str(P)+'), alpha ='+str(aplha)+', beta='+str(beta)+' the system '+system+' mass transfer sequence will not fall within the limit of observations.'
            	
			SetupOneRLOJob=0
			return SetupOneRLOJob
		JobDir=MTGridsDir+'/'+str(M2)+'_'+str(Mbh)+'_'+str(P)
		OutFileS=str(M2)+'_'+str(Mbh)+'_'+str(P)+'S.data'
		OutFileB=str(M2)+'_'+str(Mbh)+'_'+str(P)+'B.data'
		ModelFile = 'model'+str(int(100000+iModel))+'.mod'
		shutil.copytree(ExecutableDir, JobDir)	
		shutil.copy(SingleGridsDir+'/'+str(M2)+'/'+ModelFile+'.gz', JobDir+'/')	
		os.chdir(JobDir)
		subprocess.Popen('gunzip *.mod.gz', shell=True).wait()

		replace("inlist_project", "m2", "m2 = "+str(Mbh))
		replace("inlist_project", "m1", "m1 = "+str(M2))
		replace("inlist_project", "initial_period_in_days", "initial_period_in_days = "+str(P))
		replace("inlist1", "saved_model_name", "saved_model_name = \'"+ModelFile+"\'")

		subprocess.Popen('time ./rn', shell=True).wait()
		subprocess.Popen('mv ./LOGS1/history.data ../'+OutFileS, shell=True).wait()
		subprocess.Popen('mv ./binary_history.data ../'+OutFileB, shell=True).wait()
		os.chdir("../")
		subprocess.Popen('rm -rf '+str(M2)+'_'+str(Mbh)+'_'+str(P), shell=True).wait()
		subprocess.Popen('gzip -f '+OutFileS, shell=True).wait()
		subprocess.Popen('gzip -f '+OutFileB, shell=True).wait()
		os.chdir(mesa_runs_path+'/scripts')


		SetupOneRLOJob=1
		return SetupOneRLOJob
	elif w_single_track < 0:	
		JobDir=MTGridsDir+'/'+str(M2)+'_'+str(Mbh)+'_'+str(P)
		OutFileS=str(M2)+'_'+str(Mbh)+'_'+str(P)+'S.data'
		OutFileB=str(M2)+'_'+str(Mbh)+'_'+str(P)+'B.data'
		shutil.copytree(ExecutableDir, JobDir)	
		os.chdir(JobDir)
		subprocess.Popen('gunzip *.mod.gz', shell=True).wait()

		replace("inlist_project", "m2", "m2 = "+str(Mbh))
		replace("inlist_project", "m1", "m1 = "+str(M2))
		replace("inlist_project", "initial_period_in_days", "initial_period_in_days = "+str(P))

		subprocess.Popen('time ./rn', shell=True).wait()
		subprocess.Popen('mv ./LOGS1/history.data ../'+OutFileS, shell=True).wait()
		subprocess.Popen('mv ./binary_history.data ../'+OutFileB, shell=True).wait()
		os.chdir("../")
		subprocess.Popen('rm -rf '+str(M2)+'_'+str(Mbh)+'_'+str(P), shell=True).wait()
		subprocess.Popen('gzip -f '+OutFileS, shell=True).wait()
		subprocess.Popen('gzip -f '+OutFileB, shell=True).wait()
		os.chdir(mesa_runs_path+'/scripts')

		SetupOneRLOJob=1
		return SetupOneRLOJob

if __name__ == '__main__':
	M2                      = float(sys.argv[1])
	Mbh                     = float(sys.argv[2])
	P                       = float(sys.argv[3])
	SetupScriptPathsFile    = sys.argv[4]
	w_single_track 		= float(sys.argv[5])
	
	execfile(SetupScriptPathsFile, globals())

	SetupOneRLOJob(M2, Mbh, P, w_single_track)
