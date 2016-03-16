import os, astropy, numpy
from binary import * 
import sys, fiducialMTseq
import math
from multiprocessing import Pool

result_list = []
def log_result(result):
    # This is called whenever MTseqSearchOne returns a result.
    # result_list is modified only by the main process, not the pool workers.
    result_list.append(result)

def Add2JobList(SetupScriptPathsFile, alpha, beta, system, M2, Mbh, P, w_single_track):
	#execfile(SetupScriptPathsFile,globals())
	dum = ""
	if w_single_track > 0:
		#Single track files exists.
		SingleTrackFile = SingleGridsDir+'/'+str(M2)+'/LOGS/history.data'
		if os.path.isfile(SingleTrackFile):
			model_number, star_age, star_mass, log10_R = numpy.loadtxt(SingleTrackFile,skiprows=6,usecols=(0,1,2,37),unpack=True)
			StarRadius = 10.**log10_R * u.Rsun
			fiducial = 1
			if system != "":
				fiducial = fiducialMTseq.LagrangianMTseq(alpha, beta, M2, Mbh, P, system)
			if fiducial != 1:
				return dum
			RLORadius = roche_lobe(star_mass*u.Msun,float(Mbh)*u.Msun)*period_to_separation(float(P)*u.day,star_mass*u.Msun,float(Mbh)*u.Msun)
			if (RLORadius < StarRadius).any() and (RLORadius[0] > StarRadius[0]):
				JobDir=MTGridsDir+'/'+str(M2)+'_'+str(Mbh)+'_'+str(P)
				OutFileS=JobDir+'S.data'
				OutFileB=JobDir+'B.data'
				OutFileSgz=JobDir+'S.data.gz'
				OutFileBgz=JobDir+'B.data.gz'
				if (not os.path.isdir(JobDir)) and (not os.path.isfile(OutFileS)) and (not os.path.isfile(OutFileB)) and (not os.path.isfile(OutFileSgz)) and (not os.path.isfile(OutFileBgz)):
					return M2, Mbh, P #, alpha, beta, system
		else:
			return dum
	elif w_single_track < 0 and system != 'none' and system != 'None':
		fiducial = fiducialMTseq.LagrangianMTseq(alpha, beta, M2, Mbh, P, system)
		if fiducial != 1:
			return dum
		else:
			return M2, Mbh, P
	else:
		return M2, Mbh, P

if __name__ == "__main__":
# Inputs
	ncpus 			= int(sys.argv[1])
	SetupScriptPathsFile 	= sys.argv[2]
	alpha 			= float(sys.argv[3])
	beta 			= float(sys.argv[4])
	system 			= sys.argv[5]
	w_single_track 		= float(sys.argv[6])
	if system == 'none' or systen == 'None':
		system == ""
	
	execfile(SetupScriptPathsFile,globals())
# Initial parameters at onset of RLO.
	Periods    = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, ]
	M2s        = [2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0]
	Mbhs       = [5.0, 5.5]

	bhxrb_grid = [(str(M2),str(Mbh),str(P)) for M2 in M2s for Mbh in Mbhs for P in Periods]	
	
	if ncpus > 1:
		# Creates jobserver with ncpus workers
		pool = Pool(processes=ncpus)
		for M2, Mbh, P in bhxrb_grid:
			pool.apply_async(Add2JobList, args = (SetupScriptPathsFile, alpha, beta, system, M2, Mbh, P,w_single_track,) ,callback = log_result)
		pool.close()
		pool.join()

		models = 0
		wins = 0
		for result in result_list:
			models = models + 1
			if result != "":
				print result[0], result[1], result[2]
				wins = wins + 1
	else:
		for M2, Mbh, P in bhxrb_grid:
			result = Add2JobList(SetupScriptPathsFile, alpha, beta, system, M2, Mbh, P,w_single_track)
			if result != "":
				print result[0], result[1], result[2]
