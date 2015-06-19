import os
import numpy as np
from binary import * 
import sys, fiducialMTseq


if __name__ == "__main__":
	SetupScriptPathsFile 	= sys.argv[1]
	alpha 					= float(sys.argv[2])
	beta 					= float(sys.argv[3])
	system 					= sys.argv[4]	
	execfile(SetupScriptPathsFile,globals())
	
	Periods    = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2 ]
	M2s        = [2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0, 10.2, 10.4, 10.6, 10.8, 11.0, 11.2, 11.4, 11.6, 11.8, 12.0, 12.2, 12.4, 12.6, 12.8, 13.0, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2 ]
	Mbhs       = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
 	if beta == 1.0:
		Mbhs = [7.0]
	bhxrb_grid = [(Mbh,P) for Mbh in Mbhs for P in Periods]
	for M2 in M2s:
		SingleTrackFile = SingleGridsDir+'/'+str(M2)+'/LOGS/history.data'
		if os.path.isfile(SingleTrackFile):
			model_number, star_age, star_mass, log10_R = np.loadtxt(SingleTrackFile,skiprows=6,usecols=(0,1,2,37),unpack=True)
			StarRadius = 10.**log10_R * u.Rsun
			
			for Mbh, P in bhxrb_grid:
				fiducial = fiducialMTseq.LagrangianMTseq(alpha, beta, M2, Mbh, P, system)
				RLORadius = roche_lobe(star_mass*u.Msun,Mbh*u.Msun)*period_to_separation(P*u.day,star_mass*u.Msun,Mbh*u.Msun)
				if (RLORadius < StarRadius).any() and (RLORadius[0] > StarRadius[0]):
					if fiducial == 1:
					
						JobDir=MTGridsDir+'/'+str(M2)+'_'+str(Mbh)+'_'+str(P)
						OutFileS=JobDir+'S.data'
						OutFileB=JobDir+'B.data'
						OutFileSgz=JobDir+'S.data.gz'
						OutFileBgz=JobDir+'B.data.gz'
						if (not os.path.isdir(JobDir)) and (not os.path.isfile(OutFileS)) and (not os.path.isfile(OutFileB)) and (not os.path.isfile(OutFileSgz)) and (not os.path.isfile(OutFileBgz)):
							print M2, Mbh, P #, alpha, beta, system

if __name__ == "__main__":
		SetupScriptPathsFile 	= sys.argv[1]
		alpha 					= float(sys.argv[2])
		beta 					= float(sys.argv[3])
		system 					= sys.argv[4]	
		execfile(SetupScriptPathsFile,globals())
		
		if ncpus > 1:
        # Creates jobserver with ncpus workers
        pool = Pool(processes=ncpus)
        for M2, Mbh, P in bhxrb_grid:
            pool.apply_async(MTseqSearchOne, args = (M2, Mbh, P, GridPath,),callback = log_result)
            pool.apply_async(GenerateJoblist, args = (M2, Mbh, P, GridPath,),callback = log_result)
        pool.close()
        pool.join()

        for result in result_list:
            models = models + 1
            if result != "":
                print(result)
                wins = wins + 1