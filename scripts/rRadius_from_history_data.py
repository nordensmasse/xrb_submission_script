# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 10:34:20 2015

READ THE INITIAL ZAMS RADIUS OF A SINGLE STAR.

@author: MS
"""

import os, numpy
M2s = numpy.linspace(2.0, 15, (15-2)/0.1+1)

SingleGridsDir = '/home/sorensen/mesa_runs/SINGLEgrids/Z0.006/'

for M2 in M2s:
    SingleTrackFile = SingleGridsDir+'/'+str(M2)+'/LOGS/history.data'
    if not os.path.isfile(SingleTrackFile):
        message = 'Not data exist for companion star of '+str(M2)+' solar masses'
        raise NameError(message)

    model_number, star_age, star_mass, log10_R = numpy.loadtxt(SingleTrackFile,skiprows=6,usecols=(0,1,2,37),unpack=True)
    StarRadius = 10.**log10_R[0]
    print str(M2)+'    '+str(StarRadius)
        
