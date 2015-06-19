# -*- coding: utf-8 -*-
"""
Created on Tue Feb  23 18:20:15 2015

@author: Mads Sørensen
"""
import numpy, sys, os, odespy
from numpy import where as where
#import MTseq_npy

def constant_G():
#Function transforms gravitational constant G from SI to solar; radii, masses, and Myr.
    Myr    = 365*24*3600*1e6
    Rsun  = 6.9551*1e8
    Msun  = 1.999*1e30
    G     = 6.67*1e-11*Myr**2*Msun/Rsun**3 # m**3 s**-2 kg**-1
#    print 'Loaded Gravitational constant G in'
#    print 'units [Myr**-2 Msun**-1 Rsun**3]'
    return G
G = constant_G()

#Lagrangian Formulation of Single Mass Transfer sequence
def orbit_separation(M1, alpha, beta, a, M2low,M2up,dM):
    # Integrate forward in time, hence mass flows from donor M2 onto accretor.
    # Alpha is fraction of mass lost from the donor and transferred towards the accretor.
    # Beta is the fraction of mass captured by the accretor.
    
    import odespy

    #Given input determine orbital separation    
    def int_Mdonor(u,M2):
        # Funciton to integrate over M2, the donor mass
        a, M1  = u
        q      = M2/M1
        dM1dM2 = -alpha*beta
        dadM2  = -2*a/M2*(alpha*(1-q)+q/2/(1+q)*(1-alpha*beta))
        return [dadM2, dM1dM2]
        
    # Initial conditions
    u_init = [a,M1]
    u      = numpy.array([2])    
    solver = odespy.RK4(int_Mdonor)    
    solver.set_initial_condition(u_init)
    N      = int(round(abs(M2up-M2low)/dM))    # Number of points 
    # Integrate
    mass_points = numpy.linspace(M2up,M2low,N+1)
    u, M2       = solver.solve(mass_points)
    result = numpy.vstack((M2,u[:,0],u[:,1]))
    return result


#The black hole spin
def BH_spin(Mi,M0):
    # From Thorne, K. S. (1974). Disk-accretion onto a black hole. ii. evolution of the hole.
    # The Astrophysical Journal, 191:507–520.
    M = Mi
    if (M/Mi <= numpy.sqrt(6.)): 
        #eq. 2b solved for M
        M = 3. * Mi * ( numpy.sin( M0/3./Mi + numpy.arcsin(1./3.) ) )
        #eq. 2a
        a = (2./3.)**0.5*Mi/M * (4.-(18.*Mi**2/M**2 - 2)**0.5)

    if M/Mi >= numpy.sqrt(6.):
        #eq. 2b solved for M
        M = (M0 - 3.*Mi* (numpy.sqrt(numpy.arcsin(2./3.)) - numpy.arcsin(1./3.)) ) / numpy.sqrt(3.) + numpy.sqrt(6.)*Mi
        #eq. 2a
        a = 1.
    return a

class ObservedSystem(object):
	def __init__(self, SystemName):
		self.SystemName = SystemName
		obs_data_path = '/Users/MS/mesa_runs/submission_script/systems/'
		self.datafile = obs_data_path+SystemName+'.dat'
		self.exists = os.path.exists(self.datafile)
		if self.exists:
			# All values are in solar units, except Teff (in K), P (in days), and logg (in cgs)
			self.Value={} 
			self.Error={} 
			d= numpy.loadtxt(self.datafile, dtype=[("Properties","S5"),('Values',float),('Errors',float)])
			self.Properties=d['Properties']
			for Property, Value, Error in zip(d['Properties'], d['Values'], d['Errors']):
				self.Value[Property]=Value
				self.Error[Property]=Error
		else:
			print self.datafile
			sys.exit("Obs. data for system ",SystemName," do not exist")

def LagrangianMTseq(alpha, beta, M2, Mbh, P, system):
	M2 = float(M2)
	Mbh = float(Mbh)
	P = float(P)
	
	Nsigma = 2.0 		# Number of accepted standard deviations away observed value.
	
	#Read in system values to compare with
	b = ObservedSystem(system)
	Pobs = b.Value['P']
	#if b.Value['M2'] > 4*b.Error['M2']:
	#	M2low 	= b.Value['M2']-(Nsigma+2.)*b.Error['M2']
	#else:
	M2low 	= 1.0

	# Array with accepted models
	F_M2       = 0
	F_M1       = 0
	F_a        = 0
	F_spin 	   = 0

	p          	= P/1e6/365.0 									# To Myr
	dM2        	= 0.01 										# Integration step size
	a_h  		= ((p/(2.*numpy.pi))**2*G*(Mbh+M2))**(1./3.)     			# $R_{\odot}$

	result     	= orbit_separation(Mbh , 1.0-alpha, 1.0-beta, a_h, M2low, M2, dM2) 	# Bulk integrate

	M2_cal     	= result[0,:] 									# $M_{\odot}$
	a_cal      	= result[1,:] 									# $R_{\odot}$
	M1_cal     	= result[2,:] 									# $M_{\odot}$
		
	# determine orbital period            
	omega_cal  	= numpy.sqrt(G*(M1_cal+M2_cal)/a_cal**3)
	P_cal      	= 2.*numpy.pi/omega_cal*1e6*365.0

	# Criteria for interpolation
	# Avoid extrapolation
	if min(P_cal) > Pobs or max(P_cal) < Pobs:
		final = -3
		return final
	
	# Shift index around Pobs to find 	intersection.
	iu  = where((P_cal[:-1] > Pobs) & (P_cal[1:] < Pobs)) 
	il  = where((P_cal[:-1] < Pobs) & (P_cal[1:] > Pobs))
	index = numpy.concatenate((iu,il),axis=1)
	index = index.reshape(numpy.size(index))
	
	
	for ip in range(numpy.size(index)):
		F_M2   	= 0
		F_M1   	= 0
		F_a    	= 0
		F_spin 	= 0
		a_spin  = 0.					# BH Spin initially					
		M2sim   = M2_cal[index[ip]]
		Mbh_sim = M1_cal[index[ip]]                        
		asim    = a_cal[index[ip]]
		Macc    = abs(Mbh-Mbh_sim)   	# Rest mass accreted by BH.
		#print 'Rest mass acrreted', Macc
		a_spin  = BH_spin(Mbh, Macc) 	# Spin of black hole with initial mass mbh.
						
	for property in b.Properties:
		if ((property == 'P')): # period of system
			continue
		elif (property == 'M2'): #Donor mass
			m_temp     = M2sim
			obs_temp   = b.Value['M2']
			obs_err    = b.Error['M2']
			if (m_temp > obs_temp-Nsigma*obs_err) and (m_temp < obs_temp+Nsigma*obs_err):
				F_M2 = 1
			continue
			 
		elif (property == 'Mbh'): # Mass of Black hole
			m_temp 	= Mbh_sim
			obs_temp 	= b.Value['Mbh']
			obs_err 	= b.Error['Mbh']
			if (m_temp > obs_temp-Nsigma*obs_err) and (m_temp < obs_temp+Nsigma*obs_err):
				F_M1 = 1
			continue
		
		elif (property == 'Aspin'): #BH spin
			m_temp = a_spin
			obs_temp = b.Value['Aspin']
			obs_err = b.Error['Aspin']
			#if (m_temp > obs_temp-Nsigma*obs_err) and (m_temp < obs_temp+Nsigma*obs_err):
			#The spin starts at 0, and must not exceed the upper limit, as the lower limit is negative, at least for LMC X-3.
			if (m_temp < obs_temp+Nsigma*obs_err):
				F_spin = 1
			
			continue

	if (F_M1 + F_M2 + F_spin) == 3:
		final = 1
	else:
		final = -3
  
	return final
