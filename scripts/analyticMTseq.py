# -*- coding: utf-8 -*-
"""
Created on Tue Feb  23 18:20:15 2015

@author: Mads Sørensen
"""
import numpy, sys, os
import scipy.optimize
from numpy import where as where
#import MTseq_npy

def evolve_orbital_period_RLO(alpha, beta, P0, M10, M20, M2):
    # Analytical solution of a lagrangian formulation of a binary system of two point masses M1 and M2
    # where alpha and beta are coefficients to describe fractions of mass lost and transfered.
    
    # Input ::
    # Alpha ::    Float, the mass fraction lost from the system at the donor
    # Beta  ::    Float, the mass fraction lost from system at the accretor.
    # P0    ::    Float, initial orbital period
    # M10   ::    Float, initial accretor mass
    # M20   ::    Float, initial donor mass
    # M2    ::    Array, Float, how M2 evolves.

    # The unit of P0 dictates the unit of P.
	
    M0 = M10 + M20
    M1 = -(1-alpha)*(1-beta)*(M2-M20)+M10
    M  = M1 + M2
    if beta < 1.0:
    	P  = P0*(M2/M20)**(3.*(alpha-1.)) * (M1/M10)**(3./(beta-1.))*(M0/M)**2.
    elif beta == 1.0:
    	P  = P0*(M2/M20)**(3.*(alpha-1.)) * numpy.exp( 3.*(1.-alpha)/M1 * (M2-M20) ) * (M0/M)**2.
    return P

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


def BH_spin_non_zero(Mrlo, a_natal):
#Assuming the BH has some initial/natal spin $a$ before it begins accretion, what initial mass $M_i$
#is this responding to according to the paper Thorne 1973 on BH spin from accretion.
#We have assumed $M_i0$ = 0.
    def F(Mi):
        dum = 3.* Mi * numpy.sin( (Mrlo-Mi)/(3.*Mi) + numpy.arcsin(1./3.) )
        return -a_natal+(2/3.)**0.5*Mi/dum*(4.- numpy.sqrt(18.*Mi**2./(dum**2.) -2 ))
    Mi = scipy.optimize.newton(F, Mrlo)
    Mi = float(Mi)
    return Mi

def BH_w_spin(Mrlo, arlo, M0):
#Find the spin of a BH: given some natal spin, a specific BH mass at RLO and some accreted mass M0:
    Mi = BH_spin_non_zero(Mrlo, arlo)
    M0i = Mrlo - Mi
    a = BH_spin(Mi,M0i+M0)
    print Mi, Mrlo, arlo, a, M0
    return a

class ObservedSystem(object):
	def __init__(self, SystemName):
		self.SystemName = SystemName
		obs_data_path = '/Users/MS/Dropbox/projects/xrb_submission_script/systems/'
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

def LagrangianMTseq(alpha, beta, M20, Mbh, P, system):
	M20 = float(M20)
	Mbh = float(Mbh)
	P = float(P)
	
	Nsigma = 2.0 		# Number of accepted standard deviations away observed value.
	
	#Read in system values to compare with
	b = ObservedSystem(system)
	Pobs = b.Value['P']
	M2low 	= 1.0

	# Array with accepted models
	F_M2       = 0
	F_M1       = 0
	F_spin 	   = 0
	F_q 	   = 0
	
	M2_cal 	= numpy.linspace(M20, M2low, 5001)
	M1_cal 	= -(1-alpha)*(1-beta)*(M2_cal-M20)+Mbh
	P_cal 	= evolve_orbital_period_RLO(alpha, beta, P, Mbh, M20, M2_cal)

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
		F_spin 	= 0
		F_q     = 0
		a_spin  = 0.					# BH Spin initially					
		M2sim   = M2_cal[index[ip]]
		Mbh_sim = M1_cal[index[ip]]                        
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
			if (m_temp >= obs_temp-Nsigma*obs_err) and (m_temp <= obs_temp+Nsigma*obs_err):
				F_M2 = 1
			continue
		elif (property == 'q'): #Mass ratio
			m_temp     = Mbh_sim/M2sim
			obs_temp   = b.Value['q']
			obs_err    = b.Error['q']
			if (m_temp >= obs_temp-Nsigma*obs_err) and (m_temp <= obs_temp+Nsigma*obs_err):
				F_q = 1
			continue
			 
		elif (property == 'Mbh'): # Mass of Black hole
			m_temp 	= Mbh_sim
			obs_temp 	= b.Value['Mbh']
			obs_err 	= b.Error['Mbh']
			if (m_temp >= obs_temp-Nsigma*obs_err) and (m_temp <= obs_temp+Nsigma*obs_err):
				F_M1 = 1
			continue
		
		elif (property == 'Aspin'): #BH spin
			m_temp = a_spin
			obs_temp = b.Value['Aspin']
			obs_err = b.Error['Aspin']
			#if (m_temp > obs_temp-Nsigma*obs_err) and (m_temp < obs_temp+Nsigma*obs_err):
			#The spin starts at 0, and must not exceed the upper limit, as the lower limit is negative, at least for LMC X-3.
			if (m_temp <= obs_temp+Nsigma*obs_err):
				F_spin = 1
			
			continue

	if (F_M1 + F_q + F_spin) == 3:
		final = 1
	else:
		final = -3
  
	return final