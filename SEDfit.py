from math import *
from numpy import *
import csv
from time import *
from scipy.spatial import ConvexHull
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.integrate import *
from amoeba import *
from matplotlib import pyplot as plt
import pyfits
import os
import platform
import smtplib

#Constants
NG=6.67384e-8 #Newton's Gravity in cm^3/g/s^2
R_sun=6.955e10 #Solar Radius in cm
M_sun=1.988435e33 #Solar Mass in g
L_sun=3.839e33 #Solar Luminosity in erg/s
pc=3.08567758e18 #1 parsec in cm
sigma_SB=5.6704e-5 #Stefan-Boltzmann constant in erg/cm^2/s/K^4
h=6.626e-27 #Planck's constant in cm^2*g/s
c=3e10 #Speed of light, cm/s
k=1.381e-16 #Boltzmann constant erg/K
g_scale=1.
g_points=0

#Which metallicity are we using for the PHOENIX models
use_Z='Z-0.0'

#stars=['HD37711','HD194789']
stars=['HD194789']

the_model='SED'
calc_lum=False

def fit():
	global tht,teff,logg,av,fixed,ig,star_dir,model_dir,star
	global phx_mu,phx_dict,teff_list,logg_list,str_teff_list,str_logg_list,phx_dir,phx_wav,filt_dict,zpf,phot_data,cwl
	global to_plot
	global dist,calc_lum
	global first_chi2
	global to_find_min
	starttime=time()
	to_find_min=[]
		
	#The number of iteration of amoeba we wish to run
	n_amo=1
	
	
	#This is the star we are looking at
	print 'Star: {}'.format(star)
	dist=1000/16.23
	print 'Using {} pc as the distance.'.format(dist)
	dist*=pc
	
	#This is the directory in which all relevant data files and model inputs are stored
	if platform.system()=='Windows':                                         
		star_dir='C:/Users/Jeremy/Dropbox/Programing/Astars/Stars/'+star+'/'
		model_dir='C:/Users/Jeremy/Dropbox/Programing/Astars/Stars/'+star+'/'+the_model+'/'
	if platform.system()=='Darwin':
		star_dir='/Users/jeremy/Dropbox/Programing/Astars/Stars/'+star+'/'
		model_dir='/Users/jeremy/Dropbox/Programing/Astars/Stars/'+star+'/'+the_model+'/'
	if platform.system()=='Linux':
		star_dir='/nfs/morgan/users/jones/Dropbox/Programing/Astars/Stars/'+star+'/'
		model_dir='/nfs/morgan/users/jones/Dropbox/Programing/Astars/Stars/'+star+'/'+the_model+'/'
	
	#Input File 
	phot_inp=model_dir+star+'.pi'
	#Output File
	phot_out=model_dir+star+'.po'
	
	#This reads the input file, star.pi
	data=dict()
	fixed=[]
	with open(phot_inp,'r') as input:
		input_reader=csv.reader(input,delimiter='\t')
		input_reader.next()
		for line in input_reader:
			data[line[0]]=float(line[1])
			if line[2] == 'fix':
				fixed.append(1)
			else:
				fixed.append(0)
					
	
	tht=data['tht']
	teff=data['teff']
	logg=data['logg']
	av=data['av']
	
	#Grabbing the photometry
	phot_inp=star_dir+star+'.phot'
	phot_data=dict()
	use_filts=[]
	with open(phot_inp,'r') as input:
		input_reader=csv.reader(input,delimiter='\t')
		input_reader.next()
		for line in input_reader:
			if line[0][0] != '#':
				phot_data[line[0]]=[line[1],line[2]]
				use_filts.append(line[0])
	read_phoenix() #Sets up the dictionary for use with the phoenix models
	read_filters(use_filts) #For the filters being used, this reads and stores the filter bandpass
	zpf=dict()	#Zero-point fluxes in erg/s/cm^2/cm
	zpf['2massH']=1.133e-10
	zpf['2massJ']=3.129e-10
	zpf['2massK']=4.283e-11
	zpf['cousinsI']=9.39e-10
	zpf['cousinsR']=1.92e-9
	zpf['johnsonB']=6.40e-9
	zpf['johnsonU']=4.34e-9
	zpf['johnsonV']=3.67e-9
	zpf['stromgrenb']=5.89e-9
	zpf['stromgrenu']=11.72e-9
	zpf['stromgrenv']=8.66e-9
	zpf['stromgreny']=3.73e-9
	zpf['johnsonR']=1.717e-9
	zpf['johnsonI']=8.499e-10
	zpf['johnsonJ']=3.14e-10
	zpf['johnsonH']=1.20e-10
	zpf['johnsonK']=4.12e-11
	zpf['F2740']=3.388e-9
	zpf['F2365']=3.388e-9
	zpf['F1965']=3.388e-9
	zpf['F1565']=3.388e-9
	zpf['wes15W']=3.6308e-9
	zpf['wes15N']=3.6308e-9
	zpf['wes18']=3.6308e-9
	zpf['wes22']=3.6308e-9
	zpf['wes25']=3.6308e-9
	zpf['wes33']=3.6308e-9

	r=[tht,teff,logg,av]
	ig=r
	range=[0.1,500.,0.1,0.1]
	#'''
	to_plot = False
	first_chi2=1e10
	first_chi2=1./sedfit(r)
	if sum(fixed) < 4:
		for i in arange(n_amo):
			#print 'Amoeba run {} out of {}.'.format(i+1,n_amo)
			#print '==========================================='
			r,chi2inv,nit=amoeba(r,range,sedfit)
			chi2=1./chi2inv
			#print 'Number of iterations: ',nit
	#'''
	
	to_plot = True
	calc_lum=True
	finalchi2,F_bol,L_bol=sedfit(r)
	calc_lum=False
	print '==========================================='
	print 'Final Chi^2: ',1./finalchi2
	print 'F_bol: ',F_bol
	print 'L_bol: ',L_bol
	
	fixins=['nofix','fix']
	
	with open(phot_out,'w') as output:
		output.write('\ntht\t{}\t{}'.format(tht,fixins[fixed[0]]))
		output.write('\nteff\t{}\t{}'.format(teff,fixins[fixed[1]]))
		output.write('\nlogg\t{}\t{}'.format(logg,fixins[fixed[2]]))
		output.write('\nav\t{}\t{}'.format(av,fixins[fixed[3]]))
	endtime=time()
	elapsed=(endtime-starttime)/60.
	print 'Done on {}, {} min elapsed'.format(ctime(time()),elapsed)
def get_errors():
	global tht,teff,logg,av,fixed,ig,star_dir,model_dir,star
	global phx_mu,phx_dict,teff_list,logg_list,str_teff_list,str_logg_list,phx_dir,phx_wav,filt_dict,zpf,phot_data,cwl
	global to_plot
	global dist,calc_lum
	global first_chi2
	global to_find_min
	starttime=time()
	#This is the star we are looking at
	print 'Star: {}'.format(star)
	
	#This is the directory in which all relevant data files and model inputs are stored
	if platform.system()=='Windows':                                         
		star_dir='C:/Users/Jeremy/Dropbox/Programing/Astars/Stars/'+star+'/'
		model_dir='C:/Users/Jeremy/Dropbox/Programing/Astars/Stars/'+star+'/'+the_model+'/'
	if platform.system()=='Darwin':
		star_dir='/Users/jeremy/Dropbox/Programing/Astars/Stars/'+star+'/'
		model_dir='/Users/jeremy/Dropbox/Programing/Astars/Stars/'+star+'/'+the_model+'/'
	if platform.system()=='Linux':
		star_dir='/nfs/morgan/users/jones/Dropbox/Programing/Astars/Stars/'+star+'/'
		model_dir='/nfs/morgan/users/jones/Dropbox/Programing/Astars/Stars/'+star+'/'+the_model+'/'
	
	#Input File 
	phot_inp=model_dir+star+'.po'
	#Output File
	phot_out=model_dir+star+'.perr'
	printfile = False
	
	#This reads the input file, star.pi
	data=dict()
	fixed=[]
	with open(phot_inp,'r') as input:
		input_reader=csv.reader(input,delimiter='\t')
		input_reader.next()
		for line in input_reader:
			data[line[0]]=float(line[1])
			if line[2] == 'fix':
				fixed.append(1)
			else:
				fixed.append(0)
					
	
	tht=data['tht']
	teff=data['teff']
	logg=data['logg']
	av=data['av']
	
	#Grabbing the photometry
	phot_inp=star_dir+star+'.phot'
	phot_data=dict()
	use_filts=[]
	with open(phot_inp,'r') as input:
		input_reader=csv.reader(input,delimiter='\t')
		input_reader.next()
		for line in input_reader:
			if line[0][0] != '#':
				phot_data[line[0]]=[line[1],line[2]]
				use_filts.append(line[0])
	read_phoenix() #Sets up the dictionary for use with the phoenix models
	read_filters(use_filts) #For the filters being used, this reads and stores the filter bandpass
	zpf=dict()	#Zero-point fluxes in erg/s/cm^2/cm
	zpf['2massH']=1.133e-10
	zpf['2massJ']=3.129e-10
	zpf['2massK']=4.283e-11
	zpf['cousinsI']=9.39e-10
	zpf['cousinsR']=1.92e-9
	zpf['johnsonB']=6.40e-9
	zpf['johnsonU']=4.34e-9
	zpf['johnsonV']=3.67e-9
	zpf['stromgrenb']=5.89e-9
	zpf['stromgrenu']=11.72e-9
	zpf['stromgrenv']=8.66e-9
	zpf['stromgreny']=3.73e-9
	zpf['johnsonR']=1.717e-9
	zpf['johnsonI']=8.499e-10
	zpf['johnsonJ']=3.14e-10
	zpf['johnsonH']=1.20e-10
	zpf['johnsonK']=4.12e-11
	zpf['F2740']=3.388e-9
	zpf['F2365']=3.388e-9
	zpf['F1965']=3.388e-9
	zpf['F1565']=3.388e-9
	zpf['wes15W']=3.6308e-9
	zpf['wes15N']=3.6308e-9
	zpf['wes18']=3.6308e-9
	zpf['wes22']=3.6308e-9
	zpf['wes25']=3.6308e-9
	zpf['wes33']=3.6308e-9

	r=[tht,teff,logg,av]
	itht=r[0]
	iteff=r[1]
	ilogg=r[2]
	iav=r[3]
	ig=[itht,iteff,ilogg,iav]
	to_plot = False
	first_chi2=1./sedfit(r)
	to_find_min=[]
	chi2_now=first_chi2
	diff=chi2_now-first_chi2
	#print 'Diff: ',diff
	
	#Do high error in angular radius
	print '===============Finding the upper bound of angular radius'
	if fixed[0] == 0:
		inc=0.1
		sign=1.
		while round(diff,3) != 1.000:
			r[0]+=inc*sign
			chi2_now=1./sedfit(r)
			diff=chi2_now-first_chi2
			#print 'Diff: ',diff
			if diff > 1.:
				if sign > 0.:
					inc*=0.1
				sign=-1.
			else:
				if sign < 0.:
					inc*=0.1
				sign=1.
	tht_hi=r[0]
	#Do low error in angular radius
	print '===============Finding the lower bound of angular radius'
	r=[itht,iteff,ilogg,iav]
	chi2_now=1./sedfit(r)
	diff=chi2_now-first_chi2
	if fixed[0] == 0:
		inc=0.1
		sign=-1.
		while round(diff,3) != 1.000:
			r[0]+=inc*sign
			chi2_now=1./sedfit(r)
			diff=chi2_now-first_chi2
			#print 'Diff: ',diff
			if diff > 1.:
				if sign < 0.:
					inc*=0.1
				sign=1.
			else:
				if sign > 0.:
					inc*=0.1
				sign=-1.
	tht_lo=r[0]
	#Do high error in effective temperature
	print '===============Finding the upper bound of effective temperature'
	r=[itht,iteff,ilogg,iav]
	chi2_now=1./sedfit(r)
	diff=chi2_now-first_chi2
	if fixed[1] == 0:
		inc=500.
		sign=1.
		while round(diff,3) != 1.000:
			r[1]+=inc*sign
			chi2_now=1./sedfit(r)
			diff=chi2_now-first_chi2
			#print 'Diff: ',diff
			if diff > 1.:
				if sign > 0.:
					inc*=0.1
				sign=-1.
			else:
				if sign < 0.:
					inc*=0.1
				sign=1.
	teff_hi=r[1]

	#Do low error in effective temperature
	print '===============Finding the lower bound of effective temperature'
	r=[itht,iteff,ilogg,iav]
	chi2_now=1./sedfit(r)
	diff=chi2_now-first_chi2
	if fixed[1] == 0:
		inc=500.
		sign=-1.
		while round(diff,3) != 1.000:
			r[1]+=inc*sign
			chi2_now=1./sedfit(r)
			diff=chi2_now-first_chi2
			#print 'Diff: ',diff
			if diff > 1.:
				if sign < 0.:
					inc*=0.1
				sign=1.
			else:
				if sign > 0.:
					inc*=0.1
				sign=-1.
	teff_lo=r[1]
	#Do high error in surface gravity
	print '===============Finding the upper bound of surface gravity'
	r=[itht,iteff,ilogg,iav]
	chi2_now=1./sedfit(r)
	diff=chi2_now-first_chi2
	if fixed[2] == 0:
		inc=0.1
		sign=1.
		while round(diff,3) != 1.000:
			r[2]+=inc*sign
			chi2_now=1./sedfit(r)
			diff=chi2_now-first_chi2
			#print 'Diff: ',diff
			if diff > 1.:
				if sign > 0.:
					inc*=0.1
				sign=-1.
			else:
				if sign < 0.:
					inc*=0.1
				sign=1.
			if sign > 0. and r[2] > 5.999:
				diff = 1.
				r[2] = 5.999
	logg_hi=r[2]

	#Do low error in surface gravity
	print '===============Finding the lower bound of surface gravity'
	r=[itht,iteff,ilogg,iav]
	chi2_now=1./sedfit(r)
	diff=chi2_now-first_chi2
	if fixed[2] == 0:
		inc=0.1
		sign=-1.
		while round(diff,3) != 1.000:
			r[2]+=inc*sign
			chi2_now=1./sedfit(r)
			diff=chi2_now-first_chi2
			#print 'Diff: ',diff
			if diff > 1.:
				if sign < 0.:
					inc*=0.1
				sign=1.
			else:
				if sign > 0.:
					inc*=0.1
				sign=-1.
			if sign < 0. and r[2] < 2.001:
				diff = 1.
				r[2] = 2.001
	logg_lo=r[2]
	#Do high error in reddening
	print '===============Finding the upper bound of reddening'
	r=[itht,iteff,ilogg,iav]
	chi2_now=1./sedfit(r)
	diff=chi2_now-first_chi2
	if fixed[3] == 0:
		inc=0.1
		sign=1.
		while round(diff,3) != 1.000:
			r[3]+=inc*sign
			chi2_now=1./sedfit(r)
			diff=chi2_now-first_chi2
			#print 'Diff: ',diff
			if diff > 1.:
				if sign > 0.:
					inc*=0.1
				sign=-1.
			else:
				if sign < 0.:
					inc*=0.1
				sign=1.
	av_hi=r[3]

	#Do low error in reddening
	print '===============Finding the lower bound of reddening'
	r=[itht,iteff,ilogg,iav]
	chi2_now=1./sedfit(r)
	diff=chi2_now-first_chi2
	if fixed[3] == 0:
		inc=0.1
		sign=-1.
		while round(diff,3) != 1.000:
			r[3]+=inc*sign
			chi2_now=1./sedfit(r)
			diff=chi2_now-first_chi2
			#print 'Diff: ',diff
			if diff > 1.:
				if sign > 0.:
					inc*=0.1
				sign=-1.
			else:
				if sign < 0.:
					inc*=0.1
				sign=1.
			if sign < 0. and r[3] < 0.:
				diff=1.
				r[3]=0.
				#print 'Lower bound of reddening is zero'
	av_lo=r[3]
	
	calc_lum=True
	finalchi2,F_bol,L_bol=sedfit(ig)
	all_Fs=[F_bol]
	all_Cs=[1./finalchi2]
	if fixed[0] == 0:
		chi2_thtlo,F_thtlo,L_thtlo=sedfit([tht_lo,ig[1],ig[2],ig[3]])
		all_Fs.append(F_thtlo)
		all_Cs.append(1./chi2_thtlo)
		chi2_ththi,F_ththi,L_ththi=sedfit([tht_hi,ig[1],ig[2],ig[3]])
		all_Fs.append(F_ththi)
		all_Cs.append(1./chi2_ththi)
	if fixed[1] == 0:
		chi2_tefflo,F_tefflo,L_tefflo=sedfit([ig[0],teff_lo,ig[2],ig[3]])
		all_Fs.append(F_tefflo)
		all_Cs.append(1./chi2_tefflo)
		chi2_teffhi,F_teffhi,L_teffhi=sedfit([ig[0],teff_hi,ig[2],ig[3]])
		all_Fs.append(F_teffhi)
		all_Cs.append(1./chi2_teffhi)
	if fixed[2] == 0:
		chi2_logglo,F_logglo,L_logglo=sedfit([ig[0],ig[1],logg_lo,ig[3]])
		all_Fs.append(F_logglo)
		all_Cs.append(1./chi2_logglo)
		chi2_logghi,F_logghi,L_logghi=sedfit([ig[0],ig[1],logg_hi,ig[3]])
		all_Fs.append(F_logghi)
		all_Cs.append(1./chi2_logghi)
	if fixed[3] == 0:
		chi2_avlo,F_avlo,L_avlo=sedfit([ig[0],ig[1],ig[2],av_lo])
		all_Fs.append(F_avlo)
		all_Cs.append(1./chi2_avlo)
		chi2_avhi,F_avhi,L_avhi=sedfit([ig[0],ig[1],ig[2],av_hi])
		all_Fs.append(F_avhi)
		all_Cs.append(1./chi2_avhi)
	calc_lum=False
	for i in range(len(all_Fs)):
		print all_Fs[i],all_Cs[i]
	
	print 'Low F, F_bol, High F'
	print min(all_Fs),F_bol,max(all_Fs)

	fixins=['nofix','fix']
	with open(phot_out,'w') as output:
		output.write('\ntht\t{}\t+{}\t-{}\t{}'.format(itht,tht_hi-itht,itht-tht_lo,fixins[fixed[0]]))
		output.write('\nteff\t{}\t+{}\t-{}\t{}'.format(iteff,teff_hi-iteff,iteff-teff_lo,fixins[fixed[1]]))
		output.write('\nlogg\t{}\t+{}\t-{}\t{}'.format(ilogg,logg_hi-ilogg,ilogg-logg_lo,fixins[fixed[2]]))
		output.write('\nav\t{}\t+{}\t-{}\t{}'.format(iav,av_hi-iav,iav-av_lo,fixins[fixed[3]]))
	endtime=time()
	elapsed=(endtime-starttime)/60.
	print 'Done on {}, {} min elapsed'.format(ctime(time()),elapsed)
	return to_find_min

def sedfit(p,data=None):
	global tht,teff,logg,av,fixed,ig,star_dir,model_dir,star
	global phx_mu,phx_dict,teff_list,logg_list,str_teff_list,str_logg_list,phx_dir,phx_wav,filt_dict,zpf,phot_data,cwl
	global to_plot
	global dist,calc_lum
	global first_chi2
	global to_find_min
	ministarttime=time()
	try:	
		if fixed[0] == 0:
			tht=p[0]
		else:
			tht=ig[0]
		if fixed[1] == 0:
			teff=p[1]
		else:
			teff=ig[1]
		if fixed[2] == 0:
			logg=p[2]
		else:
			logg=ig[2]
		if fixed[3] == 0:
			av=abs(p[3])
		else:
			av=ig[3]
		
		n_params=4-fixed[0]-fixed[1]-fixed[2]-fixed[3]
		
		#Convert tht from radius in mas to diameter in radians
		tht/=2.06265e8
		tht*=2.
		
		#Caulculating the unreddened, unscaled SED (sed_urus)
		if model == 'phoenix':
			sed_urus=extract_phoenix(teff,logg)
		if model == 'bb':
			sed_urus=2.*h*c**2./(phx_wav)**5.*1./(exp(h*c/k/teff/(phx_wav))-1)
		#Calculating the unreddened SED (sed_ur)
		sed_ur=sed_urus*tht**2.*pi/2.*tht/abs(tht) #The tht/abs(tht) part is to prevent negative thetas
		
		#Calculating the SED (sed)
		x=1./(10000.*phx_wav)	#wave-number in um^-1
		y=x-1.82
		r_v=3.1
		a_ext=[]
		b_ext=[]
		for i in range(len(phx_wav)):
			if phx_wav[i] < 9.091e-4:	#For lambda <~0.9 um. From O'Donnell, 1994, ApJ, 422, 158O
				a_ext.append(1.+0.104*y[i]-0.609*y[i]**2.+0.701*y[i]**3.+1.137*y[i]**4.-1.718*y[i]**5.-0.827*y[i]**6.+1.647*y[i]**7.-0.505*y[i]**8.)
				b_ext.append(1.952*y[i]+2.908*y[i]**2.-3.989*y[i]**3.-7.985*y[i]**4.+11.102*y[i]**5.+5.491*y[i]**6.-10.805*y[i]**7.+3.347*y[i]**8.)
			if phx_wav[i] >=9.091e-4:	#For lambda >=~0.9 um. From Cardelli, Clayton, Mathis, 1989, ApJ, 345, 245
				a_ext.append(0.574*x[i]**1.61)
				b_ext.append(-0.574*x[i]**1.61)
		a_ext=array(a_ext)
		b_ext=array(b_ext)
		alam=(a_ext+b_ext/r_v)*av
		sed=sed_ur*10**(-0.4*alam)
		
		sed_filts=dict()
		for f in filt_dict:
			sed_filts[f]=sed*filt_dict[f]/fwhm(phx_wav,filt_dict[f])/1e8
		filt_fluxes=dict()
		for f in filt_dict:
			filt_fluxes[f]=trapz(sed_filts[f],x=phx_wav)
		phot_dict=dict()
		for f in filt_dict:
			phot_dict[f]=-2.5*log10(filt_fluxes[f]/zpf[f])
			
		phot_diff=[]
		phot_err=[]
		for f in filt_dict:
			phot_diff.append(filt_fluxes[f]-zpf[f]*10.**(-0.4*double(phot_data[f][0])))
			phot_err.append(double(phot_data[f][1])*zpf[f]*0.4*log(10.)*10.**(-0.4*double(phot_data[f][0])))
		phot_diff=array(phot_diff)
		phot_err=array(phot_err)
		if len(phot_diff)-n_params <= 0:
			print 'It would be best if you had either fewer free parameters or more data points. Changing how we calculate chi^2...'
			chi2=sum(phot_diff**2./phot_err**2.)/(float(len(phot_diff)))
		else:
			chi2=sum(phot_diff**2./phot_err**2.)/(float(len(phot_diff))-n_params)
	
	
	
		tht*=2.06265e8
		tht/=2.
		miniendtime=time()
		minielapsed=miniendtime-ministarttime
		if chi2 < first_chi2:
			print 'Angular Radius: {}, Effective Temperature: {}, Surface Gravity: {}, Reddening: {}, Chi^2: {}'.format(tht,teff,logg,av,chi2)
		else:
			print 'Chi^2 ({}) is larger than first chi^2 ({})'.format(chi2,first_chi2)
		if chi2 < first_chi2-0.001:
			to_find_min.append([chi2,[tht,teff,logg,av]])
		#print '{} seconds elapsed.'.format(minielapsed)
		if to_plot:
			plot_sed=model_dir+star+'_disksed.pdf'
			print 'Saving plot to {}.'.format(plot_sed)
			plt.plot(phx_wav,sed/1e8,color='0.4',linestyle='-')
			print 'Filter   Model Phot   Measured Phot'
			print '============================'
			plot_fluxes=[]
			for f in filt_dict:
				print f+':',phot_dict[f],phot_data[f][0]
				this_flux=zpf[f]*10.**(-0.4*double(phot_data[f][0]))
				this_flux_err=double(phot_data[f][1])*zpf[f]*0.4*log(10.)*10.**(-0.4*double(phot_data[f][0]))
				plt.plot([cwl[f],cwl[f]],[this_flux+this_flux_err,this_flux-this_flux_err],'k-')
				plt.plot(cwl[f], filt_fluxes[f],'bs')
				plt.plot(cwl[f], this_flux,'ro')
				plot_fluxes.append(filt_fluxes[f])
				plot_fluxes.append(this_flux)
			plt.xlabel('Wavelength (cm)')
			plt.xscale('log')
			plt.ylabel('Flux (erg/s/cm^2/A)')
			plt.yscale('log')
			#plt.show()
			plt.xlim(1e-5,2.6e-4)
			'''
			sed_inrange=sed[(phx_wav >= 1e-5).nonzero()]
			new_wave=phx_wav[(phx_wav >= 1e-5).nonzero()]
			sed_inrange=sed_inrange[(new_wave <= 2.6e-4).nonzero()]
			sed_inrange/=1e8
			#plt.ylim(1e-12,5e-10)
			plt.ylim(amin(sed_inrange)*0.9,amax(sed_inrange)*1.1)
			'''
			plt.ylim(min(plot_fluxes)*0.9,max(plot_fluxes)*1.1)
			plt.savefig(plot_sed)
			plt.close()
		if calc_lum:
			F_bol=trapz(sed,x=phx_wav)
			#print 'F_bol: {} erg/s/cm^2'.format(F_bol)
			L_bol=4.*pi*dist**2.*F_bol/L_sun
			#print 'L_bol: {} L_sun'.format(L_bol)
			return 1./chi2,F_bol,L_bol
		
		if chi2 > 0:
			return 1./chi2
		else:
			print 'Angular Radius: {}, Effective Temperature: {}, Surface Gravity: {}, Reddening: {}, Chi^2: {}'.format(tht,teff,logg,av,chi2)
			return 1e-20
	except:
		if calc_lum:
			return 1e-20,0.,0.
			#print 'High Chi^2'
		return 1e-20
	
def read_phoenix():
	global phx_mu,phx_dict,teff_list,logg_list,str_teff_list,str_logg_list,phx_dir,phx_wav,filt_dict,cwl
	
	if platform.system()=='Windows':
		phx_dir='C:/Users/Jeremy/Phoenix_Spectra/'+use_Z+'/'
	if platform.system()=='Darwin':
		phx_dir='/Users/jeremy/Phoenix_Spectra/'+use_Z+'/'
	if platform.system()=='Linux':
		phx_dir='/nfs/morgan/users/jones/Phoenix_Spectra/'+use_Z+'/'
	file='lte07000-4.00'+use_Z[1:]+'.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits'
	hdulist = pyfits.open(phx_dir+file)
	first_file = hdulist[0].data
	phx_mu=hdulist[1].data
	#print hdulist[0].header
	hdulist.close()
	phx_mu=list(phx_mu)
	phx_mu=array(phx_mu)
	first_file=list(first_file)
	first_file=array(first_file)
	phx_wav=(arange(len(first_file[-1]))+500.)*1e-8
	phx_dict={file:store_phoenix(first_file)}

	teff_list=[]
	logg_list=[]
	str_teff_list=[]
	str_logg_list=[]

	teff_list_1=arange(27)*100+2300
	teff_list_2=arange(20)*100+5100
	teff_list_3=arange(25)*200+7200
	teff_list=concatenate((teff_list_1,teff_list_2,teff_list_3))
	str_teff_list=[]
	for i in range(len(teff_list)):
		if teff_list[i] >= 10000:
			str_teff_list.append(str(teff_list[i]))
		else:
			str_teff_list.append('0'+str(teff_list[i]))
	logg_list=arange(13)*0.5
	str_logg_list=[]
	for i in range(len(logg_list)):
		if i == 0:
			str_logg_list.append('+'+str(logg_list[i])+'0')
		else:
			str_logg_list.append('-'+str(logg_list[i])+'0')

def store_phoenix(inarr):
	global phx_mu,phx_dict,teff_list,logg_list,str_teff_list,str_logg_list,phx_dir,phx_wav,filt_dict,cwl
	almost_outarr=[]
	for i in range(len(phx_mu)):
		almost_outarr.append(phx_mu[i]*inarr[i])
	almost_outarr=array(almost_outarr)
	outarr=[]
	for i in range(len(phx_wav)):
		outarr.append(trapz(almost_outarr[:,i],x=phx_mu))
	return array(outarr)
	
def extract_phoenix(teff,logg):
	global phx_mu,phx_dict,teff_list,logg_list,str_teff_list,str_logg_list,phx_dir,phx_wav,filt_dict,cwl
	#print teff,logg
	path_list=os.listdir(phx_dir)
	ftp_dir='ftp://phoenix.astro.physik.uni-goettingen.de/SpecIntFITS/PHOENIX-ACES-AGSS-COND-SPECINT-2011/'+use_Z+'/'
	
	gett=teff_list[(teff_list >= teff).nonzero()]
	lett=teff_list[(teff_list <= teff).nonzero()]
	tlo=max(lett)
	thi=min(gett)
	gegg=logg_list[(logg_list >= logg).nonzero()]
	legg=logg_list[(logg_list <= logg).nonzero()]
	glo=max(legg)
	ghi=min(gegg)

	#print str_teff_list[(teff_list == teff).nonzero()[0]]
	tlo_str=str_teff_list[(teff_list == tlo).nonzero()[0]]
	glo_str=str_logg_list[(logg_list == glo).nonzero()[0]]
	thi_str=str_teff_list[(teff_list == thi).nonzero()[0]]
	ghi_str=str_logg_list[(logg_list == ghi).nonzero()[0]]

	ll_file='lte'+tlo_str+glo_str+use_Z[1:]+'.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits'
	lh_file='lte'+tlo_str+ghi_str+use_Z[1:]+'.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits'
	hl_file='lte'+thi_str+glo_str+use_Z[1:]+'.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits'
	hh_file='lte'+thi_str+ghi_str+use_Z[1:]+'.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits'

	#Low Temperature, Low Gravity
	if ll_file in phx_dict:
		store_ll = phx_dict[ll_file]
	else:
		#print 'Grabbing file {}'.format(tlo_str+glo_str)
		try:
			ll_hdulist = pyfits.open(phx_dir+ll_file)
		except:
			print 'Grabbing file {}'.format(tlo_str+glo_str)
			print '... from the Internet'
			ll_hdulist = pyfits.open(ftp_dir+ll_file)
			print 'Grabbed'
		ll=ll_hdulist[0].data
		ll=list(ll)
		ll=array(ll)
		ll_hdulist.close()
		store_ll=store_phoenix(ll)
		phx_dict[ll_file] = store_ll
	#Low Temperature, High Gravity
	if lh_file in phx_dict:
		store_lh = phx_dict[lh_file]
	else:
		#print 'Grabbing file {}'.format(tlo_str+ghi_str)
		try:
			lh_hdulist = pyfits.open(phx_dir+lh_file)
		except:
			print 'Grabbing file {}'.format(tlo_str+ghi_str)
			print '... from the Internet'
			lh_hdulist = pyfits.open(ftp_dir+lh_file)
			print 'Grabbed'
		lh=lh_hdulist[0].data
		lh=list(lh)
		lh=array(lh)
		lh_hdulist.close()
		store_lh=store_phoenix(lh)
		phx_dict[lh_file] = store_lh
	#High Temperature, Low Gravity
	if hl_file in phx_dict:
		store_hl = phx_dict[hl_file]
	else:
		#print 'Grabbing file {}'.format(thi_str+glo_str)
		try:
			hl_hdulist = pyfits.open(phx_dir+hl_file)
		except:
			print 'Grabbing file {}'.format(thi_str+glo_str)
			print '... from the Internet'
			hl_hdulist = pyfits.open(ftp_dir+hl_file)
			print 'Grabbed'
		hl=hl_hdulist[0].data
		hl=list(hl)
		hl=array(hl)
		hl_hdulist.close()
		store_hl=store_phoenix(hl)
		phx_dict[hl_file] = store_hl
	#High Temperature, High Gravity
	if hh_file in phx_dict:
		store_hh = phx_dict[hh_file]
	else:
		#print 'Grabbing file {}'.format(thi_str+ghi_str)
		try:
			hh_hdulist = pyfits.open(phx_dir+hh_file)
		except:
			print 'Grabbing file {}'.format(thi_str+ghi_str)
			print '... from the Internet'
			hh_hdulist = pyfits.open(ftp_dir+hh_file)
			print 'Grabbed'
		hh=hh_hdulist[0].data
		hh=list(hh)
		hh=array(hh)
		hh_hdulist.close()
		store_hh=store_phoenix(hh)
		phx_dict[hh_file] = store_hh
	
	if thi != tlo and ghi != glo:
		interpolated_flux=store_ll/(thi-tlo)/(ghi-glo)*(thi-teff)*(ghi-logg)+store_lh/(thi-tlo)/(ghi-glo)*(thi-teff)*(logg-glo)+store_hl/(thi-tlo)/(ghi-glo)*(teff-tlo)*(ghi-logg)+store_hh/(thi-tlo)/(ghi-glo)*(teff-tlo)*(logg-glo)
	elif thi != tlo and ghi == glo:
		interpolated_flux=store_ll+(store_hl-store_ll)*(teff-tlo)/(thi-tlo)
	elif thi == tlo and ghi != glo:		
		interpolated_flux=store_ll+(store_lh-store_ll)*(logg-glo)/(ghi-glo)
	elif thi == tlo and ghi == glo:
		interpolated_flux=store_ll
	return array(interpolated_flux)
def read_filters(filts):
	global phx_mu,phx_dict,teff_list,logg_list,str_teff_list,str_logg_list,phx_dir,phx_wav,filt_dict
	global cwl
	if platform.system()=='Windows':
		filt_dir='C:/Users/Jeremy/Dropbox/Programing/Astars/Band_Passes/'
	if platform.system()=='Darwin':
		filt_dir='/Users/jeremy/Dropbox/Programing/Astars/Band_Passes/'
	if platform.system()=='Linux':
		filt_dir='/nfs/morgan/users/jones/Dropbox/Programing/Astars/Band_Passes/'
	filt_dict=dict()
	cwl=dict()	#Central Wavelengths in cm
	cwl['stromgrenu']=3.491e-5
	cwl['johnsonU']=3.735e-5
	cwl['stromgrenv']=4.111e-5
	cwl['johnsonB']=4.443e-5
	cwl['stromgrenb']=4.662e-5
	cwl['stromgreny']=5.456e-5
	cwl['johnsonV']=5.483e-5
	cwl['cousinsR']=6.855e-5
	cwl['cousinsI']=8.060e-5
	cwl['johnsonJ']=1.25e-4
	cwl['johnsonH']=1.65e-4
	cwl['johnsonK']=2.2e-4
	cwl['2massJ']=1.235e-4
	cwl['2massH']=1.662e-4
	cwl['2massK']=2.159e-4
	cwl['F2740']=2.74e-5
	cwl['F2365']=2.365e-5
	cwl['F1965']=1.965e-5
	cwl['F1565']=1.565e-5
	cwl['wes15W']=1.545e-5
	cwl['wes15N']=1.549e-5
	cwl['wes18']=1.799e-5
	cwl['wes22']=2.2e-5
	cwl['wes25']=2.493e-5
	cwl['wes33']=3.294e-5

	for i in range(len(filts)):
		this_wav=[]
		this_tran=[]
		with open(filt_dir+filts[i]+'.txt') as input:
			try:
				input_reader=csv.reader(input,delimiter='\t')
				input_reader.next()
				for line in input_reader:
					this_wav.append(float(line[0])/1e8)
					this_tran.append(float(line[1]))
			except:
				input_reader=csv.reader(input,delimiter=' ')
				input_reader.next()
				for line in input_reader:
					this_wav.append(float(line[0])/1e8)
					this_tran.append(float(line[1]))
		int_trans=zeros(len(phx_wav))
		f=interp1d(this_wav,this_tran)
		for j in range(len(phx_wav)):
			if phx_wav[j] >= this_wav[0] and phx_wav[j] <= this_wav[-1]:
				int_trans[j]=f(phx_wav[j])
		abs_wav_minus_cwl=abs(phx_wav - cwl[filts[i]])
		cwl_trans=int_trans[(abs_wav_minus_cwl == min(abs_wav_minus_cwl)).nonzero()]
		#cwl_trans=int_trans[(phx_wav > cwl[filts[i]]-0.5e-8).nonzero()]
		#tmp_wav=phx_wav[(phx_wav > cwl[filts[i]]-0.5e-8).nonzero()]
		#cwl_trans=cwl_trans[(tmp_wav < cwl[filts[i]]+0.5e-8).nonzero()]
		int_trans=int_trans/cwl_trans[0]
		filt_dict[filts[i]]=int_trans

def unitrange(res):
	return arange(res+1)/float(res)
def fwhm(wave,transmission):
	test_trans=transmission[(transmission > 0.01).nonzero()]
	test_trans=test_trans[(test_trans < 0.99).nonzero()]
	if len(test_trans) == 0:
		ones_wave=wave[(transmission > 0.01).nonzero()]
		fwhm_high=max(ones_wave)
		fwhm_low=min(ones_wave)
		fwhm=fwhm_high-fwhm_low
	else:
		close_to_one=min((transmission-1)**2.)
		max_wave=wave[((transmission-1)**2.==close_to_one).nonzero()]
		max_wave=max_wave[len(max_wave)/2]
		wave_high=wave[(wave > max_wave).nonzero()]
		trans_high=transmission[(wave > max_wave).nonzero()]
		close_to_half=min((trans_high-0.5)**2.)
		fwhm_high=wave_high[((trans_high-0.5)**2.==close_to_half).nonzero()]
		wave_low=wave[(wave < max_wave).nonzero()]
		trans_low=transmission[(wave < max_wave).nonzero()]
		close_to_half=min((trans_low-0.5)**2.)
		fwhm_low=wave_low[((trans_low-0.5)**2.==close_to_half).nonzero()]
		fwhm=fwhm_high-fwhm_low	
	return fwhm

model='phoenix'
#model='bb'
print 'Using {} models'.format(model)

for j in range(len(stars)):
	star = stars[j]
	pi_file='C:/Users/Jeremy/Dropbox/Programing/Astars/Stars/'+star+'/'+the_model+'/'+star+'.pi'
	
	on=True
	while on:
		print 'on = True'
		fit()
		to_find_min=get_errors()
		if not to_find_min:
			on=False
		else:
			chi2s=[]
			rs=[]
			pi_list=['tht','teff','logg','av']
			for i in range(len(to_find_min)):
				chi2s.append(to_find_min[i][0])
				rs.append(to_find_min[i][1])
			chi2s=array(chi2s)
			rs=array(rs)
			print min(chi2s)
			print chi2s
			print rs
			replace_r=rs[(chi2s == min(chi2s)).nonzero()][0]
			print replace_r
			
			print 'Rewriting .pi file.'
			for i in range(len(pi_list)):
				print '{}\t{}'.format(pi_list[i],replace_r[i])
				if i == 0:
					open(pi_file,'w').write('\n{}\t{}\t{}'.format(pi_list[i],replace_r[i],'nofix'))
				else:
					open(pi_file,'a').write('\n{}\t{}\t{}'.format(pi_list[i],replace_r[i],'nofix'))
		print '==================='


