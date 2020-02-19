# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 09:49:57 2016

@author: Roberto

This code wraps the executable batch file of SMARTS2 to produce a series of 
solar spectra based on given parameters such as time point, geographic location,
angle of tilt of the array, angle of azimuth of the array, air temperature, 
relative humidity, among others

Written in Python 3.5.1

"""

import datetime as dt # builtin 
from os import remove, chdir, getcwd # builtin 
from os.path import join as osjoin # builtin 
from os.path import isfile # builtin 
from subprocess import call # builtin 
import time as T # builtin 
import pandas as pd # '0.18.0'


def SMARTS2(timepoint, press, air_t, rel_h, av_t, aer_m, til_ang, 
		azh_ang, lim_i, lim_s, folder, t_file,
		spectral_only = False, scalar_only = False):
		
	'''
	This function runs the executable batch file of SMARTS2, which is available 
	in: https://www.nrel.gov/grid/solar-resource/smarts.html
	This changes a template file and excecutes the batch file, which provides 
	the possibility to produce in series several solar spectral models based on 
	a set of parameters:
	- timepoint		:	should be a datetime type data.
	- press			: 	atmospheric pressure measured in mbar.
	- air_t			: 	air temperature in centigrade degrees.
	- rel_h			: 	relative humidity with this the precipitable water is 
						calculated, its measured as a percentage
	- av_t			: 	average temperature of the region.
	- aer_m			:	model of aerosol these are the ones available in 
						SMARTS2, see the instructive or read the paper for more 
						information.
	- til_ang		:	tilt angle of the array receiving the modelled solar 
						irradiance.
	- azh_ang		:	azimuth angle of the array receiving the modelled solar 
						irradiance.
	- lim_i			: 	lower limit of the modelled spectrum.
	- lim_s			:	upper limit of the modelled spectrume
	- folder		: 	folder in which the batch SMARTS2 file and the template 
						file are found.
	- t_file		: 	name of the template file.
	- spectral_only	: 	True for only spectral outcomes of the model, the default 
						is False
	- scalar_only	:	True to return as well scalar outcomes of the model, the 
						default is False.
	'''
			
	smarts_ext = 'smarts295.ext.txt'
	smarts_out = 'smarts295.out.txt'
	c_fol = getcwd() # Variable to save the current working folder
	chdir(folder) # Sets the working directory to folder
	path = osjoin(folder, t_file)
	
	with open(path, 'r+') as param_template: # Opens the template input file
		
		# Converts the templte file from bits to string
		content = param_template.read() 
		
		# Extracting all time data from the datetime data
		yr = timepoint.year 
		mt = timepoint.month 
		dy = timepoint.day 
		hr = timepoint.hour
		mn = timepoint.minute 
		sn = timepoint.second
		
		# Produces a float with the time in decimal form
		time =  hr + (mn / 60) + (sn / 3600)
		
		# These are the parameters to change
		params = [['Press', press], ['Air_temp', air_t], ['Rel_hum', rel_h], 
			['Av_temp', av_t], ['Aer_mod', aer_m], ['Til_ang',  til_ang], 
			['Azh_ang', azh_ang], ['Lim_inf', lim_i], ['Lim_sup', lim_s], 
			['Year', yr], ['Month', mt], ['Day', dy], 
			['Time', time], ['Spc_rsl', '2'], ['Spc_out', '1 8']] 
		
		for par, val in params: 
			# This is changing the parameters in the template 
			content = content.replace(par, str(val))
		
		for n in range(10):
			try:
				# This is saving the edited template as inputfile
				with open('smarts295.inp.txt', 'w') as inputfile: 
					inputfile.write(content)
				break
			
			except Exception: # in case of error
				print(
				'Possible Error at: {}, {}, due not able to read \
				input file. Attempt {}'.format(
				timepoint, aer_m, n))
				T.sleep(1)
				continue
	
	try:
		# Delete output files to prevent misplaced data
		path = osjoin(folder, smarts_ext)
		remove(path) 
			
	except FileNotFoundError:
		print("File: {} couldn't be deleted".format(smarts_ext))
		pass
	
	try:
		path = osjoin(folder, smarts_out)
		remove(path)
	
	except FileNotFoundError:
		print("File: {} couldn't be deleted".format(smarts_out))
		pass
	
	ext_file = False
	slp = 0.1
	
	while not ext_file:
		txt_end = '\r'
		# This run the executable batch file
		batch = 'smarts295bat.exe'
		path = osjoin(folder, batch)
		call(path, shell = True) 
		T.sleep(slp) # gives time for the process to be completed before next step
		path = osjoin(folder, smarts_ext)
		ext_file = isfile(path) # This confirms the creation of the output file
		
		if not ext_file: 
			print('\rOutput file for {} not found reprocessing and \
			waiting {} secs'.format(timepoint, slp), ' ' * 20, end = 
												txt_end)
		slp += 0.1
		if slp > 5.0: break
	
	try: # This step is to avoid error of overwriting files
		path = osjoin(folder, smarts_ext)
		df_sp = pd.read_csv(path, sep = ' ', index_col = 0)
		df_sp.index.name = 'wavelength'
		df_sp.rename(columns = {'Extraterrestrial_spectrm': 'ETR', 
				'Global_tilted_irradiance': 'smt2_irr'}, inplace = True)
			
	except PermissionError:
		if n == 9 : txt_end = '\n'
		print('\rPossible Error at: {}, due to PermissionError.'.format(
									timepoint), end = txt_end)
	
	except FileNotFoundError:
		if n == 9 : txt_end = '\n'
		print('\rPossible Error at: {}, due to FileNotFound.'.format(
									timepoint), end = txt_end)
		df_sp = None
	
	except OSError:
		if n == 9 : txt_end = '\n'
		print('\rPossible Error at: {}, due to OSError.'.format(
									timepoint), end = txt_end) 
	
	except ValueError:
		if n == 9 : txt_end = '\n'
		df_sp = None
	
	# If spectral_only is True, the function ends here just returning 
	# extra-terrestrial solar irradiance and global tilted irradiance.
	if spectral_only: return df_sp
	
	rtrv0 = lambda f, t, n : [f.index(t) + len(t), f.index(t) + len(t) + n]
	rtrv1 = lambda f, t, n : [f.index(t) + n, f.index(t)]

	fds = [['Pressure (mb) = ', 8, 'Pressure (mb)'], 
		['Ground Altitude (km) = ', 8, 'Altitude (m)'], 
		['Height above ground (km) = ', 8, 'Height (km)'], 
		['Relative Humidity (%) = ', 6, 'Rel_humidity (%)'], 
		['Precipitable Water (cm) = ', 7, 'Prec_water (cm)'], 
		['Ozone (atm-cm) = ', 6, 'Ozone (atm-cm)'], 
		['Optical Depth at 500 nm = ', 6, 'Optical Depth at 500 nm'], 
		['Optical depth at 550 nm = ', 6, 'Optical depth at 550 nm'], 
		["Angstrom's Beta = ", 6, "Angstrom's Beta"], 
		["Schuepp's B = ", 6, "Schuepp's B"], 
		['Meteorological Range (km) = ', 6, 'Meteorological Range (km)'], 
		['Visibility (km) = ', 6, 'Visibility (km)'], 
		['Alpha1 = ', 6, 'Alpha1'], ['Alpha2 = ', 6, 'Alpha2'], 
		["Mean Angstrom's Alpha = ", 6, "Mean Angstrom's Alpha"], 
		['Season = ', 8, 'Season'], 
		["Instantaneous at site's altitude = ", 5, 
									'Instantaneous_temp (K)'], 
		["Daily average (reference) at site's altitude = ", 5, 
										'Daily_av_temp (K)'], 
		['Stratospheric Ozone and NO2 (effective) = ', 5, 
									'Stratospheric_temp (K)'], 
		['Zenith Angle (apparent) = ', 6, 'Sol_Zth_Angle'], 
		['Azimuth (from North) = ', 7, 'Sol_Azh_Angle'], 
		['- Rayleigh = ', 6, 'OM_Rayleigh'], 
		['- Water Vapor = ', 6, 'OM_W_vapor'], 
		['- Ozone = ', 6, 'OM_Ozone'], ['- NO2 = ', 6, 'OM_NO2'], 
		['- Aerosols = ', 6, 'OM_Aerosols'], 
		['Year = ', 4, 'Year'], 
		['Month = ', 2, 'Month'], ['Day = ', 2, 'Day'], 
		['Hour (LST) = ', 6, 'Hour (LST)'], 
		['Day of Year = ', 3, 'Day_year'], 
		['Day (UT) = ', 2, 'Day (UT)'], 
		['Hour (UT) = ', 6, 'Hour (UT)'], 
		['Julian Day = ', 12, 'Julian_day'], 
		['Declination = ', 7, 'Declination'], 
		['Radius vector = ', 7, 'Radius_vec'], 
		['Equation of Time (min) = ', 7, 'EoT'], 
		['Local Apparent Time (or Solar Time): ', 7, 'Solar_time'], 
		['CO2 Mixing Ratio (ppmv): ', 6, 'CO2 (pmmv)'], 
		['Surface Tilt = ', 7, 'Surf_tilt'], 
		['Surface Azimuth (from North) = ', 7, 'Surf_Azh'], 
		['Incidence Angle = ', 7, 'Incidence_Angle'], 
		['(isotropic approximate conversion--for reference)', -8, 
										'Diff_irr_ratio_iso'], 
		['(anisotropic conversion model--used here)', -8, 
										'Diff_irr_ratio_ani']]
		
	try: # This step is to avoid error of overwriting files
		path = osjoin(folder, smarts_out)
		with open(path, 'r') as file:
			file = file.read()
			list_ = []
			for txt, N, head in fds[:-2]:
				i0, i1 = rtrv0(file, txt, N)
				try:
					list_ += [[head, float(file[i0: i1])]]
				except ValueError:
					list_ += [[head, file[i0: i1]]]
			for txt, N, head in fds[-2:]:
				i0, i1 = rtrv1(file, txt, N)
				list_ += [[head, float(file[i0: i1])]]
		list_ = dict(list_)
		df_sc = pd.DataFrame(list_, index = [timepoint])
		
	except PermissionError:
		if n == 9 : txt_end = '\n'
		print('\rPossible Error at: {}, due to PermissionError.'.format(
									timepoint), end = txt_end)
#		continue
	
	except FileNotFoundError:
		if n == 9 : txt_end = '\n'
		print('\rPossible Error at: {}, due to FileNotFound.'.format(
									timepoint), end = txt_end)
		df_sc = None
	
	except OSError:
		if n == 9 : txt_end = '\n'
		print('\rPossible Error at: {}, due to OSError.'.format(
									timepoint), end = txt_end) 
	
	except ValueError:
		if n == 9 : txt_end = '\n'
		df_sc = None

	chdir(c_fol) # restabilsh the previous working directory

	# if scalar_only is true it will come till this point and return scalar 
	# values as only
	if scalar_only: return df_sc



if __name__ == '__main__':

	# For demonstarion purposes it is only shown how this wrapping works for a 
	# single timepoint for a series of timepoints a loop running through a list 
	# of datetime data could be used.
	
	# If both spectral_only and scalar_only are False, the function returns both
	# outcomes.
	return df_sp, df_sc
	
		
	# Template file name
	t_file = 'SheffieldTemplate_INP_AllDays.txt'
	
	# Folder where template and batch executable files are found
	folder = r'.'
	
	lim_i = 300
	lim_s = 1050
	av_t = 10.0
	
	# Aerosol model is just an example, for more information see instruction 
	# manual of SMARTS2
	aer_m = "'B&D_C'"
	timepoint  = pd.to_datetime(15, 6, 2019, 12, 30, 30) 
	press = 1150
	air_t = 15
	rel_h = 85
	
	# smt2_sp is a CSV file containing the spectral outcomes, while smt2_sc is a 
	# CSV file with the scalar outcomes, pandas' function pd.to_csv() can be 
	# used to save the series of modelled spectra in a single CSV file
	smt2_sp, smt2_sc = SMARTS2(timepoint, press, air_t, rel_h, av_t, 
								aer_m, til_ang, azh_ang, lim_i, lim_s,	
								folder, t_file)


