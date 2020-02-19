# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 09:49:57 2016

@author: Roberto

This is used for:
- The final calibration of the spectra of 2047 points (only within the 450 to 750 nm range)
- Use SMARTS2 to extrapolate the spectra from 450-750 to 300-1050 nm
- Calculate total irradaiance
- Calculate the spectral photon flux to:
	- Calculate APE 
	- The Jsc of PV techs
-  USE SMARTS2 to calculate ETR to calculate Kt in base to the total irradiance

writen in Python 3.5.1

"""

import datetime as dt # builtin 
from os import remove, chdir, getcwd # builtin 
from os.path import join as osjoin # builtin 
from os.path import isfile # builtin 
from subprocess import call # builtin 
import time as T # builtin 
from functools import reduce # builtin 
from math import cos, sin, radians # builtin 

import pandas as pd # '0.18.0'
import numpy as np # '1.15.1'
from scipy.interpolate import interp1d # '0.17.0'


def Calibration(df, calcurve, noise):
	
	idx = pd.IndexSlice
	
	merge_on = 'wavelength'
	lim_a = 384.45
	lim_b = 750.0
	
	dfs = [df, calcurve, noise]
	df = reduce(lambda dfl, dfr: pd.merge(dfl.reset_index(), 
				dfr.reset_index(), on = merge_on).set_index(
								['datetime', merge_on]), dfs)
	df.sort_index(inplace = True)
	df.interpolate(inplace = True)
	df['cal_irr'] = (df['counts'] - df['noise']) * df['CF']
	df = df.loc[idx[:, lim_a:lim_b], :] # This range is safe of noise and 
								# IR interference
	#raise Exception
	df = df.drop(['CF', 'noise', 'counts'], axis = 1)
	
	timepoints = df.index.levels[0]
	
	return df, timepoints
	
		
def SMARTS2(timepoint, press, air_t, rel_h, av_t, aer_m, til_ang, 
		azh_ang, lim_i, lim_s, folder, t_file,
		spectral_only = False, scalar_only = False):
			
	smarts_ext = 'smarts295.ext.txt'
	smarts_out = 'smarts295.out.txt'
	c_fol = getcwd()
	chdir(folder)
	path = osjoin(folder, t_file)
	
	with open(path, 'r+') as param_template: # template input file
		# converts the file from bits to string
		content = param_template.read() 
		
		yr = timepoint.year; mt = timepoint.month; 
		dy = timepoint.day; hr = timepoint.hour; 
		mn = timepoint.minute ; sn = timepoint.second
		time =  hr + (mn / 60) + (sn / 3600)
		
		params = [['Press', press], ['Air_temp', air_t], ['Rel_hum', rel_h], 
			['Av_temp', av_t], ['Aer_mod', aer_m], ['Til_ang',  til_ang], 
			['Azh_ang', azh_ang], ['Lim_inf', lim_i], ['Lim_sup', lim_s], 
			['Year', yr], ['Month', mt], ['Day', dy], 
			['Time', time], ['Spc_rsl', '2'], ['Spc_out', '1 8']] # parameters to change
		
		for par, val in params: # changing the parameters in the template 
			content = content.replace(par, str(val))
		
		for n in range(10):
			try:
				with open('smarts295.inp.txt', 'w') as inputfile: # saving the edited template as inputfile
					#print(content)
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
		# Delete output files to prevent shifting data from timepoint
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
		# This run the executable
		batch = 'smarts295bat.exe'
		path = osjoin(folder, batch)
		call(path, shell = True) 
		T.sleep(slp) # I am not sure if this is necessary
		path = osjoin(folder, smarts_ext)
		ext_file = isfile(path) # To confirm the creation of the output file
		
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

	if scalar_only: return df_sc
	
	return df_sp, df_sc
	
		
def SMARTS2_extrapolation(df, smt2):

	sample = 100
	lim_a = 384.45
	lim_b = 750.0
	
	d = pd.to_datetime(df.index.levels[0], format = '%Y-%m-%d %H:%M:%S')[0]
	n_index = df.index.droplevel(level = 0)
	df.index = n_index
	
	B = pd.DataFrame(df.iloc[:sample])
	Bmin, Bmax = B.index[0], B.index[-1]
	
	R = pd.DataFrame(df.iloc[-sample:])
	Rmin, Rmax = R.index[0], R.index[-1]
	
	df_factor_B = B.join(smt2.loc[Bmin: Bmax, 'smt2_irr'], how = 'outer'
											).interpolate()
	df_factor_R = R.join(smt2.loc[Rmin: Rmax, 'smt2_irr'], how = 'outer'
											).interpolate()
	
	df_factor_B['factor'] = (df_factor_B[df_factor_B.columns[0]] / 
							df_factor_B[df_factor_B.columns[1]])
	df_factor_R['factor'] = (df_factor_R[df_factor_R.columns[0]] / 
							df_factor_R[df_factor_B.columns[1]])
	
	factor_B = df_factor_B['factor'].mean()
	factor_R = df_factor_R['factor'].mean()
	
	dfB = (smt2.loc[:lim_a, ['smt2_irr']] * factor_B)
	dfB.rename(columns = {'smt2_irr': 'cal_irr'}, inplace = True)
	dfR = (smt2.loc[lim_b:, ['smt2_irr']] * factor_R)
	dfR.rename(columns = {'smt2_irr': 'cal_irr'}, inplace = True)
	
	df = dfB.append(df).append(dfR)
	df['datetime'] = d
	df.set_index('datetime', append = True, inplace = True)
	df = df.reorder_levels([1, 0], axis = 0)
	
	return df

def APE(df, lim_i, lim_s, E = False, Ph = False, l_eff = False, 
										ph_flux = False):

	lim_i = float(lim_i)
	lim_s = float(lim_s) 
	df0 = df.copy(deep = True)
	n_index = df0.index.droplevel(level = 0)
	df0.index = n_index  
			
	q = 1.600217662E-19 #coulombs = (A * s)
	h = 6.62607004E-34 #m^2 * Kg * s^-1
	c = 299792458 #m * s^-1
	
	df1 = df0.copy(deep = True)
	df2 = df0.copy(deep = True)
	df1.rename(columns = {df1.columns[0]: 'ph_influx'}, inplace = True)
	df1['ph_influx'] = (df1['ph_influx'] * df1.index / 1E9) / (h * c)
	df2.rename(columns = {df2.columns[0]: 'energy'}, inplace = True)

	energy = np.trapz(df2.loc[lim_i:lim_s, 'energy'], df2.loc[
								lim_i:lim_s, :].index, axis = 0)
	if E: return energy
	
	photon_flux = np.trapz(df1.loc[lim_i:lim_s, 'ph_influx'], df1.loc[
								lim_i:lim_s, :].index, axis = 0)
	if Ph: return photon_flux
    
	ape = energy / (q * photon_flux)
	
	if l_eff:
		lambda_eff = (h * c) / (ape * 1E9 * q)
		return lambda_eff
	
	if ph_flux:
		return ape, df1.loc[:, ['ph_influx']]
    
	return ape
	
def SR2Jsc(df, PVtech, folder, file):
	# Originally this function would take the photon flux spectral curve,
	# multipy it times the electron charge (q) and integrate it.
	# This was a bad approach since Dirnberger data was not EQE but SR,
	# SR is already the (A/W) of irradiance
	
#	q = 1.600217662E-19 #coulombs = (A * s)
#	h = 6.62607004E-34 #m^2 * Kg * s^-1
#	c = 299792458 #m * s^-1
	
	path  = osjoin(folder, file) # path for SR file
	
	# SR spectrum
	SR = pd.read_csv(path, index_col = 0).sort_index().interpolate(
										method = 'index')
	SR.replace(to_replace = np.nan, value = 0, inplace = True)
	SR.sort_index(inplace =  True)
	
	SR_func = interp1d(SR[PVtech].index.values, SR[PVtech].values) 
	SR_df = df.copy(deep = True)
	SR_df[SR_df.columns[0]] = SR_df[SR_df.columns[0]].mul(SR_func(
							df.index), axis = 0, level = 0) 
	SR_Jsc = np.trapz(SR_df[SR_df.columns[0]], SR_df.index)
	
	return SR_Jsc
	
def Atm_functions(folder, file):
	
	path = osjoin(folder, file)	
	try:
		parameters_gdays = pd.read_csv(path, index_col = 0, 
										encoding = "utf-8")
	except UnicodeDecodeError: # To avoid problems with encoding
		parameters_gdays = pd.read_csv(path, index_col = 0, 
									encoding = 'ISO-8859-1')
										
	datetime_to_float = parameters_gdays.index.astype('datetime64'
										).values.astype(float)
	
	func_press = interp1d(datetime_to_float, 
							parameters_gdays['pressure (mbar)'])
	func_airtemp = interp1d(datetime_to_float, 
							parameters_gdays[r'temperature (°C)'])
	func_humidity = interp1d(datetime_to_float, 
								parameters_gdays['humidity (%)'])
	
	return func_press, func_airtemp, func_humidity
	

def Deltax(folder, file):
	
	path = osjoin(folder, file)
	Delta_x = pd.read_csv(path, index_col = 2, header = None, 
										nrows = 2047).index
	
	return Delta_x

def CalCurve(folder, file, Delta_x):
	
	path = osjoin(folder, file)
	C_factor = pd.read_csv(path, index_col = 0).dropna()
	C_factor.rename(columns = {'calcurve': 'CF'}, inplace = True)
	CalCurv = interp1d(C_factor.index, C_factor['CF'], bounds_error = None, 
									fill_value = "extrapolate")
	calcurve = pd.DataFrame(index = Delta_x, data = CalCurv(Delta_x), 
										columns = ['CF'])
	calcurve.index.name = 'wavelength'
	
	return calcurve

def Noise(folder, file, sensor, Delta_x):
	
	sensor = str(sensor)
	path = osjoin(folder, file)
	noise_mean = pd.read_csv(path, index_col = 0, usecols = [0, sensor]
				).dropna() # change from 4 to 5 according to Sensor
	noise_mean.rename(columns = {sensor : 'noise'}, inplace  = True)
	noise = interp1d(noise_mean.index, noise_mean['noise'], 
					bounds_error = None, fill_value = "extrapolate")
	noise = pd.DataFrame(index = Delta_x, data = noise(Delta_x), 
										columns = ['noise'])
	noise.index.name = 'wavelength'
	
	return noise

def SensorVals(sensor):
	
	raw_files = {'4' : 'RAW2047_s4.csv', '5' : 'RAW2047_s5.csv'}
	cal_files = {'4' : 'SMARTS2_CalCurve_S4_03.csv', 
								'5' : 'SMARTS2_CalCurve_S5.csv'}
	dx_files = {'4' : 'RAW2047_s4.csv', '5' : 'RAW2047_s5.csv'}
	til_angs = {'4' : 45.0, '5' : 0.0}
	azh_angs = {'4' : 180.0, '5' : 0.0}
	tps = {'4' : 156555, '5' : 152210}
	
	return raw_files[sensor], cal_files[sensor], dx_files[sensor], \
	til_angs[sensor], azh_angs[sensor], tps[sensor]
	
def AM(df, timepoint):
	
	# From Kasten, F.; Young, A. T. (1989). "Revised optical air mass tables
	# and approximation formula"
	z = df.at[timepoint, 'Sol_Zth_Angle']
	am = 1 / (cos(radians(z)) + 0.50572 * (96.07995 - z) ** -1.6364)
	
	return am

def CosAOI(df, timepoint):
	
	# From Sandia National Laboratories (https://pvpmc.sandia.gov)
	list_lbs = ['Sol_Zth_Angle', 'Sol_Azh_Angle', 'Surf_tilt', 'Surf_Azh']	
	angs = []
	
	for lb in list_lbs:
		angs += [df.at[timepoint, lb]]
		
	S_zth, S_azh, A_til, A_azh = list(map(radians, angs))
	
	cos_aoi = cos(S_zth) * cos(A_til) + \
	sin(S_zth) * sin(A_til) * cos(S_azh - A_azh)
	
	return cos_aoi

def FinalAnalyses():
	
	t_file = 'SheffieldTemplate_INP_AllDays.txt'
	ns_file = 'darknoise_mean.csv'
	atm_file = 'AllDays_AtmData.csv'
	sr_file = 'Dirnberger_SpectralResponse.csv'
	
	# for all files but darknoise
	folder = \
	r'C:\Users\Roberto\Documents\PhD 2013\Spectral\SMARTS_295_PC\SMARTS_295_PC'
	
	# for darknoise and SR
	n_folder = r'C:\Users\Roberto\Documents\PhD 2013\Experiments\Roof'
	
	lim_i = 300
	lim_s = 1050
	av_t = 10.0
	aer_m = "'B&D_C'"
	PVtechs = ['a-Si', 'c-Si', 'CdTe', 'CIGS', 'HE c-Si']
	sensors = ['4', '5']
	s0 = sensors[0]
	
	
	func_press, func_airtemp, func_humidity = Atm_functions(folder, atm_file)
	
	for sensor in sensors[1:]:
		
		start_time = dt.datetime.now()
		
		file, cal_file, dx_file, til_ang, azh_ang, tps = SensorVals(sensor)
	
		deltax = Deltax(folder, dx_file)
		calcurve = CalCurve(folder, cal_file, deltax)	
		
		noise = Noise(n_folder, ns_file, sensor, deltax)
		
		cols = [0, 2, 3] 
		indices = [0, 1]
		names = ['datetime', 'wavelength', 'counts']
		chunk = 2047
		
		new = True
	
		path = osjoin(folder, file)
		reader =  pd.read_csv(path, chunksize = chunk, 
					index_col = indices, header = None, 
					usecols  = cols, names = names, parse_dates = [0])
	
		for n, df in enumerate(reader, 1):
						
			df_cal, timepoints = Calibration(df, calcurve, noise)
			
			if len(timepoints) == 1:
				
				timepoint  = timepoints[0]
				press = func_press(timepoint.value)
				air_t = func_airtemp(timepoint.value)
				rel_h = func_humidity(timepoint.value)
				
				
				smt2_sp, smt2_sc = SMARTS2(timepoint, press, air_t, 
								rel_h, av_t, aer_m, til_ang,
								azh_ang, lim_i, lim_s,	folder, 
												t_file)
				if smt2_sp is None: continue

				df_cal = SMARTS2_extrapolation(df_cal, smt2_sp)
				ape, df_ph = APE(df_cal, lim_i, lim_s, ph_flux = True)
				func_etr = interp1d(smt2_sp.index, smt2_sp['ETR'])
				func_ITot = interp1d(smt2_sp.index, smt2_sp['smt2_irr'])
				func_PhIfx = interp1d(df_ph.index, df_ph['ph_influx'])

				df_cal = pd.merge(df_cal.reset_index(), 
								df_ph.reset_index(), 
								on = 'wavelength', how = 'outer')
				
				df_cal['ETR'] = func_etr(df_cal['wavelength'])
				df_cal['smt2_irr'] = func_ITot(df_cal['wavelength'])
				df_cal['ph_influx'] = func_PhIfx(df_cal['wavelength'])
				df_cal.set_index(['datetime', 'wavelength'], 
											inplace = True)

				am = AM(smt2_sc, timepoint)
				cosaoi = CosAOI(smt2_sc, timepoint)
				
				values = []
				to_integrate = ['ETR', 'cal_irr', 'smt2_irr', 
											'ph_influx']
				for col in to_integrate:			
					values += [np.trapz(df_cal[col].values, 
									df_cal.index.levels[1])]
				etr, cal_irr, smt2_irr, ph_influx = values
				
				etr_srf = etr * cosaoi
				kt = cal_irr / etr_srf
				
				scal_labels = to_integrate + ['cosAOI', 'ETR_srf', 'Kt',
											'AM', 'APE']
				scal_values = values + [cosaoi, etr_srf, kt, am, ape]
				
				for val, lb in zip(scal_values, scal_labels):
					smt2_sc.loc[timepoint, lb] = val
				
				for PVtech in PVtechs:
					df_irr = df_cal.loc[:, ['cal_irr']]
					Jsc = SR2Jsc(df_irr, PVtech, n_folder, sr_file)
					smt2_sc.at[timepoint, PVtech] = Jsc
				
				if n > 20: return df_cal, smt2_sc; raise Exception
				raise Exception
				if new:
					md = 'w'
					new = False
					header = True
				
				else:
					md = 'a'
					header = False
				if s0 != sensor: end_txt = '\n'
				else: end_txt = '\r'
					
				df_cal.to_csv('CalSpectra_S{}_Final.csv'.format(sensor), 
									mode = md, header = header)
				smt2_sc.to_csv('CalScalars_S{}_Final.csv'.format(sensor), 
									mode = md, header = header)
				current_time = dt.datetime.now()
				elapsed_time = current_time - start_time
				hours, remdr = divmod(elapsed_time.seconds, 3600)
				mins, secs = divmod(remdr, 60)
				pct = (n / tps) * 100
				tx0 = '\rSensor:{}, date: {}, '.format(sensor, 
											timepoint)
				tx1 = 'time elapsed: {:.0f}:{:.0f}:{:.0f}, '.format(
										hours, mins, secs)
				tx2 = '{:.2f}% of timepoints processed'.format(pct)
				print('\r{}{}{}'.format(tx0, tx1, tx2), 
									' ' * 20,  end = end_txt)
				s0 = sensor

if __name__ == '__main__':
#	df_cal, smt2_sc = FinalAnalyses()
	FinalAnalyses()