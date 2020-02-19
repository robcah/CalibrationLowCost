# -*- coding: utf-8 -*-
"""

Created on Sun Apr 12 23:53:05 2015
@author: R Cahuantzi
In languague: Python V3.5.1, with built-in modules math and os.
With modules: numpy V1.10.4, pandas V0.18.0.
Files required: ETR_Gueymard2003.csv, A_H2O_Bird&Riordan1983.csv, 
A_O3_Bird&Riordan1983.csv, A_UG_Bird&Riordan1983.csv

Pysolar V0.6 geometry based on: I. Reda and A. Andreas, “Solar Position Algorithm 
for Solar Radiation Applications,” National Renewable Energy Laboratory, 
NREL/TP-560-34302, revised November 2005.

"""

from math import acos, cos, sin, pi, log, exp, floor
from math import degrees as deg
from math import radians as rad

import numpy as np
import pandas as pd
from datetime import timedelta
import os

from Pysolar.solar import * 
from Pysolar import julian


"""
Code based on SPECTRAL2 irradiance model. It is able to calculate solar 
position and spectral irradiance, some modifications were made, specially 
the inclusion of solar geometry calculations. The equations used in this 
code come from:
[1] R. Bird and C. Riordan, "Simple solar spectral model for direct and 
diffuse irradiance on horizontal and tilted planes at the Earth's surface 
for Cloudless Atmospheres", 1984.
[2] T. K. Van Heuklon, "Estimating Atmospheric Ozone for Solar Radiation 
Models", 1978.
[3] J.W. Spencer, "Fourier Series Representations of the Position of the Sun", 
1971. From J. Pickard's email (corrected by M. Oglesby), 1998.
[4] B. Leckner, "The Spectral Distribution of Solar Radiation at the Earth's 
Surface--Elements of a Model", 1978.
[5] F. X. Kneizys, et al., "Atmospheric Transmittance/Radiance: 
Computer Code LOWTRAN5".
[6] A. Angstrom, "Technique of Determining Turbidity of the Atmosphere", 1961.
[7] A. McEvoy, T. Markvart and Luis Castaner, "Practical Handbook of 
Photovoltaics: Fundamental and applications", 2013.
[8] J. A. Duffie and W.A. Beckman, "Solar Engineering of Thermal Processes", 
2013.
[9] M. Jacobson, "Fundamentals of atmospheric modelling", 2005.
[10] C. Gueymard, "SMARTS2, A simple model of the atmosphere radiative 
transfer of sunshine", 1995.
"""

# Defining default values
# For time:
year = 2015.0;		month = 6.0;		day = 21.0;		
hour = 12.0;		mins = 0.0;			secs = 0.0; 
consider_bst = False;

# For PV/acquisition system:
Ang_slope_sys = 45.0;		Ang_tilt_sys = 180.0;		tracking = False;

# Geographycal location:
lat = 53.380813; 			long = -1.485708;			elev = 0

# Atmospheric parameters:
AOD = 0.27; 		alpha = 1.14;		rho_g = 0.10;
O3 = 0.3438; 		H2O = 1.4164;		pr = 1013.25;
omega = 0.945; 		omega_p = 0.095;		asym = 0.65; 
O3_h = 22.;			t = 25

alt_return = ''
	
def Spctral(year = year, 	month = month, 	day = day, 		hour = hour, 
		mins = mins, 	secs = secs, 	consider_bst = consider_bst, 
		
		Ang_slope_sys = Ang_slope_sys, 	Ang_tilt_sys = Ang_tilt_sys,
		tracking = tracking, 

		lat = lat, 		long = long, 	elevation = elev,	
		
		AOD = AOD,		alpha = alpha,	rho_g = rho_g,	O3 = O3, 	
		H2O = H2O, 		pr = pr, 		omega = omega, 	
		omega_p = omega_p, asym = asym, 	O3_h = O3_h,
		
		alt_return = alt_return):

	"""
	Inputs:
	consider_bst : Activates the option to recognise the dates within British 
		Summer Time and make the arithmetic correction of -60 minutes.
	Ang_slope_sys : Angle of inclination of PV system's plane subtended to 
		ground surface, in degrees.
	Ang_tilt_sys : Angle of inclination of the projection of the normal of 
		the plane of the system on PV system's plane of the surface 
		subtended to the north axis. nown as well as azimuth angle 
		(180 south, 0 north, 90 east, 270 west).
	tracking : Activates the calculation as if the PV system has a tracking 
		base.
	lat : Geographical latitud of PV system, in degrees.
	long : Geographical longitud of PV system, in degrees.
	AOD : Aerosol Optical Dept: Range for clear skies is 0.05 to 0.55. 
		Higher optical depths result from clouds, smoke, and larger 
		particles in the atmosphere. 
	alpha : Exponent alpha of Angstrom's [6] expression of turbidity.
	rho_g : Albedo, surface reflectance used to calculate the diffuse 
		irradiance, does not affect Direct Beam computation. Used for 
		computing diffuse sky and reflected diffuse irradiance.
	O3 : The total column of ozone, as if it were condensed on the surface, 
		in cm.
	H2O : Total precipitable water if it were condensed on the surface, 
		in cm.
	pr : Atmospheric pressure. The standard at sea level is 1013.25 mbar. 
		Lower pressures mean less atmosphere to absorb the radiation, in 
		millibars (mbar).
	omega : Single scattering albedo factor at 0.4 micrometers (400 nm) of 
		wavelength.
	omega_p : Prime single scattering albedo factor for wavelength variation.
	asym : Rural aerosol scattering asymmetry factor, forward to total 
		backwards scattering ratio.
	O3_h : Height of ozone, in Km.
	
	Output:
	DataFrame with spectra as columns (units W/m^2/nm):
	ETR : Extraterrestrial spectral irradiance.
	ETR_rv_tilt :  Extraterrestrial spectral irradiance on PV system's plane.
	I_dir :  Direct solar irradiance on horizontal.
	I_dir_tilt : Direct solar irradiance on PV system's plane.
	I_dif_tilt : Diffuse solar irradiance on PV system's plane.
	I_tot_tilt : Total solar irradiance on PV system's plane.
	"""
	
	Eq_T, Ang_hr, Ang_dcl, AM, \
	AM_pr, r_vec, Ang_sol_az, Ang_slp, \
	Ang_zth, cos_AOI_sys = SunPos(year = year, 	month = month, 	day = day,
						hour = hour, 	mins = mins, 	secs = secs, 	
						consider_bst = consider_bst, 
		
						Ang_slope_sys = Ang_slope_sys, 	
						Ang_tilt_sys = Ang_tilt_sys,
						tracking = tracking, 

						lat = lat, 		long = long,
						elevation =  elev)
	
	list_alt_ret = [Eq_T, Ang_hr, Ang_dcl, AM, AM_pr, r_vec, Ang_sol_az, \
			Ang_slp, Ang_zth, cos_AOI_sys]
	
	list_str = []
	for var in list_alt_ret:
		for k, v in list(locals().items()):
			if v is var and k != 'var':
				list_str += [[k, v]]
	
	for a_ret in list_str:
		if alt_return == a_ret[0]: return a_ret[1]
	
	O3_m = O3AM(O3_h, Ang_zth)
	alog = Alog(asym)
	aF_s = AFs(alog)
	bF_s = BFs(alog)
	F_s = Fs(aF_s, bF_s, Ang_zth)
	F_s_p = Fsp(aF_s, bF_s)
	
	ETR, A_H2O, A_O3, A_UG, WL_range = EetrAbs()
	omega_wl = OmegaWl(WL_range, omega, omega_p)
	T_R = Tr(WL_range, AM_pr)
	T_O3 = To3(WL_range, A_O3, O3_m, O3)
	T_UG = Tug(WL_range, A_UG, AM_pr)
	T_H2O = Th2o(WL_range, A_H2O, H2O, AM)
	tau = Tau(WL_range, AOD, alpha)
	T_a = Ta(WL_range, tau, AM)
	T_as = Tas(omega_wl, tau, AM)
	T_aA = Taa(omega_wl, tau, AM)
	I_dif_R = Edifr(ETR, T_O3, T_H2O, T_UG, T_aA, T_R, r_vec, Ang_zth)
	I_dif_a = Edifa(ETR, T_O3, T_H2O, T_UG, T_aA, T_R, T_as, Ang_zth, 
											F_s, r_vec)
	T_R_p = Trp(WL_range)
	T_H2O_p = Th2op(A_H2O, WL_range, H2O)
	T_UG_p = Tugp(A_UG, WL_range)
	T_as_p, T_aA_p = TaspTaap(omega_wl, tau)
	rho_s = Rhos(T_UG_p, T_H2O_p, T_aA_p, T_R_p, T_as_p, F_s_p)
	I_dir = Edir(ETR, T_R, T_O3, T_UG, T_H2O, T_a, r_vec)
	I_dif_g = Edifg(I_dir, I_dif_R, I_dif_a, rho_s, rho_g, Ang_zth)
	global C_s # Declare C_s as a global variable
	C_s = Cs(WL_range)
	I_dif_h = EdifhWL(I_dif_R, I_dif_a, I_dif_g, C_s)
	I_dif_c_g = EdifcgWL(I_dir, I_dif_h, Ang_zth, rho_g, Ang_slp)
	I_dif_c_sc = EdifcscWL(I_dif_h, I_dir, ETR, r_vec, cos_AOI_sys, Ang_zth)
	I_dif_c_ssi = EdifcssiWL(I_dif_h, I_dir, ETR, r_vec, Ang_slp)
	I_dif_srf = EdifslpWL(I_dif_c_g, I_dif_c_sc, I_dif_c_ssi)
	I_dir_srf = EdirslpWL(I_dir, cos_AOI_sys)
	ETR_rv_srf = EetrrvslpWL(ETR, r_vec, cos_AOI_sys)
	I_tot_srf = Etotslp(I_dir_srf, I_dif_srf)
	ModIrr = Etot(I_dir, ETR, I_dif_srf, ETR_rv_srf, I_tot_srf, I_dir_srf)

	return ModIrr
	
def SunPos(year = year, 	month = month, 	day = day, 		hour = hour, 
		mins = mins, 	secs = secs, 	consider_bst = consider_bst, 
		
		Ang_slope_sys = Ang_slope_sys, 	Ang_tilt_sys = Ang_tilt_sys,
		tracking = tracking, 

		lat = lat, 		long = long,	elevation =  0):
	
	"""
	Inputs:
	consider_bst : Activates the option to recognise the dates within British 
	Summer Time and make the arithmetic correction of -60 minutes.
	Ang_slope_sys : Angle of inclination of PV system's plane subtended to 
		ground surface, in degrees.
	Ang_tilt_sys : Angle of inclination of the projection of the normal of 
		the plane of the system on PV system's plane of the surface 
		subtended to the north axis. nown as well as azimuth angle 
		(180 south, 0 north, 90 east, 270 west).
	tracking : Activates the calculation as if the PV system has a tracking 
		base.
	lat : Geographical latitud of PV system, in degrees.
	long : Geographical longitud of PV system, in degrees.
	
	Output:
	AM : Air mass
	AM_pr : Pressured corrected air mass 
	r_vec : Vector correction of the Earth-Sun distance
	Ang_slp : Slope angle of the system, considered the option of tracking
	Ang_zth :  Zenith angle of the sun, called altitude as well
	cos_AOI_sys: Cosine of the angle of incidence of the system
	"""
	
	date_time = Datetime(year, month, day, hour, mins, secs)
	bst = BST(date_time, consider_bst)
	datetime_corr = date_time if not bst else date_time - \
										timedelta(hours = 1)
								 
	Ang_dy = AngDay(datetime_corr) 
	Eq_T = EquationOfTime(GetDayOfYear(datetime_corr))
	Ang_hr = -GetHourAngle(datetime_corr, long)
	Ang_dcl = GetDeclination(GetDayOfYear(datetime_corr))
	Ang_zth = 90 - GetAltitude(lat, long, datetime_corr)
	AM = GeoAM(Ang_zth)
	AM_pr = AMp(AM, pr)
	r_vec = Crv(Ang_dy)
	az = GetAzimuth(lat, long, datetime_corr)
	Ang_sol_az = az if az > -180 else az + 360
	Ang_slp  = AngSlope(Ang_slope_sys, Ang_zth, tracking)
	Ang_tilt_sys = AngTilt(Ang_sol_az, Ang_tilt_sys, tracking)
	
	"""
	Incidence angle calculation based on pysolar
	"""
	jd = julian.GetJulianDay(datetime_corr)
	jde = julian.GetJulianEphemerisDay(jd, 65)
	jce = julian.GetJulianEphemerisCentury(jde)
	jme = julian.GetJulianEphemerisMillenium(jce)
	geocentric_longitude = GetGeocentricLongitude(jme)
	nutation = GetNutation(jde)
	radius_vector = GetRadiusVector(jme)
	aberration_correction = GetAberrationCorrection(radius_vector)
	apparent_sun_longitude = GetApparentSunLongitude(geocentric_longitude, 
								nutation, aberration_correction)
	true_ecliptic_obliquity = GetTrueEclipticObliquity(jme, nutation)
	geocentric_latitude = GetGeocentricLatitude(jme)
	geocentric_sun_declination = GetGeocentricSunDeclination(
				apparent_sun_longitude, true_ecliptic_obliquity, 
										geocentric_latitude)
	latitude_deg = latitude = lat
	projected_axial_distance = GetProjectedAxialDistance(elevation,																								latitude_deg)
	equatorial_horizontal_parallax = GetEquatorialHorizontalParallax(
											radius_vector)
	projected_radial_distance = GetProjectedRadialDistance(elevation, 
											latitude_deg)
	equatorial_horizontal_parallax = GetEquatorialHorizontalParallax(
											radius_vector)
	apparent_sidereal_time = GetApparentSiderealTime(jd, jme, nutation)
	longitude_deg = long
	geocentric_sun_right_ascension = GetGeocentricSunRightAscension(
					apparent_sun_longitude, true_ecliptic_obliquity, 
										geocentric_latitude)
	local_hour_angle = GetLocalHourAngle(apparent_sidereal_time, 
					longitude_deg, geocentric_sun_right_ascension)
	parallax_sun_right_ascension = GetParallaxSunRightAscension(
			projected_radial_distance, equatorial_horizontal_parallax, 
						local_hour_angle, geocentric_sun_declination)
	topocentric_sun_declination = GetTopocentricSunDeclination(
				geocentric_sun_declination, projected_axial_distance, 
				equatorial_horizontal_parallax, 
					parallax_sun_right_ascension, local_hour_angle)
	topocentric_local_hour_angle = GetTopocentricLocalHourAngle(
					local_hour_angle, parallax_sun_right_ascension)
	pressure_millibars = pr
	temperature_celsius = t
	topocentric_zenith_angle = GetTopocentricZenithAngle(latitude, 
							topocentric_sun_declination,
								topocentric_local_hour_angle, 
										pressure_millibars, 
										temperature_celsius)
	slope = Ang_slope_sys
	# This seems to be necesary to adapt to spctral method
	slope_orientation = Ang_tilt_sys 
	topocentric_azimuth_angle = GetTopocentricAzimuthAngle(
		topocentric_local_hour_angle, latitude, topocentric_sun_declination)
	cos_AOI_sys = cos(rad(GetIncidenceAngle(topocentric_zenith_angle, 
								slope, slope_orientation, 
									topocentric_azimuth_angle)))
		
	return Eq_T, Ang_hr, Ang_dcl, AM, AM_pr, r_vec, Ang_sol_az, Ang_slp, \
										Ang_zth, cos_AOI_sys

def Datetime(year, month, day, hour, mins, secs):
	
	dt = [year, month, day, hour, mins, secs]
	dt = list(map(int, dt))
	
	return pd.datetime(*dt)
	
	
def BST(date_time, bst):
	
	if bst:
	
		# British Summer Time Limited up to year 2021.
		BST_lim = [[[2014, 3, 30, 1], [2014, 10, 26, 2]],
			[[2015, 3, 29, 1], [2015, 10, 25, 2]],
			[[2016, 3, 27, 1], [2016, 10, 30, 2]],
			[[2017, 3, 26, 1], [2017, 10, 29, 2]],
			[[2018, 3, 25, 1], [2018, 10, 28, 2]],
			[[2019, 3, 31, 1], [2019, 10, 27, 2]],
			[[2020, 3, 29, 1], [2020, 10, 25, 2]],
			[[2021, 3, 28, 1], [2021, 10, 31, 2]]]
			
		#for i in range(len(BST_lim[0])):
		for dt in BST_lim:
			if pd.datetime(*dt[0]) < date_time < pd.datetime(*dt[1]):
				bst = True
				break
			else:
				bst = False
					
	return bst
		
def AngDay(date_time):	
	
	'''
	Number of the day, some sources call it julian number day
	From J.A. Duffie and W.A. Beckman [8]. 
	Yield result in radians but converted to degrees at 
	the end of the function.
	'''
	d = pd.Series(date_time).dt.dayofyear[0] # to get the day of year
	year = date_time.year
	
	# Formula to identify leap years
	year_days = 366 if ((year % 4 == 0 and year % 100 != 0) or 
					(year % 400 == 0 and year % 100 != 0)) else 365
	
	Ang_dy  = 2 * pi * ((d - 1) / year_days)
	
	return deg(Ang_dy)

def EoT(Ang_dy):	
	
	'''	
	Equation of Time in minutes = 
				True Solar time - Local Standard time - Long correction.
	From: J.W. Spencer [3].
	Yields the result in MINUTES, according to J.A. Duffie and 
	W.A. Beckman [8].
	Multiplier 229.18 converts the result to minutes, without this 
	the equation would yield in radians (4 minutes = 1 degree).
	'''
	Eq_T = (0.0000075 + 0.001868 * cos(rad(Ang_dy)) - 0.032077 \
		* sin(rad(Ang_dy)) - 0.014615 * cos(2 * rad(Ang_dy)) - 0.040849 \
							* sin(2 * rad(Ang_dy))) * (229.18)
	
	return Eq_T
	
def AngHour(hour, mins, secs, Eq_T, long, bst):

	'''
	From Solar engineering of thermal processes, J.A. Duffie and 
	W.A. Beckman (2013) [8].
	Solar time - standard time = (Long_ST - Long) * 4 + Eq_T
	Ang_hr = (15 * Solar time) - 180
	4 min = 1 deg, 1 hour = 15 deg
	Yields the results in degrees.
	'''
	
	if bst:
		mins = mins - 60
		
	Ang_hr = (15 * ((hour + (mins / 60) + (secs / (60 * 60))) \
				+ ((floor(long) - long) * 4) / 60 + Eq_T / 60)) - 180
	
	return Ang_hr
	
def AngDcl(Ang_dy): 
	
	'''	
	Declination angle (degrees), above or below ecliptic [3].
	Yields result in DEGREES
	'''
	Ang_dcl = ((0.006918 - 0.399912 * cos(rad(Ang_dy)) + 0.070257 \
		* sin(rad(Ang_dy)) - 0.006758 * cos(2 * rad(Ang_dy)) + 0.000907 \
			* sin(2 * rad(Ang_dy)) - 0.002697 * cos(3 * rad(Ang_dy)) \
					+ 0.00148 * sin(3 * rad(Ang_dy))) * (180 / pi))
		
	return Ang_dcl
		
def AngZth(Ang_dcl, lat, Ang_hr):
	
	'''
	Angle of Solar Zenith, it is the complement of Solar elevation 
	(or altitude angle, e), From D L Hartman, 1994.
	Z= 90 - e, cos(Z)=cos(d)cos(L)cos(H)+sin(d)sin(L) 
	(*this last one should be L), 
	e= solar elevation, L= Latitude, d= Declination, H= hour angle
	For solar altitude the arcsin would be required in stead of arccosine. 
	Yields result in RADIANS converted to degrees.
	'''
	
	Ang_zth = (acos(cos(rad(Ang_dcl)) * cos(rad(lat)) * cos(rad(Ang_hr)) 
							+ sin(rad(Ang_dcl)) * sin(rad(lat))))
	
	return deg(Ang_zth)
		
def GeoAM(Ang_zth): 
	
	'''
	Geometrical Air Mass (path length), called "M" in Bird's paper [1]. 
	Eq. 2-5, pg 4
	'''
	if Ang_zth > 91.013177 :
		AM = np.inf
	else:
		AM = ((cos(rad(Ang_zth))) + 0.15 * (93.885 - Ang_zth) ** (-1.253))\
													** -1
	if type(AM) is complex : AM = np.inf
	
	return AM
		
def AMp(AM, pr):
	
	'''
	Pressure corrected Air Mass. AM_pr being = AM * P / P0. Seen in 
	Bird's paper [1]. Eq. 2-5
	P0 = 1013 milibars
	'''
	AM_pr = AM * (pr/1013)
	
	return AM_pr
	
def O3AM(O3_h, Ang_zth):
	
	'''
	Effective O3 Air mass, h_o = 22 Km which is the height of maximum ozone 
	concentration [1], Eq 2-10 pg 5.
	'''
	O3_m = ((1 + O3_h / 6370) / (cos(rad(Ang_zth)) ** 2 + 2 * (O3_h / 6370))
												** 0.5)
	
	return O3_m
	
def Crv(Ang_dy):
	
	'''
	Radius vector correction for Earth-Sun distance, called "D" in 
	Bird's paper [1]. Eq 2-2 pg 2, also appears in M. Jacobson [9].
	It is the calculation (1 / r ** 2) allows the correction of the solar 
	constant to be applied in the given date.
	Yields the result in Astronomical Units (AU),
	AU = 149,597,870,700 m (the mean distance from the sun to earth).
	'''
	r_vec = (1.00011 + 0.034221 * cos(rad(Ang_dy)) + 0.00128 
			* sin(rad(Ang_dy)) + 0.000719 * cos(2 * rad(Ang_dy)) 
							+ 0.000077 * sin(2 * rad(Ang_dy)))
	
	return r_vec
	
def AngSolAz(Ang_dcl, Ang_hr, Ang_zth, lat):
	
	'''
	Approximate Solar Azimuth angle (0=N, 90=E, 180=S, 270=W) computed 
	from 180 + HA
	From Solar engineering of thermal processes, J.A. Duffie and 
	W.A. Beckman [8].
	Yield result in radians to be converted at the end into degrees.
	'''
	Ang_sol_az = (np.sign(Ang_hr) * abs(acos((cos(rad(Ang_zth)) 
			* sin(rad(lat)) - sin(rad(Ang_dcl))) / (sin(rad(Ang_zth)) 
										* cos(rad(lat))))))
		
	return deg(Ang_sol_az)
		
def AngSlope(Ang_slope_sys, Ang_zth, tracking):
	
	'''
	This part calculates the slope angle of the system. If 'tracking' is True
	it is considered that the slope and tilt of the system are always facing 
	directly to the sun, therefore the slope angle is equal to the angle 
	of solar zenith and tilt equals to solar azimuth.
	'''
	if tracking:
		Ang_slope_sys = Ang_zth
	
	return Ang_slope_sys
	
def AngTilt(Ang_sol_az, Ang_tilt_sys, tracking):
	
	'''
	This part intent to calculate the angle of surface azimuth in degrees.
	If tracking == True it is considered to be in a tracking module 
	the azimuth of the system is always facing directly to the sun.
	Also translate the value of south from 180 deg to 0 deg, with negative 
	values for east facing angles.
	'''
	if tracking:
		Ang_tilt_sys = Ang_sol_az 
	else:
		Ang_tilt_sys = (Ang_tilt_sys - 180) 
	
	return Ang_tilt_sys

def CosAOI(Ang_hr, Ang_dcl, lat, Ang_tilt_sys, Ang_slope_sys = 0):
		
	'''
	This function to limits the ETR on tilt to the time range of ETR on 
	horizontal. Cosine of incidence angle based on Ang_slope_sys of surface and 
	solar geometry. 
	From: J. A. Duffie and W. A. Beckman [8] and A. McEvoy, T. Markvart 
	and L. Castaner [7] formulas in p. 626 and p.624: 13, 9, 10 and 12.
	'''
	cos_AOI = (sin(rad(Ang_dcl)) * sin(rad(lat)) * cos(rad(Ang_slope_sys)) 
			- sin(rad(Ang_dcl)) * cos(rad(lat)) * sin(rad(Ang_slope_sys)) 
			* cos(rad(Ang_tilt_sys)) + cos(rad(Ang_dcl)) * cos(rad(lat)) 
			* cos(rad(Ang_slope_sys)) *cos(rad(Ang_hr)) + cos(rad(Ang_dcl)) 
			* sin(rad(lat)) * sin(rad(Ang_slope_sys)) * cos(rad(Ang_tilt_sys)) 
			* cos(rad(Ang_hr)) + cos(rad(Ang_dcl)) * sin(rad(Ang_slope_sys)) 
						* sin(rad(Ang_tilt_sys)) * sin(rad(Ang_hr)))
	
	return cos_AOI
	
def CosAOI_srf(Ang_hr, Ang_dcl, lat, Ang_tilt_sys, tracking, Ang_slope_sys): 
	
	'''
	Cosine of incidence angle based on Ang_slope_sys of surface and solar 
	geometry. 
	From: Solar engineering of thermal processes, J.A. Duffie and W.A. 
	Beckman [8].
	'''
	if CosAOI(Ang_hr, Ang_dcl, lat, Ang_tilt_sys) <= 0:
		cos_AOI_sys = 0
	else:
		if not tracking:
			cos_AOI_sys = CosAOI(Ang_hr, Ang_dcl, lat, Ang_tilt_sys, 
												Ang_slope_sys)
		else:
			cos_AOI_sys = 1
	
	return cos_AOI_sys
	
def Alog(asym):
	
	alog = log(1 - asym)
	
	return alog
		
def AFs(alog):
	
	aF_s = alog * (1.459 + alog * (0.1595 + alog * 0.4129))
	
	return aF_s
		
def BFs(alog):
	
	bF_s = alog * (0.0783 + alog * (-0.3824 - alog * 0.5874))
	
	return bF_s
		
def Fs(aF_s, bF_s, Ang_zth): 
	
	'''
	Fraction of the forward aerosol scatter ratio dependant on the solar 
	zenith angle. [1] Eq: 3-11, pg 7.
	'''
	F_s = 1 - 0.5 * exp((aF_s + bF_s * cos(rad(Ang_zth))) * 
										cos(rad(Ang_zth)))
	
	return F_s
		
def Fsp(aF_s, bF_s):
	
	'''
	Prime fraction of the downward aerosol scattering ratio dependant on the 
	solar zenith angle,  (evaluated at Air Mass = 1.8). [1] Eq: 3-15, pg 7.
	'''
	F_s_p = 1 - 0.5 * exp((aF_s + bF_s / 1.8) / 1.8)
	
	return F_s_p
	
def EetrAbs(ETR = 'ETR_Gueymard2003_ndelta.csv', 
			A_H2O = 'A_H2O_Bird&Riordan1983_ndelta.csv',
				A_O3 = 'A_O3_Bird&Riordan1983_ndelta.csv', 
					A_UG = 'A_UG_Bird&Riordan1983_ndelta.csv'):

	'''
	Calling the values of Extraterrestrial Irradiation and Absorption of 
	Water, Ozone and Mixed Gases. Each in separate files, these are arrays 
	with the absorption and values, obtained from the tables in Bird and 
	Riordan [1]
	This is an attempt to do the simulation more resoluted, using the 
	ETR from Gueymard, 1985 [10].
	'''
	try:
		# in case of being run as code in console.		
		path = os.path.dirname(os.path.abspath(__file__)) 
		
	except Exception:
		# in case of being called as a module.		
		path = os.path.dirname(os.path.abspath('__file__')) 
	
	dfs = []
	for file in [ETR, A_H2O, A_O3, A_UG]:
		df = pd.read_csv(r'{}\{}'.format(path, file), index_col = 0, 
											header = 0)
		dfs += [df.iloc[:, 0]]
	
	WL_range = pd.Series(index = np.geomspace(0.3, 4.0, 2500))
	
	for i in dfs:
		i.index = WL_range.index
	
	ETR, A_H2O, A_O3, A_UG = dfs
	
	return ETR, A_H2O, A_O3, A_UG, WL_range
	
def OmegaWl(WL_range, omega, omega_p):
	
	# Aerosol single scattering albedo as a function of wavelength [1]. 
	# Eq. 3-16 pg 7.
	omega_wl = omega * np.exp(-omega_p * (np.power(np.log(WL_range.index 
											/ 0.4), 2)))
	omega_wl = pd.Series(omega_wl, index = WL_range.index, name = 'omega_wl')
	
	return omega_wl
	
def Tr(WL_range, AM_pr):
	
	# Rayliegh Transmission as function of wavelength, From: F. X. 
	# Kneizys [5]. Eq. 2-4, pg 2.
	T_R = (np.exp(-AM_pr / (np.power(WL_range.index , 4) * (115.6406 - 1.335 
								/ np.power(WL_range.index, 2)))))
	T_R = pd.Series(T_R, index = WL_range.index, name = 'T_R')
	
	return T_R
	
def To3(WL_range, A_O3, O3_m, O3):
	
	# Ozone Transmission as function of wavelength.  From: B. Leckner [4]. 
	# Eq. 2-9, pg 5.
	T_O3 = np.exp(-A_O3 * O3_m * O3)
	T_O3 = pd.Series(T_O3, name = 'T_O3')
		
	return T_O3
	
def Tug(WL_range, A_UG, AM_pr):
	
	# Uniform Mixed Gases Transmission as function of wavelength. 
	# From: B. Leckner [4]. Eq. 2-11, pg 5.
	#print(A_UG.head(), AM_pr)
	T_UG = (np.exp(-1.41 * A_UG * AM_pr / (np.power((1 + 118.93 * A_UG * 
											AM_pr), 0.45))))
	#return T_UG
	#print(T_UG.head())
	T_UG = pd.Series(T_UG, name = 'T_UG')
	
	return T_UG
		
def Th2o(WL_range, A_H2O, H2O, AM):
	
	# Water Vapor Transmission as function of wavelength. 
	# From: B. Leckner [4]. Eq. 2-8, pg 4.
	T_H2O = (np.exp(-0.2385 * A_H2O * H2O * AM / (np.power((1 + 20.07 * 
									A_H2O * H2O * AM), 0.45))))
	T_H2O = pd.Series(T_H2O, name = 'T_H2O')
		
	return T_H2O
	
def Tau(WL_range, AOD, alpha):
	
	# Aerosol turbity (tau) in function of wavelength [1]. Eq. 2-7 pag 4.
	tau = AOD * (np.power((WL_range.index / 0.5), -alpha))
	tau = pd.Series(tau, index = WL_range.index, name = 'tau')
	
	return tau
		
def Ta(WL_range, tau, AM):

	# Transmittance for Aerosol as function of wavelength. 
	# From A. Angstrom [6]. Eq. 2-6, pg 4.
	T_a = np.exp(-tau * AM)
	T_a = pd.Series(T_a, index = WL_range.index, name = 'T_a')

	return T_a
	
def Tas(omega_wl, tau, AM):
	
	# Transmittance for aerosol scattering [1]. Eq. 3-9, p. 7.
	T_as = pd.concat([omega_wl, tau], axis = 1).interpolate(
										method = 'linear')
	T_as['T_as'] = np.exp(-T_as.omega_wl * T_as.tau * AM)
	T_as = T_as.iloc[:, -1]
	
	return T_as
		
def Taa(omega_wl, tau, AM):
	
	# Transmittance for aerosol absorption [1]. Eq. 3-10, pg 7.
	T_aA = pd.concat([omega_wl, tau], axis = 1).interpolate(
										method = 'linear')
	T_aA['T_aA'] = np.exp(-(1 - T_aA.omega_wl) * T_aA.tau * AM)
	T_aA = T_aA.iloc[:, -1]

	return T_aA
		
def Edifr(ETR, T_O3, T_H2O, T_UG, T_aA, T_R, r_vec, Ang_zth):
	
	'''
	Rayleigh scattering component for diffuse irradiation, 
	from Bird [1]. Eq. 3-5 pag 7.
	'''
	I_dif_R = pd.concat([ETR, T_O3, T_H2O, T_UG, T_aA, T_R], axis = 1
								).interpolate(method = 'linear')

	I_dif_R['I_dif_R'] = (I_dif_R.ETR * r_vec * cos(rad(Ang_zth))
					* I_dif_R.T_O3 * I_dif_R.T_H2O * I_dif_R.T_UG 
					* I_dif_R.T_aA * (np.power((1 - I_dif_R.T_R), 
											0.95) * 0.5))
	I_dif_R = I_dif_R.iloc[:, -1]
	
	return I_dif_R

def Edifr2(ETR, T_O3, T_H2O, T_UG, T_aA, T_R, r_vec, Ang_zth):
	
	'''
	Rayleigh scattering component for diffuse irradiation, 
	from Bird [1]. Eq. 3-5 pag 7.
	'''
	I_dif_R = pd.DataFrame((ETR.values * r_vec * cos(rad(Ang_zth)) * T_O3.values * 
			T_H2O.values * T_UG .values * T_aA.values * (np.power(
				(1 - T_R.values), 0.95) * 0.5)), index = ETR.index, 
										columns = ['I_dif_R'])

		return I_dif_R
	
def Edifa(ETR, T_O3, T_H2O, T_UG, T_aA, T_R, T_as, Ang_zth, F_s, r_vec):

	# Aerosol scattering component for diffuse irradiation, from Bird [1]. 
	# Eq. 3-6 pag 7.
	I_dif_a = pd.concat([ETR, T_O3, T_H2O, T_UG, T_aA, T_R, T_as], axis = 1
							). interpolate(method = 'linear')
	I_dif_a['I_dif_a'] = (I_dif_a.ETR * cos(rad(Ang_zth)) * I_dif_a.T_O3 
					* I_dif_a.T_H2O * I_dif_a.T_UG * I_dif_a.T_aA 
								* (np.power(I_dif_a.T_R, 1.5)) 
							* (1 - I_dif_a.T_as) * F_s * r_vec)
	I_dif_a = I_dif_a.iloc[:, -1]
	
	return I_dif_a
		
def Trp(WL_range):	
	
	'''
	Atmospheric transmittance after Rayleigh scattering, from Bird [1]. 
	Eq. 2-4 pag 2.
	'''
	T_R_p = (np.exp(-1.8 / (np.power(WL_range.index, 4) * (115.6406 - 1.335 
								/ np.power(WL_range.index, 2)))))
	T_R_p = pd.Series(T_R_p, index = WL_range.index, name = 'T_R_p')
	
	return T_R_p
		
def Th2op(A_H2O, WL_range, H2O):
	
	'''
	Water vapour transmittance, from Bird [1]. Eq. 2-8 pag 4 
	(It's different equation due the terms evaluated at AM = 1.8)
	'''
	T_H2O_p = (np.exp(-0.2385 * A_H2O * H2O * 1.8 / (np.power((1 + 20.07 * 
									A_H2O * H2O * 1.8), 0.45))))
	T_H2O_p = pd.Series(T_H2O_p, index = WL_range.index, name = 'T_H2O_p')
	
	return T_H2O_p
		
def Tugp(A_UG, WL_range):	

	'''
	Uniformly mixed gas transmittance, from Bird [1]. Eq. 2-11 pag 5 
	(It's different equation due the terms evaluated at AM = 1.8) paper 
	indicates a value of 118.93 but remarks that the value 
	from Leckner [4] 118.3.
	'''
	T_UG_p = (np.exp(-1.41 * A_UG * 1.8 / (np.power((1 + 118.3 * A_UG * 1.8), 
												0.45))))
	T_UG_p = pd.Series(T_UG_p, index = WL_range.index, name = 'T_UG_p')
	
	return T_UG_p
		
def TaspTaap(omega_wl, tau):
	
	# Tasp (Prime transmittance terms evaluated at AM = 1.8).
	T_as_p = np.exp(-omega_wl * tau * 1.8)
	T_as_p = pd.Series(T_as_p, name = 'T_as_p')
	
	# Taap (Prime transmittance terms evaluated at AM = 1.8).
	T_aA_p = np.exp(-(1 - omega_wl) * tau * 1.8)
	T_aA_p = pd.Series(T_aA_p, name = 'T_aA_p')
	
	return T_as_p, T_aA_p
		
def Rhos(T_UG_p, T_H2O_p, T_aA_p, T_R_p, T_as_p, F_s_p):
	
	# Sky reflectivity (rho_g = ground albedo in function of wavelength), [1]
	# Eq. 3-8, pg 7.
	rho_s = pd.concat([T_UG_p, T_H2O_p, T_aA_p, T_R_p, T_as_p], axis = 1
								).interpolate(method = 'linear')
	rho_s['rho_s'] = (rho_s.T_UG_p * rho_s.T_H2O_p * rho_s.T_aA_p * (0.5 * 
		(1 - rho_s.T_R_p) + (1 - F_s_p) * rho_s.T_R_p * (1 - rho_s.T_as_p)))
	rho_s = rho_s.iloc[:, -1]
	
	return rho_s
	
def Edir(ETR, T_R, T_O3, T_UG, T_H2O, T_a, r_vec):
	
	# Direct solar irradiance.
	I_dir = pd.concat([ETR, T_R, T_O3, T_UG, T_H2O, T_a], axis = 1
								).interpolate(method = 'linear')
	I_dir['I_dir'] = (I_dir.ETR * r_vec * I_dir.T_R * I_dir.T_O3 
							* I_dir.T_UG * I_dir.T_H2O * I_dir.T_a)
	I_dir = I_dir.iloc[: , -1]
	
	return I_dir
		
def Edifg(I_dir, I_dif_R, I_dif_a, rho_s, rho_g, Ang_zth):

	# Diffuse solar irradiance fraction from the ground reflection, 
	# from Bird [1]. Eq. 3-7, pag 7.
	I_dif_g = pd.concat([I_dir, I_dif_R, I_dif_a, rho_s], axis = 1
								).interpolate(method = 'linear')
	I_dif_g['I_dif_g'] = ((I_dif_g.I_dir * cos(rad(Ang_zth)) 
							+ I_dif_g.I_dif_R + I_dif_g.I_dif_a) 
									* I_dif_g.rho_s * rho_g 
								/ (1 - I_dif_g.rho_s * rho_g))
	I_dif_g = I_dif_g.iloc[:, -1]

	return I_dif_g
	
def Cs(WL_range): 
	
	# Correction factor [1]. Eq. 3-17 p. 7.
	C_s = pd.DataFrame(WL_range)
	C_s['wl'] = C_s.index
	C_s['C_s'] = np.power((C_s.loc[C_s.index <= 0.45, 'wl'] + 0.55), 1.8)
	C_s['C_s'].fillna(value = 1)
	
	C_s = C_s.iloc[:, -1]
	
	return C_s
	
def EdifhWL(I_dif_R, I_dif_a, I_dif_g, C_s):
	
	'''
	Diffuse Irradiance in horizontal surface, UV correction factor 
	for > 0.45 is (wavelength + 0.55) ^ 1.8 for <= is 1. 
	This factor is called C_s in the paper.
	'''
	I_dif_h = pd.concat([I_dif_R, I_dif_a, I_dif_g, C_s], axis = 1
								).interpolate(method = 'linear')
	I_dif_h['I_dif_h'] = ((I_dif_h.I_dif_R + I_dif_h.I_dif_a 
								+ I_dif_h.I_dif_g) * I_dif_h.C_s)
	I_dif_h = I_dif_h.iloc[:, -1]

	return I_dif_h
	
def EdifcgWL(I_dir, I_dif_h, Ang_zth, rho_g, Ang_slp):
	
	# Ground reflected isotropic component on tilted surface.
	I_dif_c_g = pd.concat([I_dir, I_dif_h,], axis = 1
								).interpolate(method = 'linear')
	I_dif_c_g['I_dif_c_g'] = ((I_dif_c_g.I_dir * cos(rad(Ang_zth)) 
		+ I_dif_c_g.I_dif_h) * rho_g * (1 - 	cos(rad(Ang_slp))) / 2)
	I_dif_c_g = I_dif_c_g.iloc[:, -1]
	
	return I_dif_c_g
		
def EdifcscWL(I_dif_h, I_dir, ETR, r_vec, cos_AOI_sys, Ang_zth):

	# Solar Circumsolar component on tilted surface.
	I_dif_c_sc = pd.concat([I_dif_h, I_dir, ETR], axis = 1).interpolate(
										method = 'linear')
	I_dif_c_sc['I_dif_c_sc'] = (I_dif_c_sc.I_dif_h * (I_dif_c_sc.I_dir 
			/ I_dif_c_sc.ETR * r_vec) * cos_AOI_sys / cos(rad(Ang_zth)))
	I_dif_c_sc = I_dif_c_sc.iloc[:, -1]

	return I_dif_c_sc
	
def EdifcssiWL(I_dif_h, I_dir, ETR, r_vec, Ang_slp):

	# Solar Sky Isotropic component on tilted surface.
	I_dif_c_ssi = pd.concat([I_dif_h, I_dir, ETR], axis = 1).interpolate(
										method = 'linear')
	I_dif_c_ssi['I_dif_c_ssi']= ((I_dif_c_ssi.I_dif_h * (1 - (
				I_dif_c_ssi.I_dir / I_dif_c_ssi.ETR * r_vec)) 
								* (1 + cos(rad(Ang_slp)))) / 2)
	I_dif_c_ssi = I_dif_c_ssi.iloc[:, -1]
	
	return I_dif_c_ssi
	
def EdifslpWL(I_dif_c_g, I_dif_c_sc, I_dif_c_ssi):
	
	# Diffuse irradiation on tilted surface.
	I_dif_srf = pd.concat([I_dif_c_g, I_dif_c_sc, I_dif_c_ssi], 
				axis = 1, join = 'outer').interpolate(method = 'linear')
	I_dif_srf['I_dif'] = (I_dif_srf.I_dif_c_g + I_dif_srf.I_dif_c_sc 
									+ I_dif_srf.I_dif_c_ssi)
	# Dismisses all non-positive values.
	I_dif_srf['I_dif_srf'] = I_dif_srf.loc[I_dif_srf['I_dif'] >= 0, 'I_dif'] 
	I_dif_srf['I_dif_srf'] = I_dif_srf['I_dif_srf'].fillna(value = 0) 
	I_dif_srf = I_dif_srf.iloc[:, -1]
	
	return I_dif_srf
	
def EdirslpWL(I_dir, cos_AOI_sys):
	
	if cos_AOI_sys < 0: cos_AOI_sys = 0
	I_dir_srf = I_dir * cos_AOI_sys
	I_dir_srf.name = 'I_dir_srf'
	
	return I_dir_srf
	
def EetrrvslpWL(ETR, r_vec, cos_AOI_sys):
	
	if cos_AOI_sys < 0: cos_AOI_sys = 0
	ETR_rv_srf = ETR * r_vec * cos_AOI_sys
	ETR_rv_srf.name = 'ETR_rv_srf'
	
	return ETR_rv_srf
		
def Etotslp(I_dir_srf, I_dif_srf):
	
	I_tot_srf = pd.concat([I_dir_srf, I_dif_srf], axis = 1, join = 'outer'
								).interpolate(method = 'linear')
	I_tot_srf['I_tot_srf'] = (I_tot_srf.I_dir_srf + I_tot_srf.I_dif_srf)
	I_tot_srf = I_tot_srf.iloc[:, -1]
	
	return I_tot_srf
		
def Etot(I_dir, ETR, I_dif_srf, ETR_rv_srf, I_tot_srf, I_dir_srf):
	
	# Concatenation of irradiance spectra.
	ModIrr = pd.concat([ETR, ETR_rv_srf, I_dir, I_dir_srf, I_dif_srf, 
							I_tot_srf], axis = 1, join = 'outer'
								).interpolate(method = 'linear')
	ModIrr.index = ModIrr.index * 1000 # Converts from W/m2/um to W/m2/nm.
	ModIrr.index.names = ['nm']
	ModIrr = ModIrr / 1000
	
	return ModIrr