# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 11:36:00 2017

@author: Roberto

After losing the priviledges of the previous access. I need to download Jsc in stead of energy production

To access by browser the url (vpn):
ssfdb2(.shef.ac.uk)/phpmyadmin 
is required
username: pma
password: pmapw
"""

import pymysql as msl
import pandas as pd
import datetime as dt

def InfoTestBed(unit):

	start = dt.datetime(2015, 3, 15)
	end = dt.datetime(2016, 7, 31)
	rng = pd.date_range(start, end, freq = 'W')
	#id_fields = [[r'Total irradiance (W/$m^2$)', 1], [r'Diffuse irradiance (W/$m^2$)', 2], ['Power_54 (W)', 54],	['Power_55 (W)', 55], ['Power_58 (W)', 58], ['Power_64 (W)', 64], ['Power_125 (W)', 125]] ;  id_fields = [[j, i] for i, j in id_fields]
	DB_d = pd.read_csv(r'C:\Users\Roberto\Documents\PhD 2013\Experiments\Pyranometer\tblSSF_type (Descriptions of DB).csv', header = None, usecols = [0,1,2,3,4,5], names = ['id', 'name', 'detail', 'min', 'max', 'units'])
	#unit = 'degC'
	id_fields = list(DB_d.loc[(DB_d['detail'].str.contains('power') == False) & (DB_d['units'] == unit), ['id', 'detail']].values)
	
	try:
		connection = msl.connect(host='ssfdb2', user='rcahuantzi', password='RCahuantzi.123', db='SSFtestbed', cursorclass=msl.cursors.DictCursor)
	except Exception:									
		connection = msl.connect(host='ssfdb2.sheffield.ac.uk', user='rcahuantzi', password='RCahuantzi.123', db='SSFtestbed', cursorclass=msl.cursors.DictCursor)
		
 
		with connection.cursor() as cursor:
			print('Connection Established')
			TestBedDataRetrieving(cursor, start, end, rng, id_fields, unit)
				
def TestBedDataRetrieving(cursor, start, end, rng, id_fields, unit):	
	for i in range(len(rng))[:18:2]:
		yr1 = rng[i].year
		mt1 = rng[i].month
		dy1 = rng[i].day

		try:
			yr2 = rng[i+2].year
			mt2 = rng[i+2].month
			dy2 = rng[i+2].day
		
		except Exception:
			yr2 = rng[-1].year
			mt2 = rng[-1].month
			dy2 = rng[-1].day
						
		df0 = pd.DataFrame()
		#error_count = 1
		
		for j in id_fields:
			error_count = 1
			while True:	
				try:
					print('Retriving {}, for dates {:02.0f}-{:02.0f}-{:02.0f} to {:02.0f}-{:02.0f}-{:02.0f}, unit: {}'.format(	j[1], yr1, mt1, dy1, yr2, mt2, dy2, unit), end = '')
					#quit()
					sel = 'SELECT `type_id`, `data`, `timestamp`'
					frm = 'FROM `tblSSF_data`'
					where0 = 'WHERE `type_id`= {0}'.format(j[0])
					where1 = "&& `timestamp` BETWEEN date('{}-{}-{}')".format(yr1, mt1, dy1)
					where2 = "AND date('{}-{}-{}')".format(yr2, mt2, dy2)

					sql = '{} {} {} {} {}'.format(sel, frm, where0, where1, where2)
										
					cursor.execute(sql)
					result = cursor.fetchall()
					
					if len(result) == 0: 
						print('\t\x1b[4;31;48mNo data for {} found.\x1b[0m'.format(j[1])) ; break
					df = pd.DataFrame(result)
					#print(df.head())
					df.set_index('timestamp', drop = True, inplace = True)
					df = df.drop(['type_id'], axis = 1)
					#print(df.head())
					df.rename(columns = {df.columns[0]:j[1]}, inplace = True)
					
					df0 = df0.join(df, how = 'outer')
					print('\t\x1b[3;32;48mOK!\x1b[0m')
				except Exception:
					print('{} error(s) encountered, process continues...'.format(error_count), end = '\r')
					if error_count >= 10: print('\t\x1b[4;31;48mAfter {} attempts, no data for {} found.\x1b[0m'.format(error_count, j[1])) ; break # only 10 errors is the 
					error_count += 1
					
					continue # this makes ignore the "break" instance and go back to the "while True"
				break
				print('\x1b[3;32;48mProcess for {} finished.\x1b[0m'.format(j[1]))
										
		df0.to_csv('{0:02.0f}-{1:02.0f}-{2:02.0f}_to_{3:02.0f}-{4:02.0f}-{5:02.0f}_TestBed_of_unit_{6:}.csv'.format(yr1, mt1, dy1, yr2, mt2, dy2, unit))
		
		print('\t\x1b[4;32;48mProcess for {0:02.0f}-{1:02.0f}-{2:02.0f} to {3:02.0f}-{4:02.0f}-{5:02.0f} unit {6:} finished.\x1b[0m'.format(yr1, mt1, dy1, yr2, mt2, dy2, unit))

if __name__ == '__main__':
	#units = ['A', 'degC', 'mA', 'W', 'V']	
	#units = ['degC', 'mA', 'W', 'V', 'ms-1', 'mms-1']	
	#units = ['A', 'W', 'V', 'ms-1', 'mms-1'] # Problems with the retriving of A and W units
	units = ['A'] # Problems after date 21-06-2015; for i in range(len(rng))[16::2]:
	for unit in units:
		InfoTestBed(unit)

