#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 23:38:07 2021

@author: jchavezolalla
"""
import numpy as np

def geomfactor(xa, xb, xm, xn):
	ram = (xa - xm) ** 2
	ran = (xa - xn) ** 2
	rbm = (xb - xm) ** 2
	rbn = (xb - xn) ** 2
	rapm = (xa - xm) ** 2
	rapn = (xa - xn) ** 2
	rbpm = (xb - xm) ** 2
	rbpn = (xb - xn) ** 2

	H = 1 / ram - 1 / ran - 1 / rbm + 1 / rbn + 1 / rapm - 1 / rapn - 1 / rbpm + 1 / rbpn
	if H == 0.0:
		K = 0.0
	else:
		K = 4 * np.pi / H
# 	print(K)
	return K


def dataformatsyscal(pos, abmn, name='suveyfile.txt'):
	"""function to convert a data file into a format readable by electepro"""
	f = open(name, 'w')
	
	f.write('#  X   Y   Z \r\n')
	
	for i in range(len(pos)):
		x = '%.2f'%(pos[i])
		y = '%.2f'%(0.0)
		z = '%.2f'%(0.0)

		f.write(str(i + 1) + '  ' + str(x) + '  ' + str(y) + '   ' + str(z) + '\r\n')

	f.write('#  A   B   M   N \r\n')
	
	for i in range(np.shape(abmn)[0]):
		f.write(str(i + 1) + '  ' + str(int(abmn[i, 0]) + 1) + '  ' + str(int(abmn[i, 1]) + 1) + '  ' + str(int(abmn[i, 2]) + 1) + '  ' + str(int(abmn[i, 3]) + 1) + '\r\n')

	f.close()


def create_data_dd(x_ini, x_end, ne, max_geom_factor=5000):
	# dipole-dipole data generation
	# x_ini: position first electrode
	# x_end: position last electrode
	# ne: number of electrodes
	    
	pos = np.linspace(x_ini, x_end, ne)
	    
	a_list = []
	b_list = []
	m_list = []
	n_list = []
	
	d_AB_max = ne
	d_AB = 1
	while d_AB <= d_AB_max:
		d_BM = d_AB
		ne_temp = int(np.ceil(ne / d_AB))
		for j in range(ne_temp):
			d_BM = d_AB
			ne_temp = ne_temp - 1
			d_BM_max = (ne_temp - 3) * d_AB
			while d_BM < d_BM_max:
				a = j * d_AB 
				b = j * d_AB + d_AB
				m = j * d_AB + d_AB + d_BM
				n = j * d_AB + d_AB + d_BM + d_AB
				
				x_a = pos[a]
				x_b = pos[b]
				x_m = pos[m]
				x_n = pos[n]
				
				geo_factor = geomfactor(x_a, x_b, x_m, x_n)
				if np.abs(geo_factor) <= max_geom_factor:
					print(geo_factor)
					a_list.append(a)
					b_list.append(b)
					m_list.append(m)
					n_list.append(n)
				else:
					break
						
				d_BM = d_BM + d_AB

		d_AB = d_AB + 1

	abmn = np.zeros((len(a_list), 4))
	abmn[:, 0] = np.array(a_list)
	abmn[:, 1] = np.array(b_list)
	abmn[:, 2] = np.array(m_list)
	abmn[:, 3] = np.array(n_list)
	    
	dataformatsyscal(pos, abmn, name='dipole-dipole.txt')

def create_data_slm(x_ini, x_end, ne, max_geom_factor=5000):
	# It needs to be finished
	    
	pos = np.linspace(x_ini, x_end, ne)
	    
	a_list = []
	b_list = []
	m_list = []
	n_list = []
	
	d_AM_max = ne
	d_AM = 1
	while d_AM <= d_AM_max:
		d_MN = d_AM
		ne_temp = int(np.ceil(ne / d_AM))
		for j in range(ne_temp):
			d_MN = d_AM
			ne_temp = ne_temp - 1
			d_MN_max = (ne_temp - 3) * d_AM
			while d_MN < d_MN_max:
				a = j * d_AM 
				m = j * d_AM + d_AM
				n = j * d_AM + d_AM + d_MN
				b = j * d_AM + d_AM + d_MN + d_AM
				
				x_a = pos[a]
				x_b = pos[b]
				x_m = pos[m]
				x_n = pos[n]
				
				geo_factor = geomfactor(x_a, x_b, x_m, x_n)
				if np.abs(geo_factor) <= max_geom_factor:
					print(geo_factor)
					a_list.append(a)
					b_list.append(b)
					m_list.append(m)
					n_list.append(n)
				else:
					break
						
				d_MN = d_MN + d_AM

		d_AM = d_AM + 1

	abmn = np.zeros((len(a_list), 4))
	abmn[:, 0] = np.array(a_list)
	abmn[:, 1] = np.array(b_list)
	abmn[:, 2] = np.array(m_list)
	abmn[:, 3] = np.array(n_list)
	    
	dataformatsyscal(pos, abmn, name='schlumberger.txt')


def count_electrodes(filename):
	with open(filename, 'r') as f:
		count = 0
		for line in f:
			if line[0] == '#':
				count = count + 1
			if count == 2:
				break
			else:
				ne = line.split(sep=' ')[0]

	return int(ne)

def count_config(filename):
	with open(filename, 'r') as f:
		lines = f.read().splitlines()
		last_line = lines[-1]
		
	nconfig = int(last_line.split(sep=' ')[0])
	return nconfig

def check_indexes(list_of_files, ne, path=''):
	for file_i in list_of_files:
		count_line = 0
		with open(path + file_i, 'r') as f:
			for line in f:
				if line[0] ==  '#' or line.strip() == '':
					pass
				else:
					count_line = count_line + 1
					
		if count_line != ne:
			print('Index file does not coincide with the number of electrodes')
			
def largest_index(list_of_files, path=''):
	max_index = 0
	for file_i in list_of_files:
		with open(path + file_i, 'r') as f:
			for line in f:
				if line[0] == '#' or line.strip() == '':
					pass
				else:
					val_index = int(line.split(sep='\t')[0])
					if val_index > max_index:
						max_index = val_index 
	
	return max_index

def get_abmn(filename, path=''):
	count = 0
	a_list = []
	b_list = []
	m_list = []
	n_list = []
	with open(path + filename, 'r') as f:
		for line in f:
			if line[0] ==  '#':
				count = count + 1
			if count == 2:
				if line[0] == '#' or line.strip() == '':
					pass
				else:
					a_list.append(int(line.split()[1]))
					b_list.append(int(line.split()[2]))
					m_list.append(int(line.split()[3]))
					n_list.append(int(line.split()[4]))
					
	abmn = np.zeros((len(a_list), 4))
	abmn[:, 0] = np.array(a_list)
	abmn[:, 1] = np.array(b_list)
	abmn[:, 2] = np.array(m_list)
	abmn[:, 3] = np.array(n_list)
				
	return abmn
					
def find_new_abmn(a_ref, b_ref, m_ref, n_ref, index_ref, index_new):
	new_a = index_new[index_ref.index(a_ref)]
	new_b = index_new[index_ref.index(b_ref)]
	new_m = index_new[index_ref.index(m_ref)]
	new_n = index_new[index_ref.index(n_ref)]

	return new_a, new_b, new_m, new_n


def swap_indexes(list_of_files, filename_base, path_index='', path_base=''):
	ne_line = count_electrodes(filename_base)
	nconfig_line = count_config(filename_base)
	check_indexes(list_of_files, ne_line, path=path_index)
	ne_tot = largest_index(list_of_files, path=path_index)
	
	abmn_ref = get_abmn(filename_base, path_base)
	
	pos_out = []
	
	abmn_out = np.zeros((np.shape(abmn_ref)[0] * 4, 4))
	count = 0
	for file_i in list_of_files:
		pos_new = []
		pos_ref = []
		with open(path_index + file_i, 'r') as f:
			abmn_new = np.zeros_like(abmn_ref)
			for line in f:
				if line[0] ==  '#' or line.strip() == '':
					pass
				else:
					pos_out.append(int(line.split(sep='\t')[0]))
					
					pos_new.append(int(line.split(sep='\t')[0])) 
					pos_ref.append(int(line.split(sep='\t')[1].strip()))
					
		for ipos in range(np.shape(abmn_ref)[0]):
			a_ref = int(abmn_ref[ipos, 0])
			b_ref = int(abmn_ref[ipos, 1])
			m_ref = int(abmn_ref[ipos, 2])
			n_ref = int(abmn_ref[ipos, 3])
			new_a, new_b, new_m, new_n = find_new_abmn(a_ref, b_ref, m_ref, n_ref, pos_ref, pos_new)
			
			
			
			abmn_new[ipos, 0] = new_a
			abmn_new[ipos, 1] = new_b
			abmn_new[ipos, 2] = new_m
			abmn_new[ipos, 3] = new_n
						
		for i_abmn in range(np.shape(abmn_new)[0]):
			abmn_out[count, 0] = abmn_new[i_abmn, 0]
			abmn_out[count, 1] = abmn_new[i_abmn, 1]
			abmn_out[count, 2] = abmn_new[i_abmn, 2]
			abmn_out[count, 3] = abmn_new[i_abmn, 3]
			count = count + 1
						
	# Unused electrodes
	for i in range(1, ne_tot+1):
		print(i)
		if i not in pos_out:
 			pos_out.append(i)
	
	return pos_out, abmn_out
	
	
	
 def pseudo_section_with_depth(data, other_array=False, array=None, array_type='dd', interp_method='linear', scatter=True, cmin=None, cmax=None):
    '''function to plot pseudo sections at depth
    depth of investigation was obtained from Loke's practical guide
    array_type can be :
    - wa = wenner alpha
    - dd
    interp_method can be:
    - nearest
    - cubic
    - linear'''
    nconfig = len(data('a'))
    pos = np.array(data.sensorPositions())
    xe_list = []
    ze_list = []
    rhoa_list = []
    try:
        if other_array:
            rhoa_from_data = array
        else:
            rhoa_from_data = data.get('rhoa').array()
        scatter_only = False
        if all(res_v == 0 for res_v in rhoa_from_data):
            scatter_only = True
    except:
        print('No apparent resistivity values found')
        rhoa_from_data = data.get('a').array() * 0
        scatter_only = True

    for i in range(nconfig):
        a = int(data('a')[i])
        b = int(data('b')[i])
        m = int(data('m')[i])
        n = int(data('n')[i])

        x_c1 = pos[a, 0]
        y_c1 = pos[a, 1]
        x_c2 = pos[b, 0]
        y_c2 = pos[b, 1]
        x_p1 = pos[m, 0]
        y_p1 = pos[m, 1]
        x_p2 = pos[n, 0]
        y_p2 = pos[n, 1]

        if array_type == 'dd':
            a_dist = x_c2 - x_c1
            n_dist = (x_p1 - x_c2) / a_dist
            n_vals = np.array([1, 2, 3, 4, 5, 6, 7, 8])
            ze_over_a_vals = np.array([0.416, 0.697, 0.962, 1.220, 1.476, 1.730, 1.983, 2.236])
            ze_over_a = np.interp(n_dist, n_vals, ze_over_a_vals)
            ze = -ze_over_a * a_dist
            xe = (x_c1 + x_c2 + x_p1 + x_p2) / 4
            rhoa = rhoa_from_data[i]
        elif array_type == 'wa':
            a_dist = x_p1 - x_c1
            ze_over_a = 0.519
            ze = -ze_over_a * a_dist
            xe = (x_c1 + x_c2 + x_p1 + x_p2) / 4
            rhoa = rhoa_from_data[i]

        else:
            print('Array cannot be recognized')

        xe_list.append(xe)
        ze_list.append(ze)
        rhoa_list.append(rhoa)

    # Contour plot
    xe_array = np.array(xe_list).reshape(len(xe_list), 1)
    ze_array = np.array(ze_list).reshape(len(ze_list), 1)
    xze = np.hstack((xe_array, ze_array))
    rhoa_array = np.array(rhoa_list)
    npoints = 1000
    # npoints_x = int(np.abs((np.max(pos[:, 0]) - np.min(pos[:, 1]))) / (pos[1, 0] - pos[0, 0]))
    # npoints_z = int(np.abs((np.max(ze_array)) - np.min(ze_array)) / (pos[1, 0] - pos[0, 0]))

    if scatter_only:
        plt.figure()
        plt.scatter(xe_array, ze_array, c='k', alpha=1, marker='.')
        plt.title('Scatter plot')

    else:
        xi = np.linspace(np.min(xe_array), np.max(xe_array), npoints).reshape(npoints, 1)
        zi = np.linspace(np.min(ze_array), np.max(ze_array), npoints).reshape(npoints, 1)
        X, Z = np.meshgrid(xi, zi)

        grid_rhoa = griddata(xze, rhoa_array, (X, Z), method=interp_method)

        v = np.linspace(np.min(rhoa_array), np.max(rhoa_array), 100, endpoint=True)
        plt.figure()
        CS = plt.contourf(X, Z, grid_rhoa, v, vmin=cmin, vmax=cmax, cmap='jet')
        plt.colorbar(CS)
        plt.scatter(xe_array, ze_array, c='k', alpha=0.1, marker='.')
        plt.title('pseudo section ' + array_type)

        if scatter:
            plt.figure()
            plt.plot(rhoa_array, 'x')
            plt.title('Apparent resistivity')
	
	
def dd_generation(scheme, max_geom_factor):
    '''function to generate dipole-dipole scheme with varying dipole spacing
    geometric factor is approximated with linear topography'''

    ne = scheme.sensorCount()
    scheme.resize(1)
    a_list = []
    b_list = []
    m_list = []
    n_list = []

    db = 1
    while db < ne:
        dbm = db
        for i in range(ne):
            dbm = db
            for j in range(ne):
                a = i
                b = i + db
                m = i + db + dbm
                n = i + db + dbm + db

                scheme('a')[0] = a
                scheme('b')[0] = b
                scheme('m')[0] = m
                scheme('n')[0] = n
                try:
                    geo_factor = geomfactor(scheme, 0, printfactor=False)
                except:
                    geo_factor = max_geom_factor + 1


                if (np.abs(geo_factor) < max_geom_factor) and (n < ne):
                    a_list.append(a)
                    b_list.append(b)
                    m_list.append(m)
                    n_list.append(n)

                else:
                    break

                # dbm = dbm + 1 # the separation is to short
                dbm = dbm + db

        db = db + 1

    scheme.resize(len(a_list))
    for i in range(len(a_list)):
        scheme('a')[i] = a_list[i]
        scheme('b')[i] = b_list[i]
        scheme('m')[i] = m_list[i]
        scheme('n')[i] = n_list[i]
        scheme('valid')[i] = 1
	

# Remove overlap configurations
for iroll in range(nconfig_roll):
    print(iroll, ' configurations analyzed out of: ', nconfig_roll)
    xa_roll = scheme_roll.sensorPositions()[int(scheme_roll.get('a')[iroll])][0]
    xb_roll = scheme_roll.sensorPositions()[int(scheme_roll.get('b')[iroll])][0]
    xm_roll = scheme_roll.sensorPositions()[int(scheme_roll.get('m')[iroll])][0]
    xn_roll = scheme_roll.sensorPositions()[int(scheme_roll.get('n')[iroll])][0]
    for iref in range(nconfig_ref):
        xa_ref = scheme_ref.sensorPositions()[int(scheme_ref.get('a')[iref])][0]
        xb_ref = scheme_ref.sensorPositions()[int(scheme_ref.get('b')[iref])][0]
        xm_ref = scheme_ref.sensorPositions()[int(scheme_ref.get('m')[iref])][0]
        xn_ref = scheme_ref.sensorPositions()[int(scheme_ref.get('n')[iref])][0]

        if xa_roll == xa_ref and xb_roll == xb_ref and xm_roll == xm_ref and xn_roll == xn_ref:
            scheme_roll('valid')[iroll] = 0
            break

 
