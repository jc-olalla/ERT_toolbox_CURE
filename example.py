#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 22:10:28 2021

@author: jchavezolalla
"""
import toolbox as tb

# Input
x_i = 0.0
x_e = 39.0
ne = 40
filename_in = 'dipole-dipole.txt'
indexes_folder = 'indexes/'
indexes_files = ['line1_indexs_no_borehole.txt', 'line2_indexs_no_borehole.txt', 'line3_indexs_no_borehole.txt', 'line4_indexs_no_borehole.txt']
filename_out = 'kragge.txt'

# Swap indexes
tb.create_data_dd(x_i, x_e, ne, max_geom_factor=5000)
pos, abmn = tb.swap_indexes(indexes_files, filename_in, path_index='indexes/', path_base='')
tb.dataformatsyscal(pos, abmn, name=filename_out)










