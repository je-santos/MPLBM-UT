#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: JIA Haowei, South West Petroleum University
E-mail contact: nv4dll@outlook.com

This script is to convert .dat files to .mat for Skeleton3D to read.

'''
import time
import numpy as np
import scipy.io
import copy

def datreader(filename,reversebool):
	'''
	Read file and return 2 numpy arrays. 
    Parameters
    ----------
	filename : str
	**reversebool : bool
		By changing reversebool=True, the user can reverse the meaning of True and False in bool_array
    Returns
    ----------
	bool_array : numpy array, in which True is represent for matrix and False is for porous
	'''
	file = open(filename) #读取dat文件
	count = 0
	for index, line in enumerate(open(filename,'r')):
		count += 1
	dataslines = file.readlines()
	row = [ [] for i in range(count) ]
	for i in range (0,count):
		row[i] = dataslines[i].split()
	if reversebool is False:
		for i in range(len(row)):
			for j in range(len(row[i])):
				if row[i][j] == '1' or row[i][j] == '2':
					row[i][j] = True
				else:
					row[i][j] = False
	else:
		for i in range(len(row)):
			for j in range(len(row[i])):
				if row[i][j] == '1' or row[i][j] == '2':
					row[i][j] = False
				else:
					row[i][j] = True
	bool_array = np.array(row) #储存row为bool
	return bool_array

def reshapearray(data,x,y,z):

	arrayreshaped = data.reshape(x,y,z) #reshape data to 3-dimentional array
	return arrayreshaped

def savematfile(data,filenames,key='testvol'):
	'''
	Save data in mat format as 'filenames.mat'.
    Parameters
    ----------
	data : numpy array
	filenames : str
	key : str, set as 'testvol' for matlab program to read
	'''
	scipy.io.savemat(filenames+'.mat',{key:data}) 

if __name__ == "__main__" :

	x,y,z=10,10,10
	#x=y=z=315
	filename = '1'
	binorigin = datreader(filename+'.dat',True)
	savematfile(reshapearray(binorigin,x,y,z),filename)