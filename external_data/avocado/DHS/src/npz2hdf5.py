#!/usr/bin/env python

#- TRANSFORM npz files into hdf5 files so usable in R

import numpy as np
import h5py
import os 
import fnmatch

for root, dirs, files in os.walk('../data/'):
	for name in fnmatch.filter(files, '*.npz'):
		fname = '../data/'+name
		print(fname)
		scores = np.load(fname)
		key    = scores.keys()[0] #- only arr_0 in each file as it looks like
		print(fname + '.hdf5')
		h5f    = h5py.File(fname + '.hdf5', 'w')
		h5f.create_dataset(key, data=scores[key])
		h5f.close()



