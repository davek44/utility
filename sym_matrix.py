#!/usr/bin/env python
import numpy as np

################################################################################
# sym_matrix.py
#
# Space-efficient symmetric matrix class.
#
# Indexing adapted from here: http://stackoverflow.com/a/24563079/4114434
################################################################################

class sym_matrix:
	def __init__(self, n):
		self.n = n
		self.M = np.zeros(self.n*self.n)

	def get(self, i, j):
		if i < j:
			i, j = j, i
		return self.M[j*self.n - (j+1)*j/2 + i]

	def set(self, i, j, v):
		if i < j:
			i, j = j, i
		self.M[j*self.n - (j+1)*j/2 + i] = v
