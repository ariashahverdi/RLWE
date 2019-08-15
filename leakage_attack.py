#!/usr/bin/env sage -python

import sys
import numpy
import operator
import time
import os.path
import numpy as np
from collections import Counter
import multiprocessing as mp
import multiprocessing
import scipy.stats
import math 
from scipy.stats import binom

# The following is needed to interact with RLWE challenges
import ChallInstParser

from sage.all import *
from sage.modules.free_module_integer import IntegerLattice
from fpylll import *

from sage.matrix.constructor import random_unimodular_matrix
matrix_space = sage.matrix.matrix_space.MatrixSpace(ZZ, 4)
A = random_unimodular_matrix(matrix_space)


# Place holder for some of the values needed by every function
n = 1024
q = 12289
leakRate = 0.25

R = Integers(q)
Q = PolynomialRing(R, 'y')
S = Q.quo(Q.gen() ** (n) + 1, 'x')

# Place holder for Threshold 
# ** Important: Set Therehold according to table in paper ** #
Threshold = 0

# Default Omega for NTT transform used in NewHope
Omega = R(7)

# Parameter used in Binomial Sampling
BinomParam = 8

# Subtract two list
multisub = lambda a, b: map(operator.sub, a, b)

# Sampling Parameter of Error
mu = 0
sigma = math.sqrt(8)

# Sampling Distribution which is supported
SAMPLINGS = ['Gaussian', 'Binomial']
sampling = 'Gaussian'

# The Leakage patterns that are supported
PATTERNS = ['1mod8', '1-7mod16', '1-15mod16', '1mod16']
leakPattern = '1mod8'

# We look at RLWE Challenges and parameters by NewHope
MODES = ['Challenges','NewHope']
mode = 'Challenges'

# setting verbose = 1 will print all the details 
# when running the attack
verbose = 0
if verbose:
	def verboseprint(*args):
		# Print each argument separately so caller doesn't need to
		# stuff everything to be printed into a single string
		for arg in args:
			print arg,
		print
else:
	verboseprint = lambda *a: None  # do-nothing function

def find_omega():
	# Find all the 2n-th root of unity and return the list of them
	res = []
	for i in range(q):
		if sq_mult(i, n, q) == q - 1 and sq_mult(i, 2 * n, q) == 1:
			res.append(i)
	return res


def binom_sample():
	# sample n element which are from Binomial Distribution with given 
	# parameter in BinomParam
	return list(np.random.binomial( BinomParam * 2 , 0.5, n ) - BinomParam )


def gauss_sample():
	# Sample n element from Gaussian Distribution with given 
	# parameter of mu and sigma, mu is 0 in our application
	# We first sample from Gaussian distribution and then 
	# round it to closest integer
	return np.rint(np.random.normal(mu, sigma, n)).tolist()


def sq_mult(base, exponent, modulus):
	# Raising a number to power might take a long time,
	# so we use square and multiply algorithm to raise
	# a number to a power and take modulus at the end
	
	# Converting the exponent to its binary form
	binaryExponent = []
	while exponent != 0:
		binaryExponent.append(exponent % 2)
		exponent = exponent / 2
	
	# Applying of the square and multiply algorithm
	result = 1
	binaryExponent.reverse()
	for i in binaryExponent:
		if i == 0:
			result = (result * result) % modulus
		else:
			result = (result * result * base) % modulus
	return result


def create_a(a, size):
	# Create a curculant matrix from polynomial version of a
	l = [[0 for _ in range(size)] for _ in range(size)]
	l[size - 1] = a.list()[::-1]
	for i in range(size - 2, -1, -1):
		for j in range(size - 1):
			l[i][j] = l[i + 1][j + 1]
		l[i][size - 1] = -l[i + 1][0]
	return matrix(l)


def leak(sk, amount):
	# This function emulates the leakage
	# Given input secret key (sk) or error (e)
	# Based on given pattern it leaks to corresponding
	# corrdinates of input sk and for each leaked
	# coordinate output their compute what root of
	# unity they are evaluated on and output that
	# as well. The root of unity part is needed in
	# the later step for Lagrange Polynomial

	if leakPattern == '1-7mod16' or leakPattern == '1-15mod16':
		# The following is either 7 or 15
		idx2 = int(leakPattern.split('-')[1].split('mod')[0])
		i = log(int(1 / amount), 2)  
		step = int(1 / amount) * 2
		leakcor = int(n * amount)
		leak1 = []; leak2 = []; leak = []
		for idx in range(leakcor / 2):
			leak1.append((mod(sq_mult(Omega, 2 * idx * step + 1, q), q), sk[idx * step]))
			leak2.append((mod(sq_mult(Omega, 2 * idx * step + idx2, q), q), sk[idx * step + ((idx2-1)/2) ]))
		leak.append(leak1); leak.append(leak2)
		return leak

	elif leakPattern == '1mod8':
		i = log(int(1/amount),2)  
		step = int(1/amount)
		leakcor = int(n * amount) 
		leak = []
		for idx in range(leakcor):
			leak.append( (mod(sq_mult(Omega, 2 * idx * step + 1, q), q),sk[idx*step]) )
		return leak

	elif leakPattern == '1mod16':
		i = log(int(1/amount),2) 
		step = int(1/amount)*2
		leakcor = int(n * amount) 
		leak = []
		for idx in range(leakcor/2):
			leak.append( (mod(sq_mult(Omega, 2 * idx * step + 1, q), q),sk[idx*step]) )
		return leak	

	else:
		# Other leakage patterns are not covered
		print("*** INVALID LEAKAGE PATTERN (5) ***")

	return


def coeff_fix(coeff, amount):
	# Keep in mind that we have to change from mod type to int type
	# Also the polynomial might have less coeffiecent we append 0
	target_len = int(amount * n)
	zero_pad = [0 for _ in range(target_len - len(coeff))]
	return coeff + zero_pad


def value_fix2(val):
	# Simply create a list of coefficients where each
	# coefficient is from Field R_q
	res = []
	for i in val:
		res.append(R(i))
	return res


def value_fix(value):
	# Just create a list where change elements in the range
	# of [0,q] to elements in the range of [-q/2,q/2]
	temp = vector(ZZ, value).list()
	result = [0 for _ in range(len(value))]
	for idx, val in enumerate(temp):
		if val > q / 2:
			result[idx] = val - q
		else:
			result[idx] = val
	return result


def calc_prob(val):
	# For a given list which is called val we compute the probability 
	# of each element of that list based on either Binomial of Gaussian
	# Distribution and since elements of the list are independetent of 
	# each other the total probability of the list is the multiplication
	# of probability of each element
	res = 1.0
	if sampling == 'Gaussian':
		for i in val:
			if i > q / 2:
				res *= gauss_prob(int(i) - q)
			else:
				res *= gauss_prob(int(i))
	elif sampling == 'Binomial':
		for i in val:
			if i > q / 2:
				res *= binom_prob(int(i) - q)
			else:
				res *= binom_prob(int(i))
	return res


def gauss_prob(i):
	# Compute the probability of certain number appearing in Gaussian distribution
	return scipy.stats.norm(mu, sigma).cdf(float(i)+0.5) - scipy.stats.norm(mu, sigma).cdf(float(i)-0.5)


def binom_prob(i):
	# Compute the probability of certain number appearing in Binomial distribution
	if abs(i) > 16:
		return 0
	return binom.pmf( i, BinomParam * 2, 0.5, -BinomParam ) 


def ntt_worker(idx, nttMat):
	# This function will compute one row of NTT matrix <--------------------------------------
	temp = [0 for _ in range(n)]
	for j in range(n):
		temp[j] = mod(sq_mult(Omega, (2 * idx + 1) * j, q), q)
	nttMat[idx] = temp


def ntt_gen():
	# Create NTT transform which is just a matrix
	# Here we use parallel processing to speed up the function
	print('\n\nNTT Generation\n\n')
	l = []
	manager = mp.Manager()
	stride = 32
	for totidx in range(n / stride):
		nttMat = manager.dict()
		jobs = []
		for i in range(stride):
			p1 = multiprocessing.Process(target=ntt_worker, args=(totidx * stride + i, nttMat))
			jobs.append(p1)
			p1.start()

		for proc in jobs:
			proc.join()

		nttMat = nttMat.values()
		for idx in range(stride):
			l.append(nttMat[idx])
		print("{} out of {} ".format(totidx, n / stride))

	return matrix(l)


def load_ntt_matrix(challID):
	# Each challenge has it's own unique ID and for each challenge ID
	# We create different NTT transform. In this function we try to
	# load NTT transform of the challenge ID or if challenge ID is 
	# unique pattern 5050, it is just NewHope case. If the file does
	# not exist we create the NTT transform and save it for future
	ntt_filename = "ntt" + str(n) + "_ChallID" + str(challID) + ".txt"

	if os.path.exists(ntt_filename):
		# NTT transform exist so just load it
		ntt_file = open(ntt_filename, "r")
		nttList = []
		for i in range(n):
			temp_list = []
			for j in range(n):
				val = ntt_file.readline()
				temp_list.append(mod(int(val[:-1]), q))
			nttList.append(temp_list)
		nttMatrix = matrix(nttList)
		ntt_file.close()
	else:
		# Generate NTT transform and save it
		nttMatrix = ntt_gen()
		ntt_file = open(ntt_filename, "w")
		nttList = nttMatrix.list()
		for item in nttList:
			ntt_file.write("%s\n" % item)
		ntt_file.close()

	verboseprint('NTT Matrix is loaded\n')
	return nttMatrix


def create_basis():
	# In this subroutine given the mat which is in 
	# the form of
	# [1 ,\omega^{u} ,\omega^{2 \cdot u}, \cdots, \omega^{(n'-1) \cdot u]
	# We create a basis for it to be used in CVP subroutine.
	# Ths significance of this basis are that every point x
	# belongs to lattice generated by this basis, the 
	# result of mat . x = 0

	if leakPattern == '1-7mod16' or leakPattern == '1-15mod16':
		# The following is either 7 or 15
		idx2 = int(leakPattern.split('-')[1].split('mod')[0])

		# we just compute mat as explained above
		mat1 = []; mat2 = []; mat = []
		for i in range(8):
			mat1.append(sq_mult(Omega, 1 * i * n / 8, q))
			mat2.append(sq_mult(Omega, idx2 * i * n / 8, q))
		mat.append(mat1)
		mat.append(mat2)
		mat = matrix(mat)

		# Here we create basis, the importance of it is
		# explained above
		A1 = mat[0:2, 0:2]
		A2 = mat[0:2, 2:]
		A3 = A1.inverse() * (-A2)
		A4 = q * matrix(matrix.identity(2))
		A5 = matrix.zero(6, 2)
		A6 = matrix.identity(6)
		A4 = A4.augment(A3)
		A5 = A5.augment(A6)
		basis = A4.transpose().augment(A5.transpose())

		B = IntegerMatrix(8, 8)
		
		# It should give something as following for NewHope parameters
		'''
		basis = [[12289, 0, 0, 0, 0, 0, 0, 0],
			[0, 12289, 0, 0, 0, 0, 0, 0],
			[12288, 5439, 1, 0, 0, 0, 0, 0],
			[5439, 9190, 0, 1, 0, 0, 0, 0],
			[9190, 392, 0, 0, 1, 0, 0, 0],
			[392, 3099, 0, 0, 0, 1, 0, 0],
			[3099, 5439, 0, 0, 0, 0, 1, 0],
			[5439, 1, 0, 0, 0, 0, 0, 1]]
		basis1 =[[12289,     0, 12288,  5439,  9190,   392,  3099,  5439],
			[    0, 12289,  5439,  9190,   392,  3099,  5439,     1],
			[    0,     0,     1,     0,     0,     0,     0,     0],
			[    0,     0,     0,     1,     0,     0,     0,     0],
			[    0,     0,     0,     0,     1,     0,     0,     0],
			[    0,     0,     0,     0,     0,     1,     0,     0],
			[    0,     0,     0,     0,     0,     0,     1,     0],
			[    0,     0,     0,     0,     0,     0,     0,     1]]
		'''		

	elif leakPattern == '1mod8':
		# Similar to previous case
		mat = []
		for i in range(int(1/leakRate)):
			mat.append((Omega**1)**(i*int(leakRate*n)))
		mat = matrix(mat)
		verboseprint("mat is ", mat)

		A1 = mat[0,0]
		A2 = mat[0,1:]
		# In this case A1 is just 1
		A3 = A1*(-A2)
		A4 = q * matrix(matrix.identity(1))
		A5 = matrix.zero(3, 1)
		A6 = matrix.identity(3)
		A4 = A4.augment(A3)
		A5 = A5.augment(A6)
		basis = A4.transpose().augment(A5.transpose())

		B = IntegerMatrix(4, 4)
		
		# It should give something like this for NewHope
		'''
		basis = [[12289, 0, 0, 0],
					[5146, 1,   0, 0],
					[1479, 0,  1, 0],
					[8246, 0,   0, 1]]


			basis1 = [[12289,  5146, 1479, 8246],
				 [    0, 	 1,    0,      0],
				 [    0,     0,     1,     0],
				 [    0,     0,     0,     1]]
		'''

	elif leakPattern == '1mod16':
		# Similar to previous case
		mat = []
		for i in range(8):
			mat.append(sq_mult(Omega, 1 * i * n / 8, q))
		mat = matrix(mat)
		verboseprint("mat is ", mat)

		A1 = mat[0,0]
		A2 = mat[0,1:]
		# In this case A1 is just 1
		A3 = A1 * (-A2)
		A4 = q * matrix(matrix.identity(1))
		A5 = matrix.zero(7, 1)
		A6 = matrix.identity(7)
		A4 = A4.augment(A3)
		A5 = A5.augment(A6)
		basis = A4.transpose().augment(A5.transpose())

		B = IntegerMatrix(8, 8)		

	else:
		print("*** Invalid Leakage Pattern (1) ***")
		return

	# In case for 1 mod 8, sage some time gave an error so,
	# we thought CVP in fplll might have a bug
	# So we randomly change the basis using random unimodular matrix
	# The following two lines are for that
	#basis1 = basis1 * A
	#basis = basis1.transpose()

	B.set_matrix(basis)
	return [mat, B]



def worker(sysIdx, status, aMat, uMat, countStatus, mat, B, eCoeff1, eCoeff2, eList, uTarget, aMatrix):
	# This function is designed to be potentially called in parallel
	# In this subroutine we find the correct error using << CVP >>

	# We initialize target which is either two values
	# or one value given the leakage pattern
	# The numVariables is simply how how many variables
	# the system we are solving for has
	target = []
	if leakPattern == '1-7mod16' or leakPattern == '1-15mod16':
		target.append(eCoeff1[sysIdx])
		target.append(eCoeff2[sysIdx])
		numVariables = 8

	elif leakPattern == '1mod8':
		target.append(eCoeff1[sysIdx])
		numVariables = 4 

	elif leakPattern == '1mod16':
		target.append(eCoeff1[sysIdx])
		numVariables = 8		
		
	else:
		# Other leakage patterns are not supported
		print("*** Invalid Leakage Pattern (2) ***")	

	# Create a vector of the target and the 
	# total number of system of equations
	target = vector(target)
	numEqns = n / numVariables

	# Here as a part of simulation we assume the 
	# Correct error is given so everytime we find 
	# a candidate which passes all our checkes,
	# we check it with correct error and if they
	# were not equal we immediately end the run
	correctError = []; correctUTarget = []
	for i in range(numVariables):
		idx = (i * numEqns) + sysIdx
		correctError.append(eList[idx])
		correctUTarget.append(uTarget[idx])	

	# Solve the system of equation for X and change it to vector
	X = mat.solve_right(target)
	X = vector(ZZ, X)

	# Run CVP on X to get closest point on lattive to X
	try:
		res = CVP.closest_vector(B, X)
	except:
		res = [n for _ in range(numVariables)]
	
	# Change the res to vector and bring it to feild R_q
	res = vector(res)
	answer = list(vector(R, X - res))

	# Place holder for some values
	aTemp = []
	uTemp = []
	countTemp = 0
	flag = True

	# Uncomment this if you want to see correct answer 
	# and the recovered one
	# print(sysIdx,target,correctError,calc_prob(correctError),answer,calc_prob(answer))

	if calc_prob(answer) > Threshold:
		# As I said in simulation if the answer which is found
		# is not exaclty correct we end the simulation by setting
		# the following flag to False
		if answer != correctError:
			flag = False
			print ("break")
		
		# We increate the number of variables found so far
		countTemp += numVariables

		# For each correct answer we save the rows of matrix a
		# as well as rows of vector correctUTarget which is simply
		# u in the paper. 
		# This values are going to be used once we are solving
		# for secret key (using Gaussian elimination) in the
		# last step of algorithm 
		for i in range(numVariables):
			idx = (i * numEqns) + sysIdx
			aTemp.append(aMatrix[idx][:])
		adjustedUTarget = list(vector(correctUTarget) - vector(answer))

		for i in adjustedUTarget:
			uTemp.append(i)

	# We are given an array of ouputs and each instance of
	# worker is working on one particular system of 
	# equation so it just simply update that part of 
	# output which it assigned to work on.
	status[sysIdx] = flag
	aMat[sysIdx] = aTemp
	uMat[sysIdx] = uTemp
	countStatus[sysIdx] = countTemp

	return

def worker2(sysIdx, status, aMat, uMat, countStatus, mat, B, eCoeff1, eCoeff2, eList, uTarget, aMatrix):
	# This function is designed to be potentially called in parallel
	# In this subroutine we find the correct error using << BruteForce >>
	# The approach we use in every part of it is meet-in-the-middle approach
	# where we break the systems into two independet one and check where
	# an answer satisfies both of them

	if leakPattern == '1-7mod16' or leakPattern == '1-15mod16':
		target = []
		target.append(eCoeff1[sysIdx])
		target.append(eCoeff2[sysIdx])
		target = vector(target)

		# Our brute force only search for element upto this value
		rangeGauss = (int(2 * sigma) * 2) + 1
		mat1 = mat[0][:]; mat2 = mat[1][:]

		# first part of the meet in the middle
		# We reconstruct the partial answers for half
		# of the system and save it in the table
		mat_temp = []
		mat1_temp = mat1[0:4]; mat2_temp = mat2[0:4]
		mat_temp.append(mat1_temp); mat_temp.append(mat2_temp)
		mat_temp = matrix(mat_temp)
		val = [0 for _ in range(4)]
		table = {}
		for i0 in range(rangeGauss):
			for i1 in range(rangeGauss):
				for i2 in range(rangeGauss):
					for i3 in range(rangeGauss):
						val[0] = i0 - ((rangeGauss - 1) / 2)
						val[1] = i1 - ((rangeGauss - 1) / 2)
						val[2] = i2 - ((rangeGauss - 1) / 2)
						val[3] = i3 - ((rangeGauss - 1) / 2)
						t = mat_temp * vector(val)
						table[tuple(t.list())] = tuple(val)

		# In the following we do several things, first 
		# we complete the meet in the middle part and for
		# the rest of systems we compute the partial answers
		# We then look at the table and check if the partial
		# answer given in this step combined with an answer
		# with the previous step give us the correct result
		# If that is the case we compute the probability of 
		# this answer. Finally we get an answer with 
		# highest probability
		answer = [0 for _ in range(8)]
		totalProb = 0
		maxProb = 0
		tempProb = 0
		mat_temp = []
		mat3_temp = mat1[4:8]; mat4_temp = mat2[4:8]
		mat_temp.append(mat3_temp); mat_temp.append(mat4_temp)
		mat_temp = matrix(mat_temp)
		for i0 in range(rangeGauss):
			for i1 in range(rangeGauss):
				for i2 in range(rangeGauss):
					for i3 in range(rangeGauss):
						val[0] = i0 - ((rangeGauss - 1) / 2)
						val[1] = i1 - ((rangeGauss - 1) / 2)
						val[2] = i2 - ((rangeGauss - 1) / 2)
						val[3] = i3 - ((rangeGauss - 1) / 2)
						t = mat_temp * vector(val)
						t = target - t
						if tuple(t.list()) in table:
							val2 = table[tuple(t)]
							temp = list(val2)
							temp.append(val[0]); temp.append(val[1])
							temp.append(val[2]); temp.append(val[3])
							result = value_fix2(temp)
							tempProb = calc_prob(result)
							if tempProb > maxProb:
								maxProb = tempProb
								answer = result

		# Compute the correct value for error
		# Since we are simulating the attack, if
		# the answer we found which passes all our
		# checkes is not correct we terminate the run
		correctError = [eList[sysIdx], eList[n / 8 + sysIdx], eList[2 * n / 8 + sysIdx], eList[3 * n / 8 + sysIdx],
						eList[4 * n / 8 + sysIdx], eList[5 * n / 8 + sysIdx], eList[6 * n / 8 + sysIdx],
						eList[7 * n / 8 + sysIdx]]
		correctUTarget = [uTarget[sysIdx], uTarget[n / 8 + sysIdx], uTarget[2 * n / 8 + sysIdx],
						  uTarget[3 * n / 8 + sysIdx], uTarget[4 * n / 8 + sysIdx], uTarget[5 * n / 8 + sysIdx],
						  uTarget[6 * n / 8 + sysIdx], uTarget[7 * n / 8 + sysIdx]]


		# It is going to be used for printnig
		target1 = mat * vector(correctError)
		target2 = mat * vector(answer)

		# Place holder for some values
		aTemp = []
		uTemp = []
		countTemp = 0
		flag = True

		# Uncomment this if you want to see correct answer 
		# and the recovered one
		# print(sysIdx, target, target1, target2, correctError, calc_prob(correctError), answer, calc_prob(answer))
			
		# Here we just check if the answer which was outputted
		# from previous step is passing the threshold
		# If it passes and not equal to the correct answer
		# we raise a flag and then terminate later

		if calc_prob(answer) > Threshold:
			if answer != correctError:
				flag = False
				#print(sysIdx, target, target1, target2, correctError, calc_prob(correctError), answer, calc_prob(answer))
			countTemp += 8

			for i in range(8):
				aTemp.append(aMatrix[(i * n) / 8 + sysIdx][:])
			adjustedUTarget = list(vector(correctUTarget) - vector(answer))

			for i in adjustedUTarget:
				uTemp.append(i)

	elif leakPattern == '1mod8':
		# The following is very similar to the previous case
		# The only difference is the way the meet in the middle
		# part works
		target = []
		target.append(eCoeff1[sysIdx])
		target = vector(target)

		correctError = [ eList[sysIdx] , eList[n/4 + sysIdx], eList[2*n/4 + sysIdx],  eList[3*n/4 + sysIdx]]	
		correctUTarget = [ uTarget[sysIdx] , uTarget[n/4 + sysIdx], uTarget[2*n/4 + sysIdx],  uTarget[3*n/4 + sysIdx]]

		maxProb = 0
		bound = int(round(2 * sigma))
		matList = mat.list()
		guess = []
		for iii in range(-bound,bound):
			for jjj in range(-bound,bound):
				for kkk in range(-bound,bound):
					RES = []
					RES.append(R(iii))
					RES.append(R(jjj))
					RES.append(R(kkk))
					tempRES = ( eCoeff1[sysIdx] - ( matList[0] * RES[0] + matList[1] * RES[1] + matList[2] * RES[2]) ) / matList[3]
					RES.append(tempRES)
					# Norm can be computed as follows
					#normResult = float(vector(value_fix(RES)).norm(2))
					tempProb = calc_prob(RES)
					if tempProb > maxProb:
						maxProb = tempProb
						guess = RES
						
		answer = guess

		aTemp = []
		uTemp = []
		countTemp = 0
		flag = True

		if calc_prob(answer) > Threshold:
			if answer != correctError:
				flag = False
				#print(sysIdx,target,correctError,calc_prob(correctError),answer,calc_prob(answer))
			countTemp += 4
			for i in range(4):
				aTemp.append(aMatrix[(i * n) / 4 + sysIdx][:])
			adjustedUTarget = list(vector(correctUTarget)-vector(answer))
			for i in adjustedUTarget:
				uTemp.append(i)

	elif leakPattern == '1mod16':
		# The following is very similar to the first case
		# The only difference is the way the meet in the middle
		# part works
		target = []
		target.append(eCoeff1[sysIdx])
		target = vector(target)

		correctError = [ eList[sysIdx], eList[n/8 + sysIdx], eList[2*n/8 + sysIdx], eList[3*n/8 + sysIdx],
						 eList[4*n/8 + sysIdx], eList[5*n/8 + sysIdx], eList[6*n/8 + sysIdx],
						 eList[7*n/8 + sysIdx]]	
		correctUTarget = [ uTarget[sysIdx] , uTarget[n/8 + sysIdx], uTarget[2*n/8 + sysIdx],  uTarget[3*n/8 + sysIdx],
						   uTarget[4*n/8 + sysIdx] , uTarget[5*n/8 + sysIdx], uTarget[6*n/8 + sysIdx], 
						   uTarget[7*n/8 + sysIdx]]

		mat1_temp = []
		matList = mat.list()
		mat_temp= matList[0:4] 
		mat_temp = matrix(mat_temp)
		val = [0 for _ in range(4)]
		table = {}
		rangeGauss = max( (int(2 * sigma) * 2) + 1, 3)

		for i0 in range(rangeGauss):
			for i1 in range(rangeGauss):
				for i2 in range(rangeGauss):
					for i3 in range(rangeGauss):
						val[0] = i0-((rangeGauss-1)/2); val[1] = i1-((rangeGauss-1)/2)
						val[2] = i2-((rangeGauss-1)/2); val[3] = i3-((rangeGauss-1)/2)
						t = mat_temp * vector(val)
						table[tuple(t.list())] = tuple(val)

		totalProb = 0; correctProb = 0
		mat_temp = matList[4:8]
		mat_temp = matrix(mat_temp)
		maxProb = 0
		for i0 in range(rangeGauss):
			for i1 in range(rangeGauss):
				for i2 in range(rangeGauss):
					for i3 in range(rangeGauss):
						val[0] = i0-((rangeGauss-1)/2); val[1] = i1-((rangeGauss-1)/2)
						val[2] = i2-((rangeGauss-1)/2); val[3] = i3-((rangeGauss-1)/2)
						t = mat_temp * vector(val)
						t = target - t
						if tuple(t.list()) in table:
							val2 = table[tuple(t)]
							temp = list(val2)
							temp.append(val[0]); temp.append(val[1])
							temp.append(val[2]); temp.append(val[3])
							result = value_fix2(temp)
							tempProb = calc_prob(result)
							if (tempProb > Threshold):
								totalProb += tempProb
								if(tempProb > maxProb):
									maxProb = tempProb
									guess = result
		answer = guess

		aTemp = []
		uTemp = []
		countTemp = 0
		flag = True
		minEntMaxRes = maxProb / totalProb
		minEntCorr	= calc_prob(correctError) / totalProb
		#print(sysIdx,target,correctError,calc_prob(correctError),answer,calc_prob(answer), minEntCorr, minEntMaxRes)		
		if calc_prob(answer) > Threshold:
			if answer != correctError:
				flag = False
				#print(sysIdx,target,correctError,calc_prob(correctError),answer,calc_prob(answer))		
			countTemp += 8
			for i in range(8):
				aTemp.append(aMatrix[(i * n) / 8 + sysIdx][:])
			adjustedUTarget = list(vector(correctUTarget) - vector(answer))

			for i in adjustedUTarget:
				uTemp.append(i)

	else:
		print("*** Invalid Leakage Pattern (3) ***")
		return

	# We are given an array of ouputs and each instance of
	# worker is working on one particular system of 
	# equation so it just simply update that part of 
	# output which it assigned to work on.
	status[sysIdx] = flag
	aMat[sysIdx] = aTemp
	uMat[sysIdx] = uTemp
	countStatus[sysIdx] = countTemp


def main():
	
	# Set the path for reading the challenges
	if len(sys.argv) < 8 or len(sys.argv) > 9:
		print "Usage:", sys.argv[0], "path/to/.challenge [path/to/.instance] [path/to/.secret] result.txt"
		sys.exit(-1)
	# Each Challenge in RLWE has 3 different files,
	# challenge file which is has (a) and (as+e)
	# secret file which has (s)
	# instance file which has the parameter set for that challenge
	challenge_path = sys.argv[1]
	instance_path = sys.argv[2]
	secret_path = sys.argv[3]
	
	# Leakage Pattern which can be either of the following patterns
	# 1 mod 8 | 1 mod 16 | 1-7 mod 16 | 1-15 mod 16
	global leakPattern
	leakPattern = sys.argv[4]

	# Depending on the pattern we set the Threshold as follows
	# Threhold is being used as a factor which determines whether
	# a recovered answer should be accepted or not
	
	# ** Important: Set Therehold according to table in paper ** #
	
	global Threshold
	if leakPattern == '1-7mod16' or leakPattern == '1-15mod16': 
		Threshold = 7e-6;
	elif leakPattern == '1mod8':
		Threshold = 7e-5;
	elif leakPattern == '1mod16':
		Threshold = 1e-8;
	else:
		Threshold == 0
	
	# Either Challenges or NewHope
	global mode
	mode = sys.argv[5] 

	# Sampling is Either Gaussian or Binomial
	global sampling
	sampling = sys.argv[6]

	global n, q, sigma
	if mode == 'Challenges':
		# Getting the parameters from challenges files
		[challID, AChallenge, BChallenge, s, n, q, svar, numSample] = ChallInstParser.get_challenge(challenge_path, instance_path, secret_path)
		# The value saved in the instance path it svar 
		# In the following we convert it to sigma which
		# then can be used as a parameter for sampling
		# from Gaussian distribution
		sigma = float(sqrt(svar * n) / (sqrt( 2 * math.pi )))

	elif mode == 'NewHope':
		# Setting parameters according the NewHope (USENIX version)
		n = 1024
		q = 12289
		sigma = math.sqrt(8)
		# Some big number of RLWE instances
		numSample = 1000
		# Setting ChallID a unique value 5050 for NewHope
		# This is just to make sure it is something different
		# from challID used in RLWE Challenges
		# We save the NTT matrix according to challID name,
		# This make sure we don't mix up the NTT transforrms
		challID = 5050

	else:
		# No other mode is supported
		print("*** Unknown Mode ***")

	# Print the values
	print('Mode ', mode)
	print('Leakage Pattern ', leakPattern)
	print('Sampling', sampling)
	print("n : ", n, " q : ", q, " sigma : ", sigma)	
	print("Threshold", Threshold)
	print("ratio", str(sigma/q))
	print("ChallID", str(challID))

	# Creating the Polynial Rings and Quotient Rings	
	global R
	R = Integers(q)
	global Q
	Q = PolynomialRing(R, 'y')
	global S
	S = Q.quo(Q.gen() ** (n) + 1, 'x')

	print ("\n\nRunning Entropy Attack ...\n\n")

	# w is the list of all 2n-th root of unity
	# We set Omega to be one of the w (e.g. the first one)
	w = find_omega()
	global Omega
	Omega = R(w[0])

	# We will run each system in parallel 
	# This shows how many processor your system has
	print("Number of processors: ", mp.cpu_count())

	# We load the NTT matrix from the file
	# If it does not exist we create on the fly
	# and save it for later
	nttMatrix = load_ntt_matrix(challID)

	# Create basis 
	[mat, B] = create_basis()

	if mode == 'Challenges':
		# If mode is Challenge wehave already loaded the secret key
		# So just save represent it by S (defined above)
		sk = S(s)
	
	elif mode == 'NewHope':
		# If we are generating secret key by ourselves then 
		# we sample the secret key either from Binomial Distribution
		# or Gaussian Distribution and represent it by S
		if sampling == 'Gaussian':
			sk = S(gauss_sample())
		elif sampling == 'Binomial':
			sk = S(binom_sample())
		else:
			# No other sampling method is supported
			print("*** Unknown Sampling Method ***")

	else:
		# No other mode is supported
		print("*** Wrong Mode ***")

	# As we go over the attack we recover some error values
	# We save correponding u and a values in the following list
	uTotal = []
	aTotal = []

	# count keep tracks of how many error coordinates we 
	# have recovered so far, the attack is done and we go 
	# to gaussian elmination part once count >= n
	count = 0

	# We have access to only numSample of RLWE instances
	# For each instance we try to recover some error coordinates
	for measidx in range(numSample):

		# Just print the sample index we are considering at the moment
		print("Sample idx : ", measidx)

		# If there is an error is some of the answers we found 
		# then flag will be False and recovering that secret key 
		# will not be successful
		flag = True

		# Change the secret key to the list and then 
		# apply NTT transform on it
		skList = sk.list()
		skNTT = nttMatrix * vector(sk)

		if mode == 'Challenges':
			# In the challenge mode everything is already loaded
			# So we just represent them in S
			a = S(AChallenge[measidx])
			# The error is not reported so we recover error from 
			# b and a and sk which are reported
			b = S(BChallenge[measidx])
			e = b - a * sk

		elif mode == 'NewHope':
			# We randomly sample public value a 
			a = S.random_element()
			# We sample error from the Distribution and 
			# compute the value of b
			if sampling == 'Gaussian':
				e = S(gauss_sample())
			elif sampling == 'Binomial':
				e = S(binom_sample())
			else:
				# No other sampling method is supported
				print("*** Unknown Sampling Method ***")
			b = a * sk + e

		else:
			# No other mode is supported
			print("*** Wrong Mode ***")

		# Create the matrix form of public key a
		aMatrix = create_a(a, n)

		# Get the error in list format
		eList = e.list()

		# Apply NTT transform to error e
		eNTT = nttMatrix * vector(e)

		# Compute u = a * sk + e (not going to use it)
		u = a * sk + e

		# uPrime is same as u but in NTT form
		uPrime = aMatrix * vector(sk) + vector(e)
		uTarget = uPrime[:]

		#############################################
		verboseprint("\nSecret Key : ", sk)
		verboseprint("\nSecret Key List : ", skList)
		verboseprint("\nSk After NTT is", skNTT)
		verboseprint("\nError is : ", e)
		verboseprint("\nError List : ", eList)
		verboseprint("\nError After NTT is", eNTT)
		verboseprint("\nChallenge B is", b)
		verboseprint("\nu is", u)
		#############################################

		# We leak on error! The logic is that the error
		# on secret key in NTT format is equivalent to 
		# leaking on the error in NTT format because we
		# know the public value u
		ePoints = leak(eNTT, leakRate)

		#############################################
		verboseprint("\nError Point Leakage : ", ePoints)
		#############################################

		# Using the leaked value we reconstruct the 
		# polynomial which is equivalent to taking 
		# error modulus
		if leakPattern == '1-7mod16' or leakPattern == '1-15mod16':
			ePoly1 = Q.lagrange_polynomial(ePoints[0])
			eCoeff1 = coeff_fix(ePoly1.list(), leakRate / 2)

			ePoly2 = Q.lagrange_polynomial(ePoints[1])
			eCoeff2 = coeff_fix(ePoly2.list(), leakRate / 2)

		elif leakPattern == '1mod8' or leakPattern == '1mod16':
			ePoly1 = Q.lagrange_polynomial(ePoints)
			eCoeff1 = coeff_fix( ePoly1.list(), leakRate )

			ePoly2 = []
			eCoeff2 = [] 

		#############################################
		verboseprint("\nError P1 Poly is ", ePoly1)
		verboseprint("\nError P1 coeff are ", eCoeff1)
		verboseprint("\nError P2 Poly is ", ePoly2)
		verboseprint("\nError P2 coeff are ", eCoeff2)
		#############################################

		Start_Time = time.time()

		# For each leakage case compute how many 
		# system of equations we have
		if leakPattern == '1-7mod16' or leakPattern == '1-15mod16' or leakPattern == '1mod16':
			numEqns = n / 8
		elif leakPattern == '1mod8':
			numEqns = n / 4
		else:
			print("*** Invalid Leakage Pattern (4) ****")

		# The following are for handling the parallel processing
		# We pass dictionaries as a place holder to get the 
		# outputs from the worker functions
		manager = mp.Manager()
		aMat = manager.dict()
		uMat = manager.dict()
		status = manager.dict()
		countStatus = manager.dict()
		jobs = []

		# ******** Note ********* #
		# We have two choice to recoer the potential answers
		# either using CVP routine (worker) or 
		# brute force routine (worker2)
		# We always use the worker (CVP) as our main solver
		# with following exception of the following leakage patterns
		# [1 mod 16, 1 mod 8]
		# for the above cases we use worker2 as our main sovlver
		for sysIdx in range(numEqns):
			if leakPattern == '1mod8' or leakPattern == '1mod16':
				p = multiprocessing.Process(target=worker2, args=(
					sysIdx, status, aMat, uMat, countStatus, mat, B, eCoeff1, eCoeff2, eList, uTarget, aMatrix))
			else:
				p = multiprocessing.Process(target=worker, args=(
					sysIdx, status, aMat, uMat, countStatus, mat, B, eCoeff1, eCoeff2, eList, uTarget, aMatrix))
			
			jobs.append(p)
			p.start()

		for proc in jobs:
			proc.join()

		# Get the values computed by each worker
		status = status.values()
		aMat = aMat.values()
		uMat = uMat.values()
		countStatus = countStatus.values()

		# Check the status of each result, we said
		# in case of simulation there might be cases
		# where the recovered answer was not correct
		# so we terminate the simulation by setting flag
		idx = 0
		for st in status:
			if st == False:
				flag = False
				break
			countTemp = countStatus[idx]
			# There was actually an answer so lets read it
			if countTemp > 0:
				count += countTemp
				# We updated aTotal and uTotal to solve later 
				#for secret key using these
				aTotal.extend(aMat[idx])
				uTotal.extend(uMat[idx])
			idx += 1
			# If we have enough equations we can break
			# and solve for secret key
			if count == n:
				break;

		# if flag is set to False we terminate and 
		# generate new systems
		measidx += 1
		if flag == False:
			break;

		# If we have enough equations we can break
		# and solve for secret key 
		if count == n:
			break;

	# If flag is not False then simulation is Okay
	# and we can try to solve for secret key
	if flag == True:
		aTotal = matrix(aTotal)
		uTotal = vector(uTotal)

		# Try to solve for secret key
		# we also catch the cases where there
		# might an error in the systems 
		try:
			KEY = aTotal.solve_right(uTotal)
		except:
			print("couldn't solve")

		# we can count how many of the secret key
		# coordinates we found are actually correct
		# hopefully the following should always be n
		correct = 0
		for i in range(len(KEY)):
			if KEY[i] == skList[i]:
				correct += 1

		# If they are equal just print the recovered key
		# The orginial key and 
		# The number of RLWE samples used 

		if list(KEY) == skList:
			print("\n\nCandidate Key")
			print(S(list(KEY)).list())
			print("\n\nChallenge Key")
			print(sk.list())
			print("\n\n** Correct Key Found **\n\n")
			print('Number of RLWE samples used : ' + str(measidx) + '\n\n')
			
			# print the time it took to run the attack
			Execution_Time = time.time() - Start_Time
			print("Time it took")
			print(Execution_Time)

			print("\n\n")
			print("** Done **")
	else:
		print("Wrong entries in the candidates")

	return


if __name__ == "__main__":
	main()
