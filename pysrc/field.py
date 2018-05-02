#!../../../../usr/bin/env python
import os
import numpy as np
import pickle
import itertools

#################################
### GENERAL UTILITY FUNCTIONS ###
#################################

def Stringefy(intlist,delim):
	"""Returns a string format of a list with elements separated by delim."""

	string = str(intlist)[1:-1].replace(', ',delim)
	return string



def InterpretLine(line):
	"""This reads the common line format of 'value partition' from Jack files."""

	if line == '':
		return 'End',[]
	else:
		data = line.split(' [')
		value = float(data[0])
		try: partition = [int(n) for n in data[1][:-2].split(',')]
		except ValueError:
			print line
			raise Exception
		return value, partition




def relfac(n,m):
	"""Computes n!/m!"""

	out = 1.
	if n>m:
		for k in range(m+1,n+1): out *= k
	elif n<m:
		for k in range(n+1,m+1): out /= k
	return out



def norm(part,p0):
	"""Gives the normalisation of a given monomial, part (relative to some base monomial, p0)"""

	Zpart = np.prod([relfac(n,n0) for n,n0 in zip(part,p0)])
	Zpart /= np.prod([np.math.factorial(part.count(n)) for n in set(part)])
	return Zpart




def Partition(N, max_n=0):
	"""Creates all the integer partitions of the integer N including only integers less than or
	equal to max_n (unless max_n=0, indicating no restrictions).
	The method is to take out some integer, n, as the maximum integer in the partition, then use
	the function recursively to partition the remainder and append the results. The outputs come
	out lexicographically ordered using this method."""

	parts = []
	if(max_n==0 or max_n>N): max_n=N
	if(N==0):
		parts.append([])
	elif(N==1):
		parts.append([1])
	else:
		for n in reversed(range(1,max_n+1)):
			sub_parts = Partition(N-n, n)
			for sub in sub_parts:
				parts.append([n]+sub)
	return parts



def list_diff(p1,p2):
	"""Returns the elements in p1 that aren't in p2 and vice versa"""

	left, right = [1*l for l in p1], [1*r for r in p2]
	for l in reversed(left):
		if l in right:
			left.remove(l)
			right.remove(l)
	return left,right



def signed_sort(lst):
	"""Sorts a list and returns the sign of the permutation required to sort"""

	sign = 1
	for passes in range(len(lst)-1,0,-1):
		for i in range(passes):
			if lst[i]>lst[i+1]:
				lst[i],lst[i+1] = lst[i+1],lst[i]
				sign *= -1
	return lst,sign



def fit(xs,ys,fs,lam=0):
	"""A simple linear fit given functions fs with a weight lam for large coefficients"""

	n = len(fs)
	yf = np.matrix(np.zeros((n,1)))
	ff = np.matrix(np.zeros((n,n)))
	for i in range(n):
		yf[i,0] = sum([y*fs[i](x) for x,y in zip(xs,ys)])
		for j in range(i+1):
			ff[i,j] = sum([fs[j](x)*fs[i](x) for x in xs])
			if i==j: ff[i,j] += lam
			else: ff[j,i] = ff[i,j]
			
	a = np.linalg.inv(ff)*yf
	return [float(a[i,0]) for i in range(n)]



def teXify(strng):
	"""Turns a sting of operators into TeX"""

	strng = strng.replace("B","P_")
	strng = strng.replace("S","\psi_")
	strng = strng.replace("F","\psi_")
	while "." in strng:
		i = strng.index(".")
		top = int(strng[i-1])*2+1
		strng = strng.replace(strng[i-1:i+2],"{"+str(top)+"/2}")
	return strng




##################################
### THEORY CALCULATION METHODS ###
##################################

class state:
	"""A class for a general state of modes acting on a vacuum"""

	def __init__(self, bosons,fermions,factor=None):
		self.bosons = sorted([int(b) for b in bosons],reverse=True)
		self.fermions = sorted([float(f) for f in fermions],reverse=True)
		self.calc_norm()
		if factor==None: self.factor = 1./np.sqrt(self.norm)
		else: self.factor = factor

	def calc_norm(self):
		self.norm = 1.
		for b in set(self.bosons):
			mb = self.bosons.count(b)
			self.norm *= b**mb * np.math.factorial(mb)

	def __str__(self):
		bosons = ["B"+str(b) for b in self.bosons]
		fermions = ["F"+str(f) for f in self.fermions]
		return "["+str(self.factor)+"] "+" ".join(bosons+fermions)+" |0>"

	def __repr__(self):
		return str(self)

	def inner(self, other):
		if other.bosons==self.bosons and other.fermions==self.fermions:
			return self.factor*self.norm*other.factor
		else: return 0.

	def annihilate(self,what):
		if isinstance(what,float):
			if what %1. == 0.: boson,mode = True,int(what)
			else: boson,mode = False,what
		elif isinstance(what,int): boson,mode = True,what
		elif isinstance(what,str):
			if what[0]=="B": boson,mode = True,int(what[1:])
			elif what[0] in ["S","F"]: boson,mode = False,float(what[1:])

		if boson:
			if mode in self.bosons:
				new_factor = self.factor*mode*self.bosons.count(mode)
				i = self.bosons.index(mode)
				new_bosons = self.bosons[:i]+self.bosons[i+1:]
				return state(new_bosons,self.fermions,new_factor)
			else: return state([],[],0.)
		else:
			if mode in self.fermions:
				new_factor = self.factor*(-1)**self.fermions.index(mode)
				i = self.fermions.index(mode)
				new_fermions = self.fermions[:i]+self.fermions[i+1:]
				return state(self.bosons,new_fermions,new_factor)
			else: return state([],[],0.)



def term_matrix(blab,flab,dM, fermionsOff=False):
	"""Calculates the matrix for the term labelled by bosonic and fermionic indices at angular momentum dM"""

	if not os.path.exists("mats/"): os.makedirs("mats/")

	filename = "mats/"+str(dM)+str(blab)+str(flab)+".dat"
	if fermionsOff: filename += "b"
	if os.path.isfile(filename):
		return pickle.load(open(filename,'r'))

	if fermionsOff: states = gen_states(dM, True)
	else: states = gen_states(dM)
	print states
	d = len(states)
	M = np.matrix(np.zeros((d,d))).astype(float)
	nb,nf = len(blab),len(flab)
	bosons = range(-dM,0)+range(1,dM+1)
	fermions = [n+0.5 for n in range(-dM,dM)]

	for s in itertools.product(*[bosons]*nb+[fermions]*nf):
		if sum(s)==0:
			### s is the string of operators (note they are in reverse order, so a1a_{-1} will be given by [-1,1])
			s = list(s)

			### Calculate the factors within the sum (recalling the reversed order)
			m = 1.
			for b,i in zip(reversed(blab),range(nb)):
				for bi in range(1,b): m *= s[i]-bi
			for f,i in zip(reversed(flab),range(nb,nb+nf)):
				for fi in range(f): m *= s[i]-fi-0.5

			### Normal order (in such a way that the first list item is the first operator applied, i.e, reverse order)
			s[:nb] = sorted(s[:nb])
			s[nb:],sign = signed_sort(s[nb:])
			m *= sign

			### Apply operators (being careful to apply in the correct order)
			if m!=0.:
				for i in range(d):
					for j in range(d):
						left, right = states[i],states[j]
						for si in [sj for sj in s if sj<0]:
							right = right.annihilate(-si)
						for si in reversed([sj for sj in s if sj>0]):
							left = left.annihilate(si)
						mij = m*left.inner(right)
						M[i,j] += mij
	pickle.dump(M, open(filename,'w'))
	return M



def gen_states(dM, fermionsOff=False):
	"""Generates all the states at angular momentum up to dM (MAX 2 FERMION MODES - WOULD NOT BE COMPLETE FOR dM>7)"""

	states = [state(p,[]) for p in Partition(dM)]
	if not fermionsOff:
		for dMf in range(2,dM+1):
			bosons = Partition(dM-dMf)
			for f1 in range(dMf/2):
				fpart = [f1+0.5,dMf-f1-0.5]
				for b in bosons:
					states.append(state(b,fpart))
	return states


def get_norms(dM,N_even=True):
	"""Finds the normalisation of all states"""

	states = gen_states(dM)
	norms = [1./s.factor**2 for s in states]
	return norms



