import os
import sys
sys.path.append("pysrc/")
import field
import numpy as np
import subprocess
import time
import pickle
import random
import matplotlib.pyplot as plt
import scipy.optimize
np.set_printoptions(linewidth=1000)
import matplotlib
matplotlib.rcParams['axes.formatter.useoffset'] = False
ox = (0,0.13,0.28)

# DEPENDANCY - Need the FQHESphereJackGenerator of DiagHam to run this code.
PATH_TO_DIAGHAM = ""



def Stringefy(intlist,delim):
	"""Returns a string format of a list with elements separated by delim."""

	string = str(intlist)[1:-1].replace(', ',delim)
	return string


def relfac(n,m):
	"""Computes n!/m!"""

	out = 1.
	if n>m:
		for k in range(m+1,n+1): out *= k
	elif n<m:
		for k in range(n+1,m+1): out /= k
	return out


def BayesFit(xs,ys,fs,erys,pw=0):
	n = len(fs)
	Y = np.matrix([[y/e] for y,e in zip(ys,erys)])
	F = np.matrix([[f(x)/e for f in fs] for x,e in zip(xs,erys)])

	Sig = pw*pw*np.eye(n) + F.T*F
	errs = []
	for i in range(n):
		notI = [j for j in range(n) if j!=i]
		A = Sig[i,i]
		B = Sig[notI,i]
		C = np.delete(Sig,i,axis=0)
		C = np.delete(C,i,axis=1)
		p2 = A - B.T*np.linalg.inv(C)*B
		errs.append(1./np.sqrt(p2[0,0]))

	a = np.linalg.inv(Sig) * F.T * Y
	return [float(a[i,0]) for i in range(n)], errs


def fitWithReg(xs,ys,erys, fs, reg):
	X = np.matrix([[f(x)/e  for f in fs] for x,e in zip(xs,erys)])
	Y = np.matrix([[y/e] for y,e in zip(ys,erys)])
	I = np.eye(len(fs))
	
	F = np.linalg.inv(X.T*X + reg*I)*X.T*Y
	return F


def leastSquares(xs,ys,erys, fs, F):
	X = np.matrix([[f(x)/e  for f in fs] for x,e in zip(xs,erys)])
	Y = np.matrix([[y/e] for y,e in zip(ys,erys)])
	L = (X*F-Y).T*(X*F-Y)
	return L[0,0]
	


def simpleFit(xs,ys, d):
	extra = 1
	fs = [lambda x, k=n: x**k for n in range(d+extra)]
	xs,ys = np.array(xs), np.array(ys)
	erys = xs**(d+extra+1)
	precisions = min(erys)/erys
	F0 = fitWithReg(xs,ys,erys, fs[:d], 0)

	weight = 0.
	dev = np.zeros((d,1))
	for i in range(500):
		draw = [0]
		while sum(draw) < d+extra:
			draw = np.random.rand(len(xs)) < precisions
		F = fitWithReg(xs[draw],ys[draw],erys[draw], fs, 0.)
		w = 1./leastSquares(xs,ys,erys, fs, F)
		weight += w
		dev += w*np.square(F0-F[:d,0])
	dev = np.sqrt(dev/weight)
	return [float(F0[i,0]) for i in range(d)], [float(dev[i,0]) for i in range(d)]



class Jacks:

	def __init__(self, N, krs):
		self.N = N
		self.krs = krs


	def getRoots(self, dM):
		"""Outputs all the Bernevig Jack roots for k,r admissible Jacks of N-particle systems made up
		of particles with statistics s to which we add dM units of angular momentum."""
	
		k,r,s = self.krs
		if dM==0:
			if s=='f': unit, end = [1]*k+[0]*(r-k), [1]*(self.N%k)
			elif s=='b': unit, end = [k]+[0]*(r-1), [self.N%k]
			roots = [unit*(self.N/k)+end]
			while roots[0][-1]==0: roots[0].pop()
		else:
			roots, lower = [], self.getRoots(dM-1)
			for low in lower:
				low = low+[0]
				for i in range(len(low)-1):
					if (s=='f' and low[i]==1 and low[i+1]==0) or (s=='b' and low[i]>0):
						if sum(low[i+1:i+r+1]) < k:
							new = [l for l in low]
							new[i], new[i+1] = new[i]-1, new[i+1]+1
							if new not in roots: roots.append(new)
		return roots


	def convertRootToPartition(self, root):
		"""Takes in a root state in terms of occupations and returns angular momenta, i.e,
		110011 --> 5,4,1,0."""
	
		p = []
		for i in range(len(root)): p = [i]*root[i]+p
		return p


	def getFilename(self, root):
		try: os.makedirs("states/")
		except OSError: pass
		dM = sum(self.convertRootToPartition(root)) - sum(self.convertRootToPartition(self.getRoots(0)[0]))
		filename = "states/"+str(self.N)+"-"+",".join([str(a) for a in self.krs])
		filename += "-"+str(dM)+","+str(self.getRoots(dM).index(root))+".dat"
		return filename


	def getMatrixFilename(self, matrix, dM):
		try: os.makedirs("data/")
		except OSError: pass
		filename = "data/"+str(self.N)+"-"+",".join([str(a) for a in self.krs])
		filename += "-"+matrix+","+str(dM)+".dat"
		return filename


	def generateJacks(self, dMi, dMf):
		for dM in range(dMi, dMf+1):
			for root in self.getRoots(dM):
				self.generateSingleJack(root)


	def generateSingleJack(self, root):
		"""Generates the k,r admissible Jack with statistics s from root partition root."""
	
		k,r,s = self.krs
		outfile = self.getFilename(root)
		if os.path.isfile(outfile) and not os.path.isfile(outfile+'r'):
			pass
		else:
			if not os.path.exists("bin/"): os.makedirs("bin/")
			if not os.path.exists("states/"): os.makedirs("states/")
	
			if self.krs == (1,1,'f'):
				dM = sum(self.convertRootToPartition(root)) - sum(self.convertRootToPartition(self.getRoots(0)[0]))
				roots = sorted([self.convertRootToPartition(r) for r in self.getRoots(dM)], reverse=True)
				roots = [r for r in roots if r<=self.convertRootToPartition(root)]
				roots = [str(r) for r in roots]
				coeffs = ["0.0 " for r in roots]
				coeffs[0] = "1.0 "
	
				output = ""
				for c,r in zip(coeffs,roots):
					output += c+r+"\n"
				f = open(outfile+'r','w')
				f.write(output)
				f.close()
			else:
				reffile = 'bin/jack_ref'+str(time.time())+'.dat'
				f = open(reffile,'w')
				f.write("NbrParticles="+str(sum(root))+"\n")
				f.write("LzMax="+str(len(root)-1)+"\n")
				f.write("ReferenceState="+Stringefy(root,' ')+"\n")
				f.close()
		
				g = open(reffile,'r')
				print g.read()
				g.close()
		
				if s=='b': alpha = -(k+1.)/(r-1.)
				else: alpha = -(k+1.)/(r-k-1.)
				command = PATH_TO_DIAGHAM+"/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator"
				command += " --alpha "+str(alpha)+" --reference-file "+reffile
				if s=='f': command += " --fermion"
				command += " --txt-output "+outfile+'r'
				print "Generating "+outfile+"."
				subprocess.call(command,shell=True,stdout=open(os.devnull,'w'))
				os.remove(reffile)
			self.cleanJack(root)


	def cleanJack(self, root):
		"""The speed of the many of the routines depend on the Jack generator outputting the data in
		lexicographic order. Therefore, this check is necessary upon generation. It also cleans the
		Jack file up to makes the partitions into Fock states.
		The method is a piecewise comparison of the partitions between lines. The coefficients are
		corrected by the square root of a factorial for each integer of the partition and a term
		accounting for the number of terms each partition refers to."""
	
		f = open(self.getFilename(root)+'r','r')
		g = open(self.getFilename(root),'w')
		print "Checking the partitions are lexicographically ordered and cleaning."
		root0 = self.getRoots(0)[0]
		p0 = self.convertRootToPartition(root0)
		while True:
			v, part = self.interpretLine(f.readline())
			part += [0]*(self.N-len(part))
			if v == 'End': break
			else:
				v *= np.sqrt(self.norm(part,p0))
				g.write(str(v)+' '+str(part)+'\n')
		f.close()
		g.close()
		os.remove(self.getFilename(root)+'r')
		print "The Jack is clean."


	def interpretLine(self, line):
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


	def readLine(self, line):
		"""This reads the common line format of 'value' from files."""
	
		if line=='': return 'End'
		else: return float(line.split()[0])


	def norm(self, part,p0):
		"""Gives the normalisation of a given monomial, part (relative to some base monomial, p0)"""
	
		Zpart = np.prod([relfac(n,n0) for n,n0 in zip(part,p0)])
		Zpart /= np.prod([np.math.factorial(part.count(n)) for n in set(part)])
		return Zpart


	def applySingleInteraction(self, interaction, root):
		"""Apply the interaction to the given root. The interaction name must have some definition in
		the C++ code which it knows how to interpret in terms of pseudopotentials."""

		k,r,s = self.krs
		if not os.path.exists(self.getFilename(root)+interaction):
			command = "./Main "+str(self.N)+" "+",".join([str(k),str(r),s])
			command += " "+"".join([str(n) for n in root])+" "+interaction
			command += " "+self.getFilename(root)
			print command
			subprocess.call(command,shell=True)#,stdout=open(os.devnull,'w'))


	def applyInteractions(self, interaction, dMi, dMf):
		for dM in range(dMi, dMf+1):
			for root in self.getRoots(dM):
				self.applySingleInteraction(interaction, root)


	def singleParticleInner(self, confinement,left,right):
		"""Takes an inner product with some weight for each orbital. Effectively the
		matrix element of a single-particle operator diagonal in the orbital basis"""
	
		matrix = confinement[0]
		Uweights = confinement[1]

		root0 = self.convertRootToPartition(self.getRoots(0)[0])
		dM = sum(self.convertRootToPartition(left)) - sum(root0)
		roots = self.getRoots(dM)
		i, j = roots.index(left), roots.index(right)
		U = self.loadMatrix(matrix, dM)
	
		if U[i,j] == -1:
			ll = open(self.getFilename(left), 'r')
			lr = open(self.getFilename(right), 'r')
			vl,pl = self.interpretLine(ll.readline())
			vr,pr = self.interpretLine(lr.readline())
			Uij = 0.
			while vl!="End" or vr!="End":
				if pl==pr:
					Uij += vl*vr * sum([Uweights[pi] for pi in pl])
					vl,pl = self.interpretLine(ll.readline())
					vr,pr = self.interpretLine(lr.readline())
				elif pr<pl:
					vl,pl = self.interpretLine(ll.readline())
				elif pl<pr:
					vr,pr = self.interpretLine(lr.readline())
			ll.close()
			lr.close()
	
			print(Uij)
	
			U = self.loadMatrix(matrix,dM)
			U[i,j], U[j,i] = Uij, Uij
			f = open(self.getMatrixFilename(matrix, dM),'w')
			pickle.dump(U, f)
			f.close()


	def singleInnerProduct(self, matrix,left,right):
		"""Usual inner product between states left and right defined by the matrix."""
	
		if matrix=="overlap": ap=""
		else: ap=matrix
	
		root0 = self.convertRootToPartition(self.getRoots(0)[0])
		dM = sum(self.convertRootToPartition(left)) - sum(root0)
		roots = self.getRoots(dM)
		i, j = roots.index(left), roots.index(right)
		V = self.loadMatrix(matrix, dM)
	
		if V[i,j] == -1:
			ll = open(self.getFilename(left), 'r')
			fl = open(self.getFilename(left)+ap, 'r')
			lr = open(self.getFilename(right), 'r')
			fr = open(self.getFilename(right)+ap, 'r')
			vl,pl = self.interpretLine(ll.readline())
			vr,pr = self.interpretLine(lr.readline())
			cl = self.readLine(fl.readline())
			cr = self.readLine(fr.readline())
			Vij = 0.
			while vl!="End" or vr!="End":
				if pl==pr:
					Vij += vl*cr + cl*vr
					vl,pl = self.interpretLine(ll.readline())
					vr,pr = self.interpretLine(lr.readline())
					cl = self.readLine(fl.readline())
					cr = self.readLine(fr.readline())
				elif pr<pl:
					vl,pl = self.interpretLine(ll.readline())
					cl = self.readLine(fl.readline())
				elif pl<pr:
					vr,pr = self.interpretLine(lr.readline())
					cr = self.readLine(fr.readline())
			ll.close()
			fl.close()
			lr.close()
			fr.close()
	
			print(Vij)
	
			if matrix=="overlap": Vij /= 2.
			V = self.loadMatrix(matrix,dM)
			V[i,j], V[j,i] = Vij, Vij
			f = open(self.getMatrixFilename(matrix, dM),'w')
			pickle.dump(V, f)
			f.close()


	def loadMatrix(self, interaction, dM):
		if not os.path.exists("data/"): os.makedirs("data/")
		filename = self.getMatrixFilename(interaction, dM)
		if os.path.exists(filename):
			M = pickle.load(open(filename,'r'))
		else:
			roots = self.getRoots(dM)
			dim = len(roots)
			M = np.matrix(-np.ones((dim,dim)))
			pickle.dump(M, open(filename,'w'))
		return M


	def innerProducts(self, matrix, dMi, dMf):
		for dM in range(dMi, dMf+1):
			roots = self.getRoots(dM)
			for i in range(len(roots)):
				for j in range(i+1):
					self.singleInnerProduct(matrix, roots[i],roots[j])


	def confinement(self, confiningData, dMi, dMf):
		"""The 'confiningData' must be fed in as a list. confiningData[0] should
		be a name for the matrix and confiningData[2] is a list of all the weights
		for each orbital (and make sure you feed in enough to reach the highest angular
		momentum in your system)."""

		for dM in range(dMi, dMf+1):
			roots = self.getRoots(dM)
			for i in range(len(roots)):
				for j in range(i+1):
					self.singleParticleInner(confiningData, roots[i],roots[j])


	def getSpectrum(self, interaction, dMi, dMf, interact=True):
		"""Returns the numerically calculated spectrum from data. Calculate these by calling
		innerProducts() for you particular interaction and for the 'overlap' matrix. Make sure
		you call applyInteractions() before you try to calculate an interaction matrix. NB:
		the 'interact=True' flag no longer does anything."""

		xs, ys = [], []
		for dM in range(dMi, dMf+1):
			O = self.loadMatrix("overlap",dM)
			V = self.loadMatrix(interaction,dM)
			H = np.linalg.inv(O)*V
			eigs = list(np.linalg.eigvals(H))
			xs += [dM]*len(eigs)
			ys += eigs
		return xs, ys


	def interactingCoeffs(self, interaction):
		"""Finds the coefficients for the particular interaction."""

		if self.krs[0]==1:
			xs, ys = self.getSpectrum(interaction,3,3)
			y0 = max(ys)
			ys = sorted([y-y0 for y in ys])
			g22 = ys[1]/-12.
			T33size = ys[0]+48*g22
			if abs(T33size/y0) < 1e-10: g33 = 0
			else: g33 = (ys[0]+48*g22)/240.
			return y0,g22,g33
		elif self.krs[0]==2:
			xs,ys = self.getSpectrum(interaction,0,0)
			e0 = ys[0]
			xs,ys = self.getSpectrum(interaction,3,3)
			ys = sorted(ys)
			g22 = (ys[-2]-e0)/-12.
			g01 = (ys[0]-5*ys[1]+4*e0)/14.
			#TODO Write routine to fit this coefficient for general interaction
			g12 = 0.0021/self.N**1.5
			return e0,g01,g22,g12
			

	def interactingEffective(self, coeffs, dM):
		"""Returns the effective Hamiltonian at angular momentum dM given the coeffs."""

		if self.krs[0]==1:
			NrB = self.N*np.sqrt(self.krs[1])
			T22 = field.term_matrix([2,2],[],dM,True)
			T22 += 1./(2.*NrB) * field.term_matrix([2,2,2],[],dM,True)
			T22 += 5./(16.*NrB**2) * field.term_matrix([2,2,2,2],[],dM,True)
			T22 += 7./(32.*NrB**3) * field.term_matrix([2,2,2,2,2],[],dM,True)
			T22 += 21./(128.*NrB**4) * field.term_matrix([2,2,2,2,2,2],[],dM,True)
			
			T33 = field.term_matrix([3,3],[],dM,True)
			T33 += 5./(2*NrB) * field.term_matrix([3,3,2],[],dM,True)
			T33 += -15./(2*NrB) * field.term_matrix([2,2,2],[],dM,True)
			T33 += 35./(8*NrB**2) * field.term_matrix([3,3,2,2],[],dM,True)
			T33 += -125./(8*NrB**2) * field.term_matrix([2,2,2,2],[],dM,True)
			T33 += 105./(16*NrB**3) * field.term_matrix([3,3,2,2,2],[],dM,True)
			T33 += -791./(32*NrB**3) * field.term_matrix([2,2,2,2,2],[],dM,True)

			e0, g22, g33 = coeffs
			H = e0*np.eye(len(T22)) + g22*T22 + g33*T33
			return H
		elif self.krs[0]==2:
			NrB = self.N*np.sqrt(self.krs[1]/2)
			T01 = field.term_matrix([],[0,1],dM)
			T01 += 1./(2*NrB)*field.term_matrix([2],[0,1],dM)
			T01 += 3./(8*NrB**2)*field.term_matrix([2,2],[0,1],dM)
			T01 += 5./(16*NrB**3)*field.term_matrix([2,2,2],[0,1],dM)
			T01 += 35./(128*NrB**4)*field.term_matrix([2,2,2,2],[0,1],dM)

			T22 = field.term_matrix([2,2],[],dM)
			T22 += 1./(2.*NrB) * field.term_matrix([2,2,2],[],dM)
			T22 += 5./(16.*NrB**2) * field.term_matrix([2,2,2,2],[],dM)
			T22 += 7./(32.*NrB**3) * field.term_matrix([2,2,2,2,2],[],dM)

			T12 = field.term_matrix([],[1,2],dM)
			T12 += 3./(2*NrB)*field.term_matrix([2],[1,2],dM)
			T12 += 15./(8*NrB**2)*field.term_matrix([2,2],[1,2],dM)
			T12 += 35./(16*NrB**3)*field.term_matrix([2,2,2],[1,2],dM)

			e0, h01, h22, h12 = coeffs
			H = e0*np.eye(len(T01)) + h01*T01 + h22*T22 + h12*T12
			return H


	def confiningCoeffs(self, matrix):
		"""Finds the coefficients for a particular confining potential"""

		if self.krs[0]==1:
			xs, ys = self.getSpectrum(matrix,0,0, False)
			e0 = ys[0]

			init = [0.,0.]
			for dM in [4,3,5]:
				xs, ys = self.getSpectrum(matrix,dM,dM, False)
				def loss(gs):
					coeffs = (e0,)+tuple(gs)
					H = self.confiningEffective(coeffs,dM)
					es = list(np.linalg.eigvals(H))
					return sum([(e-y)**2 for e,y in zip(sorted(es),sorted(ys))])
				res = scipy.optimize.minimize(loss, init)
				init = res.x
			coeffs = (e0,)+tuple(res.x)
			return coeffs
		elif self.krs[0]==2:
			xs, ys = self.getSpectrum(matrix,0,0, False)
			e0 = ys[0]

			init = [-0.1850,0.1568,0.0106,-0.0386]
			kmax = 30
			for k in range(kmax+1):
				xs, ys = self.getSpectrum(matrix,4,4, False)
				ys = sorted(ys)
				d = len(ys)
				if k==kmax: samples = range(d)
				elif k%2==0: samples = [random.randint(0,d-1) for i in range(k+1)]
				else: samples = [random.randint(0,d-1) for i in range(4)]
				def loss(gs):
					coeffs = (e0,)+tuple(gs)
					H = self.confiningEffective(coeffs,4)
					es = sorted(list(np.linalg.eigvals(H)))
					return sum([(es[s]-ys[s])**2 for s in samples])
				res = scipy.optimize.minimize(loss, init)
				init = res.x
			coeffs = (e0,)+tuple(res.x)
			return coeffs


	def confiningEffective(self, coeffs, dM):
		"""Returns the effective Hamiltonian at angular momentum dM given the coeffs."""

		if self.krs[0] == 1:
			e0,v,g = coeffs
			H = v*field.term_matrix([1,1],[], dM, True)
			H += g*field.term_matrix([1,1,1],[],dM, True)
			H += e0*np.eye(len(H))
		elif self.krs[0] == 2:
			e0,v1,v2,g1,g2 = coeffs
			H = v1*field.term_matrix([],[0,1], dM, False)
			H += v2*field.term_matrix([1,1],[], dM, False)
			H += g1*field.term_matrix([1,1,1],[], dM, False)
			H += g2*field.term_matrix([1],[0,1], dM, False)
			H += e0*np.eye(len(H))
		return H





def coefficientScalings(Ns, krs, interaction, saveAs=""):
	"""Fits the coefficients for the particular interaction and finds the scaling behaviour."""

	# Find and fit the data
	xs,ys,zs = [],[],[]
	for N in Ns:
		J = Jacks(N, krs)
		if krs[0]==1:
			e0,g22,g33 = J.interactingCoeffs(interaction)
			xs.append(1./np.sqrt(N))
			ys.append(g22*N**1.5)
			zs.append(g33*N**2.5)
		elif krs[0]==2:
			e0,g01,g22,g12 = J.interactingCoeffs(interaction)
			xs.append(1./np.sqrt(N))
			ys.append(g01*N**0.5)
			zs.append(g22*N**1.5)

	if krs==(1,1,'f'): d = 5
	else: d = 1
	txs = [max(xs)*i/100. for i in range(101)]
	cs, es = simpleFit(xs,ys, d)
	tys = [sum([cs[i]*x**i for i in range(d)]) for x in txs]
	cs2, es2 = simpleFit(xs,zs, d)
	tzs = [sum([cs2[i]*x**i for i in range(d)]) for x in txs]


	# Print the tabulated data for pseudopotentials
	if krs==(1,1,'f'):
		out = " & ".join([r"${0:.3f}\pm{1:.3f}$".format(cs[i],es[i]+0.0005) for i in range(3)])
		out = " & ".join([out, r"${0:.3f}\pm{1:.3f}$".format(cs2[0],es2[0]+0.0005)])
		print r"$V_" + interaction[6:] + "$ & " + out + r" \\ \hline"
	else:
		print r"$V_" + interaction[6:] + r"$ & ${0:.3f}\pm{1:.3f}$".format(cs[0],es[0]) + r" \\ \hline"


	# Plot data
	if "pseudo" in interaction:
		scale = "v_"+interaction[6:]
		if interaction=="pseudo1": ints = r"1$^{st}$ pseudopential"
		elif interaction=="pseudo2": ints = r"2$^{nd}$ pseudopential"
		elif interaction=="pseudo3": ints = r"3$^{rd}$ pseudopential"
		else: ints = interaction[6:]+r"$^{th}$ pseudopential"
		ints = "the "+ints
	elif interaction=="exp":
		scale = "w_0"
		ints = "exponential repulsion"

	fig = plt.figure()
	ax = fig.add_axes([0.17,0.13,0.8,0.8])
	ax2 = fig.add_axes([0.5,0.28,0.4,0.3])
	
	ax.plot(txs,tys,color="orange",linestyle="dashed",zorder=-100, linewidth=2)
	ax.plot([0,max(xs)],[0,0], color="k", linewidth=0.5)
	ax.scatter(xs,ys,color=ox, s=50)
	ax.plot(xs,ys,color=ox, linewidth=2)
	ax.set_xlim(0,max(xs))
	ax.set_ylim(min(ys+[0])*1.2,max(ys)*1.2)
	ax.set_xlabel(r"$\sqrt{N}^{-1}$",fontsize=30)

	if krs==(1,1,'f') or krs==(2,2,'b'): ax2.plot(txs,tzs,color="orange",linestyle="dashed",zorder=-100, linewidth=2)
	ax2.plot([0,max(xs)],[0,0], color="k", linewidth=0.5)
	ax2.scatter(xs,zs,color=ox, s=50)
	ax2.plot(xs,zs,color=ox, linewidth=2)
	ax2.set_xlim(0,max(xs))
	ax2.set_ylim(min(zs+[0])*1.2,max(zs)*1.2)
	ax2.set_xlabel(r"$\sqrt{N}^{-1}$",fontsize=30)
	if krs[0]==1:
		ax.set_ylabel(r"$\frac{h_{22}}{"+scale+r"\sqrt{N}^{-3}}$",fontsize=30)
		ax2.set_ylabel(r"$\frac{h_{33}}{"+scale+r"\sqrt{N}^{-5}}$",fontsize=30)
	elif krs[0]==2:
		ax.set_ylabel(r"$\frac{h_{\emptyset,01}}{"+scale+r"\sqrt{N}}$",fontsize=30)
		ax2.set_ylabel(r"$\frac{h_{22,\emptyset}}{"+scale+r"\sqrt{N}^{-3}}$",fontsize=30)
		

	title = r"Coefficient scaling for "+ints+r" at filling $\nu=1"
	if krs[0]!=krs[1]: title += r"/"+str(krs[1]/krs[0])
	title += r"$"
	ax.set_title(title, fontsize=16)
	if saveAs=="": plt.show()
	else: plt.savefig(saveAs+".pdf")
	plt.clf()


def compareConfinement(N, matrix, krs, upto=5, saveAs=""):
	"""This function does the same as compareInteraction except that the matrix should
	correspond to some confinement potential."""

	# Finding the data
	J = Jacks(N, krs)
	coeffs = J.confiningCoeffs(matrix)
	if krs[0]==1: e0,v,g = coeffs
	elif krs[0]==2: e0,v1,v2,g,t = coeffs

	xs,ys = J.getSpectrum(matrix,0,upto, False)
	zs = []
	for dM in range(upto+1):
		H = J.confiningEffective(coeffs,dM)
		zs += list(np.linalg.eigvals(H))


	#Plotting the data
	yrange = max(ys+zs)-min(ys+zs)
	y0 = min(ys)

	fig = plt.figure()
	ax = fig.add_axes([0.15,0.13,0.8,0.8])
	ax.set_ylim(y0-yrange/20., y0+yrange+yrange/20.)
	ax.set_xlim(-0.4, max(xs)+0.4)

	ax.scatter(xs,ys, color=ox,zorder=100, s=50)
	for x,z in zip(xs,zs):
		ax.plot([x-0.2,x+0.2],[z,z],'orange', linewidth=2)

	y0 = min(ys)
	yrange = max(ys)-y0

	ax.scatter([0.25],[y0+0.92*yrange],color=ox)
	ax.annotate(r"are $\delta\mathcal{H}$", (0.5,y0+0.9*yrange), fontsize=20)
	ax.plot([0,0.4],[y0+0.82*yrange,y0+0.82*yrange],color="orange")
	ax.annotate(r"are $H$", (0.5,y0+0.8*yrange), fontsize=20)

	m = krs[1]/krs[0]
	if m==1: nu = "1"
	else: nu = "1/"+str(m)
	if krs[0]==1:
		ax.annotate(r"$v="+"{:.4f}".format(v)+r"U_0$", (0,y0+0.7*yrange), fontsize=20)
		ax.annotate(r"$g="+"{:.4f}".format(g)+r"U_0$", (0,y0+0.6*yrange), fontsize=20)
		ax.set_title(r"$\nu="+nu+"$ Laughlin edge spectrum given "+matrix+" confinement", fontsize=18)
	if krs[0]==2:
		ax.annotate(r"$v_1="+"{:.4f}".format(-v1*2)+r"U_0$", (0,y0+0.7*yrange), fontsize=20)
		ax.annotate(r"$v_2="+"{:.4f}".format(v2*2)+r"U_0$", (0,y0+0.6*yrange), fontsize=20)
		ax.annotate(r"$g_1="+"{:.4f}".format(g)+r"U_0$", (0,y0+0.5*yrange), fontsize=20)
		ax.annotate(r"$g_2="+"{:.4f}".format(t)+r"U_0$", (0,y0+0.4*yrange), fontsize=20)
		ax.set_title(r"$\nu="+nu+"$ Moore-Read edge spectrum given "+matrix+" confinement", fontsize=18)
	ax.set_ylabel(r"$E/U_0$", fontsize=30)
	ax.set_xlabel(r"$\Delta L$", fontsize=30)
	if saveAs=="": plt.show()
	else: plt.savefig(saveAs+".pdf")


def compareInteraction(N, matrix, krs, upto=5, simple=False, saveAs=""):
	"""This function fits the coefficients for the particular matrix, which must be stored as
	data. It then finds the numerical and theoretical spectra. It then plots the data, finding
	a suitable harmonic confinement to make the plot look nice and making annotations."""

	# Finding the data
	J = Jacks(N, krs)
	coeffs = J.interactingCoeffs(matrix)
	if simple==True: coeffs = tuple([coeffs[i] if i<2 else 0. for i in range(len(coeffs))])
	if krs[0]==1: e0,g22,g33 = coeffs
	elif krs[0]==2: e0,g01,g22,g12 = coeffs

	xs,ys = J.getSpectrum(matrix,0,upto, True)
	zs = []
	for dM in range(upto+1):
		H = J.interactingEffective(coeffs,dM)
		zs += list(np.linalg.eigvals(H))


	# Plotting the data
	yrange = max(ys+zs)-min(ys+zs)
	v = 5./3*yrange/max(xs)
	ys = [y+v*x for x,y in zip(xs,ys)]
	zs = [z+v*x for x,z in zip(xs,zs)]
	y0 = min(ys)
	yrange = max(ys)-y0

	fig = plt.figure()
	ax = fig.add_axes([0.15,0.13,0.8,0.8])
	ax.set_ylim(y0-yrange/20., y0+yrange+yrange/20.)
	ax.set_xlim(-0.4, max(xs)+0.4)

	ax.plot([-1,max(xs)+1],[e0-v, e0+v*(max(xs)+1)], color=ox, zorder=-100)
	ax.scatter(xs,ys, color=ox,zorder=100, s=50)
	for x,z in zip(xs,zs):
		ax.plot([x-0.2,x+0.2],[z,z],'orange', linewidth=2)

	if matrix=="exp": pseudo, scale = "exponential repulsion", "w_0"
	elif matrix=="pseudo2": pseudo, scale = r"a 2$^{nd}$ pseudopotential", "v_2"
	elif matrix=="pseudo3": pseudo, scale = r"a 3$^{rd}$ pseudopotential", "v_3"
	else: pseudo, scale = "a "+matrix[-1]+r"$^{th}$ pseudopotential", "v_"+matrix[-1]

	ax.scatter([0.0],[y0+0.97*yrange],color=ox, s=50)
	ax.annotate(r"are $\delta\mathcal{H}$", (0.3,y0+0.95*yrange), fontsize=20)
	ax.plot([-0.2,0.2],[y0+0.87*yrange,y0+0.87*yrange],color="orange", linewidth=2)
	ax.annotate(r"are $H$", (0.3,y0+0.85*yrange), fontsize=20)
	if krs[0]==1:
		ax.annotate(r"$h_{22}=\frac{"+"{:.5f}".format(g22*N**1.5)+r"}{\sqrt{N}^3}"+scale+"$", (-0.2,y0+0.7*yrange), fontsize=20)
		if simple==False:
			ax.annotate(r"$h_{33}=\frac{"+"{:.5f}".format(g33*N**2.5)+r"}{\sqrt{N}^5}"+scale+"$", (-0.2,y0+0.55*yrange), fontsize=20)
	elif krs[0]==2:
		ax.annotate(r"$h_{\emptyset,01}=\frac{"+"{:.4f}".format(g01*N**0.5)+r"}{\sqrt{N}}"+scale+"$", (-0.2,y0+0.7*yrange), fontsize=20)
		if simple==False:
			ax.annotate(r"$h_{22,\emptyset}=\frac{"+"{:.4f}".format(g22*N**1.5)+r"}{\sqrt{N}^3}"+scale+"$", (-0.2,y0+0.55*yrange), fontsize=20)
			ax.annotate(r"$h_{\emptyset,12}=\frac{"+"{:.4f}".format(g12*N**1.5)+r"}{\sqrt{N}^3}"+scale+"$", (-0.2,y0+0.4*yrange), fontsize=20)

	m = krs[1]/krs[0]
	if m==1: nu="1"
	else: nu="1/"+str(m)
	ax.set_title(r"$\nu="+nu+r"$ and $N="+str(N)+r"$ with "+pseudo, fontsize=18)
	ax.set_ylabel(r"$E/"+scale+"$", fontsize=30)
	ax.set_xlabel(r"$\Delta L$", fontsize=30)
	ax.annotate(r"$U="+"{0:.4f}".format(v*N*krs[1])+scale+r"\left(\frac{r}{R}\right)^2$", (max(xs)/2.,y0), fontsize=25, ha="center")
	if saveAs=="": plt.show()
	else: plt.savefig(saveAs+".pdf")

	




if __name__=="__main__":
	pass

#	N, krs = 8, (1,2,'b')
#	compareInteraction(N, "pseudo2", krs, upto=5, simple=True, saveAs="half_pseudo")
#	N, krs = 8, (1,3,'f')
#	compareInteraction(N, "pseudo3", krs, upto=5, simple=True, saveAs="third_pseudo")
#	N, krs = 8, (1,4,'b')
#	compareInteraction(N, "pseudo4", krs, upto=5, simple=True, saveAs="quarter_pseudo")
#
#	N, krs = 100, (1,1,'f')
#	compareInteraction(N, "exp", krs, upto=8, saveAs="integer")
#	N, krs = 12, (1,2,'b')
#	compareInteraction(N, "exp", krs, upto=6, saveAs="halfexp")
#	N, krs = 12, (1,3,'f')
#	compareInteraction(N, "exp", krs, upto=6, saveAs="thirdexp")
#	N, krs = 16, (2,2,'b')
#	compareInteraction(N, "exp", krs, upto=5, saveAs="winning")
#
#	N, krs = 10, (1,3,'f')
#	compareConfinement(N, "octic", krs, upto=5, saveAs="octic")
#	N, krs = 14, (2,2,'b')
#	compareConfinement(N, "octic", krs, upto=5, saveAs="MR-octic")

#	Ns, krs = range(3,21)+[25+5*i for i in range(28)], (1,1,'f')
#	coefficientScalings(Ns, krs, "pseudo1", saveAs="integerCoeffs")
#	Ns, krs = range(10,21)+[25+5*i for i in range(28)], (1,1,'f')
#	coefficientScalings(Ns, krs, "pseudo3", saveAs=" ")
#	Ns, krs = range(10,21)+[25+5*i for i in range(28)], (1,1,'f')
#	coefficientScalings(Ns, krs, "pseudo5", saveAs=" ")
#	Ns, krs = range(10,21)+[25+5*i for i in range(28)], (1,1,'f')
#	coefficientScalings(Ns, krs, "pseudo7", saveAs=" ")
#	Ns, krs = range(15,21)+[25+5*i for i in range(28)], (1,1,'f')
#	coefficientScalings(Ns, krs, "pseudo9", saveAs=" ")

#	Ns, krs = range(5,13), (1,2,'b')
#	coefficientScalings(Ns, krs, "pseudo2", saveAs="half_scaling")
#	Ns, krs = range(5,13), (1,2,'b')
#	coefficientScalings(Ns, krs, "pseudo4", saveAs=" ")
#	Ns, krs = range(5,13), (1,2,'b')
#	coefficientScalings(Ns, krs, "pseudo6", saveAs=" ")

#	Ns, krs = [6,8,10,12,14,16,18], (2,2,'b')
#	coefficientScalings(Ns, krs, "exp", saveAs="MR_scaling")
		
	
	
	
	

