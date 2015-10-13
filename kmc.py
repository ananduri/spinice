import numpy as np
from numpy.linalg import norm
from scipy.special import erfc
# import pdb
import random
import Queue
from scipy.io import loadmat,savemat
import timeit

# pre compute the J matrix for large N

# modify statevector when flipping spins (don't)


# remove all Lattice.magmoms
# need nnguide2? nope


# make this a Python learning experience
# use 'bool(random.getrandbits(1))' 
# why is realenergy so small?

# need to modify position to take into account PBCs

# 16*a^3 spins

def main():
	a = 4
	lat = Lattice(a)

	d=loadmat('J22_a4')
	J=d['J']

	start = timeit.default_timer()

	#J=lat.buildintmat(2,2)

	alist=[]
	for x in xrange(100):
		lat=Lattice(a)
		w=lat.getautocorr(J,3.,50)
		alist.append(w)

	alist=np.array(alist)
	savemat('autocorr_a4_T3',{'alist':alist})

	end = timeit.default_timer()
	print end-start

	#savemat('J22_a6',{'J':J})

		
	#alist=[]
	#for x in xrange(100):
	#	lat=Lattice(a)
	#	w=lat.getautocorr(J,1.5,100)
	#	alist.append(w)

	#alist=np.array(alist)
	#savemat('autocorr100_a4_T1.5',{'alist':alist})
	


	# lat.randinit()
	# e = lat.spinenergy((0,0,0))
	# print 'Tenergy: {0:>7}\n'.format(e[0])
	# print 'renergy: {0:>7}\n'.format(e[1])
	# print 'kenergy: {0:>7}\n'.format(e[2])
	# print 'senergy: {0:>7}'.format(e[1]+e[2])
	# print 'dotenergy: {0:>7}'.format(e[3])
	# print 'nnenergy: {0:>7}'.format(e[4])
	# print lat.brute([0,0,0],2)
	# en=0.
	# for x in lat.allspins():
	# 	q=lat.spinenergy(x)[0]
	# 	en+=q
	# 	print q
	# print 'energy: {0}'.format(en)
	return

class Lattice:
	'Represents pyrochlore lattice, for DSIM KMC'

	# in Angstroms 10**-10 meters (MUST NOT BE int)
	L = 10. #somehow 10 is hardcoded in...

	# just normalizing to unit spins, magnetic moment already includd in J and D
	magmom = 1/np.sqrt(3)
	# in Kelvin
	J = -3.72
	D = 1.41
	rnn = 2.5*np.sqrt(2)

	# controls balance between real and k space sums
	alpha = 5e-3

	# locations of fcc bravais atoms in a unit cell
	fdict={}
	fdict[0]=np.array([0,0,0])
	fdict[1]=np.array([L/2,L/2,0])
	fdict[2]=np.array([0,L/2,L/2])
	fdict[3]=np.array([L/2,0,L/2])

	# additional basis of pyrochlore
	bdict={}
	bdict[0]=np.array([0,0,0])
	bdict[1]=np.array([L/4,L/4,0])
	bdict[2]=np.array([0,L/4,L/4])
	bdict[3]=np.array([L/4,0,L/4])

	# allowed k vectors to sum over
	# allowed q (4 vectors)
	qdict={}
	qdict[0]=np.array([0,0,0])
	qdict[1]=2*np.pi*np.array([1,0,0])/L
	qdict[2]=2*np.pi*np.array([0,1,0])/L
	qdict[3]=2*np.pi*np.array([0,0,1])/L





	# orientation of spins in unit cell
	# orientation will be plus or minus, the former is the 
	# 	direction with positive dot product with <1,1,1> and the latter
	# 	is the direction with negative dot product with <1,1,1>
	originspin=np.array([1,1,1.])	

	# simulation box is (a x a x a) array of unit cells
	def __init__(self, a):
		# side length
		self.a = a 

		self.N = 16*(a**3)

		self.ucells = self.gencubepos()

		# total number of independent spins
		self.N = a*a*a*16 

		self.topoint = self.hashcoords()

		self.spindict = self.spinindex()

		# initialize all positive dir; can also initilize randomly
		self.spinstate = np.ones([a**3,4,4],dtype=bool)

		self.statevector = self.getstatevector()


		# put NNguide generation inside a function

		# for each spin, store the list of nearest neighbors
		# random initialization, with 2 in 2 out rule enforced
		# one strategy: randomly select sites, randomly assign
		# state, unless rule will be violated
		# make a list of all points of spins, randomly select
		# without replacement. terminate when empty

		# spingen = self.allspins()
		# spinlist = []
		# for point in spingen:
		# 	spinlist.append(point)


		# does NNguide return spin, or coords? 
		# it gives spin indices (all ints, in spin space)
		# and takes in spin indices as well

		self.spinvisit = np.zeros([a**3,4,4],dtype=bool)


		self.nnguide1 = np.zeros([a**3,4,4,6,3])
		spingen1 = self.allspins()
		for ipoint in spingen1:
			posi=self.position(ipoint)
			nnlist=[]
			nnori =  (Lattice.L/4)*self.ori(ipoint)
			for x in xrange(3):
				nntrans = np.copy(nnori)
				nntrans[x] = 0. 
				nnc1 = posi+nntrans
				# enforce PBCs
				for w in xrange(3):
					if(nnc1[w]>=self.a*Lattice.L):
						nnc1[w] -= self.a*Lattice.L
					if(nnc1[w]<0):
						nnc1[w] += self.a*Lattice.L
				nnlist.append(self.topoint[tuple(nnc1)])
				nnc2 = posi - nntrans
				for w in xrange(3):
					if(nnc2[w]>= self.a*Lattice.L):
						nnc1[w] -= self.a*Lattice.L
					if(nnc2[w]<0):
						nnc2[w] += self.a*Lattice.L
				nnlist.append(self.topoint[tuple(nnc2)])
			for x in xrange(len(nnlist)):
				self.nnguide1[ipoint[0],ipoint[1],ipoint[2],x] = nnlist[x]

		#PROBLEM: nbrs of first tetra are always 0
		# first three indices are coords of spin, next is which tetrahedron,
		# then the list of 3 neighbors' coords
		#check if its neighbor is a neighbor of another neighbor,
		# this is how to separate the two tetra

		# self.nnguide2 = np.zeros([a**3,4,4,2,3,3])
		# spingen2 = self.allspins()
		# for ipoint in spingen2:
		# 	nbpoint1 = self.nnguide1[ipoint[0],ipoint[1],ipoint[2],0]
		# 	self.nnguide2[ipoint[0],ipoint[1],ipoint[2],0,0] = nbpoint1
		# 	m=1
		# 	n=0
		# 	for q in xrange(5):
		# 		nbpoint2 = self.nnguide1[ipoint[0],ipoint[1],ipoint[2],q+1]
		# 		aa = self.nnguide1[nbpoint1[0],nbpoint1[1],nbpoint1[2]]
		# 		bb = self.nnguide1[nbpoint2[0],nbpoint2[1],nbpoint2[2]]
		# 		if(self.arraycont(aa,bb)>=2):
		# 			self.nnguide2[ipoint[0],ipoint[1],ipoint[2],0,m] = nbpoint2
		# 			m+=1
		# 			# print 'tetra0'
		# 		else:
		# 			self.nnguide2[ipoint[0],ipoint[1],ipoint[2],1,n] = nbpoint1
		# 			n+=1


		# allowed Kn (a**3 vectors)
		aa=4*a
		if(aa%2==1):
			self.na=np.arange((1-aa)/2,(1+aa)/2)
		else:
			self.na=np.arange((2-aa)/2,(2+aa)/2)
		self.na=self.na[::-1]
		na2 = self.na.copy()
		naind = np.argsort(np.abs(na2))
		self.naord = np.zeros([len(na2)])
		for x in xrange(len(na2)):
			self.naord[x] = self.na[naind[x]]

		self.naord2 = np.arange(aa)

	def getstatevector(self):
		statevector = np.zeros([self.N],dtype=bool)
		for x in xrange(self.N):
			spinx = self.spindict[x]
			statevector[x] = self.spinstate[spinx[0],spinx[1],spinx[2]]
		return statevector

	# The two functions below return all n-tuples which 
	# sum to a given number, with a specified maximum
	def rec_fun(self,sum,deepness,myList,Total,cap):
	    if deepness==0:
	        if sum==Total:
	            yield myList
	    else:
	        for i in xrange(min(cap,Total-sum+1)):
	            for b in self.rec_fun(sum+i,deepness-1,myList+[i],Total,cap):
	                yield b

	def getindgen(self,digits,Tot,cap):
	    indgen=self.rec_fun(0,digits,[],Tot,cap+1)
	    return indgen

	def allKind(self,cap):
		# krad = self.a-1
		# assert(cap < self.a)
		# for n in xrange((3*(self.a-1))+1):
		# 	g=self.getindgen(3,n,cap)
		# 	for m in g:
		# 		yield m

		# set cap to 2*a to sample 2nd BZ
		for n in xrange((3*(cap+1-1))+1):
			g=self.getindgen(3,n,cap)
			for m in g:
				yield m

	def allKn(self,cap):
		g=self.allKind(cap)
		for ind in g:
			da=np.double(self.a)
			kvec = np.array([self.naord2[ind[0]]/da,self.naord2[ind[1]]/da,self.naord2[ind[2]]/da])*(2*np.pi/Lattice.L)
			yield kvec

	def alln(self,cap):
		g=self.allKind(cap)
		for ind in g:
			da=np.double(self.a)
			nvec = da*Lattice.L*np.array([self.naord[ind[0]],self.naord[ind[1]],self.naord[ind[2]]])
			yield nvec

	def allKvec(self,cap):
		h=self.allKn(cap)
		for vec in h:
			for p in xrange(4):
				yield vec+Lattice.qdict[p]


	def arraycont(self,list1,list2): #lists of arrays, how many arrays do they have in common? (no repeats in either list)
		counter = 0
		for x in list1:
			for y in list2:
				if((x==y).all()):
					counter+=1
		return counter


	# for random initialization with 2 in 2 out
	# Do in style of BFS, start at the same point every time
	def randinit(self):
		w = np.array([0,0,0])
		#use a FIFO queue to iterate through spins
		setqueue = Queue.LifoQueue()
		setqueue.put(w)
		while(np.invert(self.spinvisit).any()):
			q = setqueue.get()
			if(self.spinvisit[q[0],q[1],q[2]]):
				continue

			spinset,center = self.checkspin(q)
			print spinset
			assert(spinset != -1)

			# if z is true, then spinstate[q] being true means
			# spin q points out of the tetra with center at center
			z = (np.dot(self.ori(q),self.position(q.astype(int))-center) > 0)

			if(spinset==1):
				self.spinstate[q[0],q[1],q[2]] = (True==z)
				self.spinvisit[q[0],q[1],q[2]] = True
				print 'true'
			elif(spinset==2):
				self.spinstate[q[0],q[1],q[2]] = (False==z)
				self.spinvisit[q[0],q[1],q[2]] = True
				print 'false'

			# need one more condition in here, if unconstrained, set it randomly
			if(spinset==0):
				self.spinstate[q[0],q[1],q[2]] = (random.random() > 0.5)
				self.spinvisit[q[0],q[1],q[2]] = True

			setqueue.put(self.nnguide2[q[0],q[1],q[2],0,0])
			setqueue.put(self.nnguide2[q[0],q[1],q[2],0,1])
			setqueue.put(self.nnguide2[q[0],q[1],q[2],0,2])

	# func that checks whether a spin is constrained, and if so,
	# sets it to the appropriate value. if somehow  a spin cannot
	# be in either state, contradiction, returns an error
	# not constrained = 0
	# constrained:
	# out wrt tetra 1: 1
	# in wrt tetra 1 : 2
	# out wrt tetra 2: 3 (same as 2)
	# in wrt tetra 2 : 4 (same as 1)
	# contradiction = -1
	def checkspin(self,cpoint):
		# cpoint = self.topoint[tuple(p)]
		out = 0
		constr = False
		tetra1 = self.nnguide2[cpoint[0],cpoint[1],cpoint[2],0]
		tetra2 = self.nnguide2[cpoint[0],cpoint[1],cpoint[2],1]
		inflag = 0
		outflag = 0
		center1 = np.mean(np.vstack([cpoint,tetra1]),0)
		for p in xrange(3):
			sp = tetra1[p]
			if(self.spinvisit[sp[0],sp[1],sp[2]]):
				pointing = self.ori([sp[0],sp[1],sp[2]])
				if(self.spinstate[sp[0],sp[1],sp[2]]==False):
					pointing = -pointing
				if(np.dot(cpoint-center1,pointing)>0):
					outflag += 1
				else:
					inflag += 1
		if((outflag>=2) and (inflag>=2)):
			return -1,center1
		elif(outflag>=2):
			constr = True
			out = 1
		elif(inflag>=2):
			constr = True
			out = 2

		out2 = 0
		inflag = 0
		outflag = 0
		constr2 = False
		center2 = np.mean(np.vstack([cpoint,tetra1]),0)
		for p in xrange(3):
			sp = tetra2[p]
			if(self.spinvisit[sp[0],sp[1],sp[2]]):
				pointing2 = self.ori([sp[0],sp[1],sp[2]])
				if(np.dot(cpoint-center2,pointing2)>0):
					outflag += 1
				else:
					inflag += 1
		if((outflag>=2) and (inflag>=2)):
			return -1,center1
		elif(outflag>=2):
			constr2 = True
			out2 = 2#3
		elif(inflag>=2):
			constr2 = True
			out2 = 1#4

		if((constr and constr2) and (out != out2)):
			return -1,center1
		return out,center1


	def gencubepos(self):
		# going to output ordered list of 3-tuples, which are coordinates of 
		# 	the corner of the unit cell closest to the origin
		a=self.a
		coos=[]
		for k in xrange(1,a+1):
			if(k==1):
				coos.append((0,0,0))
				continue
			for x in xrange(k):
				for y in xrange(k):
					coos.append((x,y,k-1))
			for y in xrange(k):
				for z in xrange(k-1):
					coos.append((k-1,y,z))
			for x in xrange(k-1):
				for z in xrange(k-1):
					coos.append((x,k-1,z))
		return coos

	# generate hash from coords to point
	def hashcoords(self):
		topoint = {}
		spiniter = self.allspins()
		for spin in spiniter:
			pos = self.position(spin)
			topoint[tuple(pos)] = spin
		return topoint

	def spinindex(self):
		spindict = {}
		ind = 0
		for spin in self.allspins():
			spindict[ind] = spin
			ind+=1
		return spindict

	def state(self,point): 
		return self.spinstate[point[0]][point[1]][point[2]]

	# used to return Cartesian coords of a (point described by 3-tuple)
	# MODIFY to take into account PBCs
	def position(self,point):
		pos = Lattice.L*np.array(self.ucells[point[0]])
		pos += Lattice.fdict[point[1]]
		pos += Lattice.bdict[point[2]]
		return pos

	# return moment of spin at location in unit cell
	# (which unit cell doesn't matter)

	# CHECK that the direction convention holds (seriously though)
	# ie always return <1,1,1>
	# normalized s.t. each entry is +-1
	def ori(self,point):
		# reverse sign of entry where basisvec is 0
		dir = np.copy(Lattice.originspin)
		if(point[2]==0):
			pass
		elif(point[2]==1):
			dir[2] *= -1
		elif(point[2]==2):
			dir[0] *= -1
		elif(point[2]==3):
			dir[1] *= -1
		return dir

	# return a generator that iterators over all spins in the Lattice
	# return the point array
	def allspins(self):
		for m in xrange(self.a**3):
			for n in xrange(4):
				for o in xrange(4):
					yield np.array([m,n,o])


	def cellspins(self,c):
		for n in xrange(4):
			for o in xrange(4):
				yield np.array([c,n,o])


	# apply PBCs
	# use while instead of if?
	def applypbc(self,nnc1):
		for w in xrange(3):
			while(nnc1[w]>=self.a*Lattice.L):
				nnc1[w] -= self.a*Lattice.L
			while(nnc1[w]<0):
				nnc1[w] += self.a*Lattice.L
		return nnc1

	# yields translation vectors for #layers shells
	def layerspins(self,layers):
		for x in xrange(-layers,layers+1):
			for y in xrange(-layers,layers+1):
				for z in xrange(-layers,layers+1):
					trans = Lattice.L*np.array([x,y,z])
					if(not trans.any()):
						continue
					yield trans


	def stateori(self,ipoint):
		if(self.state(ipoint)):
			return self.ori(ipoint)
		else:
			return -self.ori(ipoint)


	# calculate energy of spin via brute force, for testing
	# sums over #layers shells of images, for dipenergy
	#for each layer, generate coords of spin in nbring cell,
	#apply PBCs to translate to a positive coord, then
	#use topoint to get the corresponding spin's state
	# 0<layers
	def brute(self,ipoint,layers):
		mui = self.stateori(ipoint)
		posi = self.position(ipoint)

		dipenergy=0.
		for jpoint in self.cellspins(ipoint[0]):
			for trans in self.layerspins(layers):
				posj = self.position(jpoint)
				newposj = posj+trans
				# print posj,trans,newposj
				newposjcopy = np.copy(newposj)
				newjpoint = self.topoint[tuple(self.applypbc(newposjcopy))]
				muj = self.stateori(newjpoint)
				rij = posi - newposj
				# print posj,trans,newposj
				d1 = np.dot(mui,muj)/(norm(rij)**3)
				d2 = 3*np.dot(mui,rij)*np.dot(muj,rij)/(norm(rij)**5)
				dipenergy += Lattice.D*(Lattice.rnn**3)*(d1-d2)

			if((jpoint-ipoint).any()==False):
				continue
			rij = posi - self.position(jpoint)
			muj = self.stateori(jpoint)
			d1 = np.dot(mui,muj)/(norm(rij)**3)
			d2 = 3*np.dot(mui,rij)*np.dot(muj,rij)/(norm(rij)**5)
			dipenergy += Lattice.D*(Lattice.rnn**3)*(d1-d2)

		nnenergy=0.
		for nn in xrange(6):
			nnpoint = self.nnguide1[ipoint[0],ipoint[1],ipoint[2],nn]
			munn = self.stateori(nnpoint)
			nnenergy += -Lattice.J*np.dot(mui,munn)

		return np.array([dipenergy+nnenergy,dipenergy,nnenergy])


	# implement for one spin at a time
	# for each spin, similar to krad or layers,
	# return only those within distance L
	def realenergy(self):
		renergy=0.
		for ipoint in self.allspins():
			mui = self.stateori(ipoint)
			posi = self.position(ipoint)
			for nn in xrange(6):
				nnpoint = self.nnguide1[ipoint[0],ipoint[1],ipoint[2],nn]
				munn = self.stateori(nnpoint)
				posnn = self.position(nnpoint.astype(int))
				rij = posi-posnn
				renergy+=np.dot(mui,munn)*self.B(norm(rij))
				renergy+=-np.dot(mui,rij)*np.dot(munn,rij)*self.C(norm(rij))
		renergy*=0.5
		return renergy
	

	def reciprocalenergy(self,cap):
		kenergy=0.
		for kvec in self.allKvec(cap):
			if(kvec.any()==False):
				continue
			Qk=0.
			for ipoint in self.allspins():
				fact1=np.dot(self.stateori(ipoint),kvec)
				fact2=np.exp(1j*np.dot(self.position(ipoint),kvec))
				Qk+=fact1*fact2
			# print Qk
			# print (norm(kvec)**-2)*np.exp(-(norm(kvec)**2)/(4*Lattice.alpha))
			# print r'\n'
			# print np.abs(Qk)**2
			kenergy+=((norm(kvec)**-2)*np.exp(-(norm(kvec)**2)/(4*Lattice.alpha))*
			(np.abs(Qk)**2))
		kenergy*=(2*np.pi/(Lattice.L**3))
		return kenergy

	#is the self energy term part of this? is it part of Ewald
	def dipenergy(self,cap):
		d1=self.realenergy()
		d2=self.reciprocalenergy(cap)
		d=Lattice.D*(Lattice.rnn**3)*(d1+d2)
		return d


	def NNenergy(self):
		nnenergy = 0.
		for ipoint in self.allspins():
			mui = self.stateori(ipoint)
			for nn in xrange(6):
			    nnpoint = self.nnguide1[ipoint[0],ipoint[1],ipoint[2],nn]
			    munn = self.stateori(nnpoint)
			    nnenergy += -Lattice.J*np.dot(mui,munn)
		return nnenergy

	# independent of Ewald; b/c surrounded by vacuum
	def dotenergy(self):
		denergy=0.
		for ipoint in self.allspins():
			mui = self.stateori(ipoint)
			for jpoint in self.allspins():
				muj = self.stateori(jpoint)
				denergy+= np.dot(mui,muj)
		return denergy*2*np.pi/(3*(Lattice.L**3))


	# if alpha is large, first term (real space sum) is only for n=0.
	# optimal is alpha=5/L (?)
	# sum over 100-200 terms (total?) for k-space sum
	def B(self,r):
		alpha = Lattice.alpha
		num = erfc(alpha*r)
		num += (2*alpha*r/np.sqrt(np.pi))*np.exp(-alpha*alpha*r*r)
		den = r*r*r
		return num/den

	def C(self,r):
		alpha = Lattice.alpha
		num = 3*erfc(alpha*r)
		num += ((2*alpha*r/np.sqrt(np.pi))*(3 + 2*alpha*alpha*r*r)*
			np.exp(-alpha*alpha*r*r))
		den = r**5
		return num/den


	def spinenergy(self,ipoint,cap):
		nlimit = 2

		mui = self.stateori(ipoint)
		posi = self.position(ipoint)

		dotenergy=0.
		spiniter = self.allspins()
		for jpoint in spiniter:
		    if((jpoint-ipoint).any()==False):
		        continue
		    rij = posi - self.position(jpoint)
		    muj = self.stateori(jpoint)

		    # real space sum first, do only for n=0 
		    realenergy = 0.
		    for nx in xrange(nlimit):
		        for ny in xrange(nlimit):
		            for nz in xrange(nlimit):
		                nvec = self.a*Lattice.L*np.array([nx,ny,nz])                        
		                realenergy += (np.dot(mui,muj) *
		                    self.B(norm(rij + nvec)))
		                realenergy += (-np.dot(mui,rij+nvec) * 
		                    np.dot(muj,rij+nvec) * self.C(norm(rij + nvec)))
		    realenergy *= 0.5
		    # print 'realdone'

		    # reciprocal space sum
		    # need to sum over allowed k vectors
		    # most expensive part right here. Cython?
		    kenergy = 0.
		    for kvec in self.allKvec(cap):
		        if((kvec[0]==0) and (kvec[1]==0) and (kvec[2]==0)):
		            continue
		        kenergy += (np.dot(mui,kvec)*np.dot(muj,kvec)*
		        np.exp(-norm(kvec)*norm(kvec)/(4*Lattice.alpha))*
		        np.cos(np.dot(kvec,rij)))
		    kenergy *= (2*np.pi/(Lattice.L**3))
		    # print 'kdone'

		    dotenergy += 0# np.dot(mui,muj)*2*np.pi/(3*(Lattice.L**3))
		    # print 'dotdone'

		# NN term
		nnenergy = 0.
		for nn in xrange(6):
		    nnpoint = self.nnguide1[ipoint[0],ipoint[1],ipoint[2],nn]
		    munn = self.stateori(nnpoint)
		    nnenergy += -Lattice.J*np.dot(mui,munn)
		# print 'nndone'

		#need to still add boundary (eps=1 or Inf) term (?)
		dipenergy = Lattice.D*(Lattice.rnn**3)*(realenergy+kenergy+dotenergy)
		totalenergy = dipenergy + nnenergy
		# totalenergy=nnenergy
		return totalenergy,realenergy,kenergy,dotenergy,nnenergy

	# test this mofo out
	def isNbr(self,ipoint,jpoint): 
		n=list(self.nnguide1[jpoint[0],jpoint[1],jpoint[2]])
		return np.array([(ipoint==x).all() for x in n]).any()

	# multiprocessing? Cython?
	# calculate energy with both true, then with one false.
	def buildintmat(self,nrad,krad):
		# optimize the size of this guy later
		intmat = np.zeros([self.N,self.N],dtype=complex) 

		for x in xrange(self.N):
			xpoint = self.spindict[x]
			mux = self.ori(xpoint)
			for y in xrange(self.N):#xrange(x+1,self.N):
				if(x==y):
					continue
				ypoint = self.spindict[y]
				rij = self.position(xpoint) - self.position(ypoint)
				muy = self.ori(ypoint)

				#dipolar energy
				# real energy
				ren = 0.
				for n in self.alln(nrad):
					n = self.a*Lattice.L*np.array(n)
					term1 = np.dot(mux,muy)*self.B(norm(rij+n))
					term2 = -(np.dot(mux,rij+n)*np.dot(muy,rij+n))*self.C(norm(rij+n))
					if(abs(term1+term2)>1000):
						print y,term1,term2,n
					ren += term1+term2
				ren *= 0.5

				# k energy
				ken = 0.
				for kvec in self.allKvec(krad):
					if(kvec.any()==False):
						continue
					k2 = norm(kvec)**2
					# print k2
					fact1 = np.dot(mux,kvec)*np.dot(muy,kvec)
					fact2 = np.exp(-k2/(4*Lattice.alpha))/k2
					# print fact2*k2
					fact3 = np.exp(1j*np.dot(kvec,rij)) #see what this comes out to
					ken += fact1*fact2*fact3
				ken *= 2*np.pi/(Lattice.L**3)

				den = np.dot(mux,muy)*2*np.pi/(3*(Lattice.L**3))

				dipen = Lattice.D*(Lattice.rnn**3)*(ren+ken) #not adding den here; tin-foil condition

				nnen = Lattice.J*np.dot(mux,muy)*self.isNbr(xpoint,ypoint)

				TTen = dipen + nnen

				intmat[x,y] = np.real(TTen)
		return np.real(intmat)

	# write the function that actually gives the energy
	def getenergy(self,J,svector):
		s=svector*2
		s=s-1
		fact1 = np.dot(J,s)
		en = np.dot(s,fact1)
		return en


	def flip_single(self,spinindex,J,T):
		s=self.getstatevector()

		E0 = self.getenergy(J,s)

		s2 = s.copy()
		s2[spinindex] = ~s2[spinindex]
		E1 = self.getenergy(J,s2)

		spin = self.spindict[spinindex]

		flipped = False

		if(E1 < E0):
			self.spinstate[spin[0],spin[1],spin[2]] = ~self.spinstate[spin[0],spin[1],spin[2]]
			flipped = True
		else:
			r = np.random.random()
			prob = np.exp(-(E1-E0)/T)
			if(r < prob):
				self.spinstate[spin[0],spin[1],spin[2]] = ~self.spinstate[spin[0],spin[1],spin[2]]
				flipped = True
		return flipped

	def MCstep(self,J,T):
		counter = 0.
		for x in xrange(self.N):
			ind = np.random.randint(self.N)
			flip = self.flip_single(ind,J,T)
			counter += flip
		frac = counter/self.N
		return frac


	def autocorr(self,initc):
		ac = 0.
		for ind in xrange(self.N):
			st = self.stateori(self.spindict[ind])			
			s0 = initc[ind]
			ac += np.dot(st,s0)/3
		return ac/self.N


	def getautocorr(self,J,T,trange):
		# trange = 50
		ac = np.ones([trange])
		initc = [self.stateori([spin[0],spin[1],spin[2]]) for spin in self.allspins()]
		ac[0] = self.autocorr(initc)
		for steps in xrange(1,trange):
			self.MCstep(J,T)
			ac[steps] = self.autocorr(initc)
		return ac



	def equilibrate(self,J,T):
		initc = [self.stateori([spin[0],spin[1],spin[2]]) for spin in self.allspins()]
		while(self.autocorr(initc) > 0.1):
			print self.autocorr(initc)
			self.MCstep(J,T)
		return



if __name__ == '__main__':
	main()
