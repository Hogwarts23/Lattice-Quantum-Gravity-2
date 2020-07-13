import numpy as np
from shell import shelling
import math
from scipy.sparse import csc_matrix, linalg as sla
from scipy.sparse.linalg import inv
# import scipy.linalg as lin
# import numpy.linalg as lin1
###########################################
import matplotlib.pyplot as plt
import math
import time
import random
from pathlib import Path

class LQGConfigs():
	def __init__(self,tfname,masslist,nsources,sourcedir=None,usesource=False):
		self.masslist = masslist #the bare mass square
		self.fname = tfname
		self.nlat = None
		self.label = None
		self.contentN4inf = None
		self.matrix = None
		# self.tmp = np.diag(np.ones(length))
		self.lu = {}
		self.prop = None
		self.maxlen = 50
		self.nsources = nsources
		self.sourcedir = sourcedir
		self.sources = np.zeros(nsources)
		self.usesource = usesource

	def readpartsfromfile(self, keyword):
		file = open(self.fname,'r')
		content = file.read()
		post = content.find(keyword)
		past = content.find("\n\n",post+len(keyword))
		content = content[post+len(keyword):past]
		return content


	def findcol(self, row, collabel, length):
		l = row
		r = min(collabel,length-1)
		#Trying to bias the binary search
		ratio = (length-1)/self.label[length-1]
		while l <= r:
			mid = l + math.floor((r - l)/2)
			if self.label[mid] == collabel:
				return mid
			elif self.label[mid] < collabel:
				l = mid + 1
			else:
				r = mid - 1


	def cstructpropa(self):
		content = self.readpartsfromfile('N4inf: ').split()
		ind1 = 1
		ind2 = 0
		length = int(content[0])
		self.nlat = length
		self.label = np.zeros(length)
		matrix = np.diag((5)*np.ones(length))
		while ind2 < length:
			self.label[ind2] = int(content[ind1])
			ind1 = ind1 + 6
			ind2 = ind2 + 1
		row = 0
		ind = 1
		while row < length:
			for i in range(1,6):# range(1,6) gives you 1 2 3 4 5
				collabel = int(content[ind+i])
				if self.label[row] < collabel:
					col = self.findcol(row, collabel,length)
					matrix[row,col] = matrix[row,col] - 1
					matrix[col,row] = matrix[col,row] - 1
			row = row + 1
			ind = ind + 6
		self.contentN4inf = content
		self.matrix = matrix

	def computelu(self):
		for m in self.masslist:
			x = self.matrix
			for i in range(self.nlat):
				x[i,i] = x[i,i]+m**2
			smat = csc_matrix(x)
			lu = sla.splu(smat)
			self.lu[m] = lu
			for i in range(self.nlat):
				x[i,i] = x[i,i]-m**2

	def propagator(self,source,mass):
		# l = len(sources)
		lu = self.lu[mass]
		tmp = np.zeros(self.nlat)
		# self.props = np.zeros((self.nlat,l))
		# i=0
		# for m in sources:
		tmp[source] = 1
		x = lu.solve(tmp)
		#print(x)
		self.prop = x

	def cstructshelling(self,origin): #origin is the row number
		if self.label is None:
			self.cstructpropa()
		length = self.nlat
		#set up the first element
		tmpdistance = np.zeros(length)+self.label[length-1]+1
		shell = [] 
		tmp = shelling(origin)
		tmp.addlat(origin)
		shell.append(tmp)
		tmpdistance[origin] = 0
		#construct the matrix
		nleft = length-1 # - 1 because we have already add the origin to the shelling
		shellnum = 0
		while nleft != 0:
			shellnum = shellnum + 1
			tmp = shelling(origin,shellnum)
			tmplat = shell[shellnum-1].shelllattice
			for i in tmplat:
				for j in range(1,6):
					col = self.findcol(0,int(self.contentN4inf[int(i*6+1+j)]),length)
					if tmpdistance[col] > shellnum:
						tmp.addlat(col)
						tmpdistance[col] = shellnum
						nleft = nleft-1
			shell.append(tmp)
		return shell

	def correlators(self,two=True):
		if self.label is None:
			self.cstructpropa()
		l = len(self.masslist)
		totalcor = np.zeros((self.maxlen,self.nsources*l))
		if two:
			totalcor2 = np.zeros((self.maxlen,self.nsources*l))
		self.computelu()
		for j in range(l):
			#t2 = time.time()
			for i in range(self.nsources):
				self.generatesource(self.masslist[j])
				s = self.cstructshelling(int(self.sources[i]))
				ori = s[0].origin
				self.propagator(ori,self.masslist[j])
				cor = self.correlator(s)
				totalcor[:,int(i+self.nsources*j)] = cor
				if two:
					cor2 = self.twoparticlecorrelator(s)
					totalcor2[:,int(i+self.nsources*j)] = cor2
		if two:
			return totalcor,totalcor2
		else:
			return totalcor


	def twoparticlecorrelator(self,slist):
		pro2 = self.prop**2
		#print(pro)
		#np.savetxt('test1',pro)
		num = 0
		cor2 = np.zeros(self.maxlen)
		for shell in slist:
			if num>=self.maxlen:
				break
			s = 0
			for index in shell.shelllattice:
				s = s + pro2[int(index)]
			s = s/shell.numlat
			cor2[num] = s
			num = num + 1
		return cor2

	def correlator(self,slist):
		pro = self.prop
		#print(pro)
		#np.savetxt('test1',pro)
		num = 0
		cor = np.zeros(self.maxlen)
		for shell in slist:
			if num>=self.maxlen:
				break
			s = 0
			for index in shell.shelllattice:
				s = s + pro[int(index)]
			s = s/shell.numlat
			cor[num] = s
			num = num + 1
		return cor

	def findmax(self,lst):
		maximum = lst[0]
		index = 0
		ind = -1
		for i in lst:
			ind = ind + 1
			if i > maximum:
				maximum = i
				index = ind
		return index, maximum

	def generatesource(self,mass): # return the distance
		if self.usesource:
			self.sources = np.load(Path(self.sourcedir)/(self.fname.split('/')[-1]+'%6f.npy'%mass))
		else:
			s = self.cstructshelling(random.randint(0,self.nlat-1))
			lent = len(s)
			y = np.zeros(lent)
			for i in range(lent):
				y[i] = s[i].numlat
			dis,maximum = self.findmax(y)
			lst1 = s[dis].shelllattice
			l = len(lst1)
			for i in range(self.nsources):
				self.sources[i] = lst1[random.randint(0,l-1)]
			np.save(Path(self.sourcedir)/(self.fname.split('/')[-1]+'%6f'%mass),self.sources)


		
