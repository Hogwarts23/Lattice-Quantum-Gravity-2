import numpy as np
class shelling():
	def __init__(self, origin, distance=0):
		self.origin = origin
		self.distance = distance
		self.shelllattice = None
		self.numlat = 0

	def addlat(self,num):
		if self.shelllattice is None:
			self.shelllattice = np.zeros(1)+num
		else:
			self.shelllattice = np.append(self.shelllattice,num)
		self.numlat = self.numlat + 1