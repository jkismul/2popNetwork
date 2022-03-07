import scipy.io
import io
import pandas
import pickle
import numpy
from pylab import *
import matplotlib.pyplot as plt

  # picklelist = [spikes,spikedCells,vSoma,[x.hname() for x in h.nclist],array(h.excvec)]



filename = 'spikes_raster_test'
with open(filename+".sav", 'rb') as fb:
  unpickledlist = pickle.load(fb, encoding='latin1')
# print(unpickledlist)
Nplaced = 0
spikedCells_all = []
for j in range(0,len(unpickledlist[1])): #j goes over time steps (258)
  spikedCells = unpickledlist[1][j]
  spikedCellsUnique = unique(spikedCells)
  spikedCells2 = zeros(spikedCells.shape)
  for i in range(0,len(spikedCellsUnique)):
    spikedCells2[spikedCells == spikedCellsUnique[i]] = Nplaced + i
  Nplaced = Nplaced + len(spikedCellsUnique)
  spikedCells_all = hstack([spikedCells_all, spikedCells2])
spikes = [hstack(unpickledlist[0]),spikedCells_all]

# plt.plot(unpickledlist[0])
# plt.plot(unpickledlist[1])
print(shape(unpickledlist[0]))
print(shape(unpickledlist[1]))
print(shape(unpickledlist[2]))
print(unpickledlist[0])
print(unpickledlist[1])
print(spikes)
# print(spikes)
# print(unpickledlist[1])
# print(shape(spikes))
# print(shape(unpickledlist[0]))
print(shape(unpickledlist[1][0]))

plt.plot(unpickledlist[0],unpickledlist[1],'b.') 

plt.show()
# for i in range(spikes):
# 	print(i)
# 	plt.plot(0.01*spikes+i)