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
for j in range(0,len(unpickledlist[1])):
  spikedCells = unpickledlist[1][j]
  spikedCellsUnique = unique(spikedCells)
  spikedCells2 = zeros(spikedCells.shape)
  for i in range(0,len(spikedCellsUnique)):
    spikedCells2[spikedCells == spikedCellsUnique[i]] = Nplaced + i
  Nplaced = Nplaced + len(spikedCellsUnique)
  spikedCells_all = hstack([spikedCells_all, spikedCells2])
spikes = [hstack(unpickledlist[0]),spikedCells_all]

for i in range(len(unpickledlist[2])):
	plt.plot(0.01*unpickledlist[2][i]+i)
  # print(unpickledlist[0][i])
# plt.plot(unpickledlist[2][0])
# plt.plot(unpickledlist[0] ,unpickledlist[2])
# print(unpickledlist[3])
# print(unpickledlist)

#plt.plot(spikes[0][::5],[x+1*150 for x in spikes[1]][::5],'b.',color='#000000',markersize=0.5)
plt.show()


# stuff = unpickledlist[-1][0]
# print(stuff, type(stuff))