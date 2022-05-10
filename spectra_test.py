import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import sys
import os

filename = 'RUNS/34295cfefedc46dba14962273d826def/Si.dm.out'
broaden=float(10)
lines=open(filename,'r').readlines()
rawdata = []
for i,line in enumerate(lines):
    if '#' in line:
        for data in lines[i+1:-1]:
            if not data.isspace():
                rawdata.append([data.split()[1],data.split()[4]])
            else:
                break
print(rawdata)
mymin=100000.0
mymax=-100000.0
iy=0.0
peakcent=[]
peakintens=[]
for lines in range(0,len(rawdata)):
   freq=float(rawdata[lines][0])
   if freq < mymin:
      mymin=freq
   if freq > mymax:
      mymax=freq
   peakcent.append(freq)
   peakintens.append(float(rawdata[lines][1].strip('\n')))
if mymin < 0.0:
   mymin=0.0
points=mymax-mymin+100.0
print(peakcent)
print(peakintens)
#determines the x axis of the spectrum, and number of grid points
ix = np.linspace(mymin-50.0,mymax+50.0,points)
                            #peak center positions
for peaks in range(0,len(peakcent)):
   iy+=2.51225*broaden*peakintens[peaks]*mlab.normpdf(ix,peakcent[peaks],broaden)
# newout1=sys.argv[2].replace('raw','fin1')
# newout2=sys.argv[2].replace('raw','fin2')
# np.savetxt(newout2, iy, delimiter=" ")
# np.savetxt(newout1, ix, delimiter=" ")
# var1=open(newout1,'r').readlines()
# var2=open(newout2,'r').readlines()
# finout=sys.argv[2].replace('raw','fin')
# finwrite=open(finout,'w')
# for lines in range(0,len(var1)):
#    finwrite.write('%s %s\n' %(var1[lines].strip('\n'),var2[lines].strip('\n')))
# os.system('rm %s\n' %(newout1))
# os.system('rm %s\n' %(newout2))

fig,ax = plt.subplots()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.xlabel('cm-1')
plt.ylabel('Intensity')
ax.plot(ix,iy)
plt.savefig('spectra.png')
