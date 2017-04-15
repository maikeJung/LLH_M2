# Plot arrival time spectrum
# Set mass, dist, number of events and detector type

# imported from c-program
import spectrum

import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from numpy import transpose
from matplotlib.colors import LogNorm
from scipy import interpolate

# set grid size - must match those defined in spectrum.h
REST = 100
RESE = 600
EMAX = 60.0
TMAX = 10.0
STEPE = EMAX/RESE
STEPT = TMAX/REST

# read in event file to plot ontop of the probability distribution
#drawnEvents = np.loadtxt('DATA_CHECK_EVENTS/events1.0eV_1.0Mpc_100.txt')
#energies = []
#times = []
#for j in range(160):
#    energies.append(drawnEvents[j][1])
#    times.append(drawnEvents[j][0])

# set plot parameters
mass = 2.5
distance = 0.5
events = 650
noise_events = 0.0
energyRes = True
triggerEff = True
#noise = 0.002
noise = pow(10,-3)*STEPT
#noise = 0.0

# create spectrum
spectrum_plot = spectrum.doubleArray( (RESE - 1) * REST )
# INPUT: array, mass, distance, #of events, bool energyRes, bool triggerEff, noise

pylogTime = np.logspace(-5.0,1.0,num=REST)
logTime = spectrum.doubleArray(REST)
for i in range(len(pylogTime)):
    logTime[i] = pylogTime[i]
spectrum.createSpectrum(spectrum_plot ,mass ,distance ,events ,energyRes , triggerEff ,noise, noise_events, logTime)

timeArray = range(0, REST, 1)
energyArray = range(0, RESE-1)
myArray = [[spectrum_plot[t*(RESE-1) +e] for t in timeArray] for e in energyArray]
#myArray = [[spectrum_plot[t*(RESE-1) +e] for t in range(0, REST, 1)] for e in range(0, RESE-1)]
#X = np.arange(0, TMAX, STEPT)
#Y = np.arange(STEPE, EMAX-STEPE, STEPE)
#Y = np.arange(0, EMAX-STEPE, STEPE)

tint = [timeArray[i]*STEPT for i in timeArray]
eint = [energyArray[i]*STEPE + STEPE for i in energyArray]

#X, Y = np.meshgrid(X, Y)
X, Y = np.meshgrid(tint, eint)

Z = myArray



'''
#Interpolate
x1 = np.asarray(np.arange(0, TMAX, STEPT))
y1 = np.asarray(np.arange(STEPE, EMAX-STEPE, STEPE))
z1 = transpose(myArray)
sp = interpolate.RectBivariateSpline(x1, y1, z1, kx=3, ky=3, s=0)

print float(sp([0.6012], [0.25]))  # array([[ 7.3]])
print sp([8.0],[12.1])  # array([[ 2.66427408]])

#myArrayTest = [ [1]  for t in np.arrange(0, TMAX, REST)] for e in np.arrange(0, EMAX-RESE,RESE)]
myArrayTEST = [[float(sp([t],[e])) for t in np.arange(0, TMAX, STEPT)] for e in np.arange(STEPE, EMAX-STEPE, STEPE)]
'''
'''
diff = 0.0
for t in range(0,REST):
    for e in range(0, RESE-2):
        diff += abs(Z[e][t] -  myArrayTEST[e][t])
print diff
'''

# Surface Plot of the arrival distribution - log scale
fig = plt.figure()
ax = fig.add_subplot(111)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(top=0.87)
#surf = ax.contourf(X,Y,Z, 8, cmap=plt.cm.jet)
surf = ax.contourf(X,Y,Z,levels=[1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1],cmap=plt.cm.jet,norm = LogNorm())
#plt.plot(times,energies, 'ro')
#surf = ax.contourf(X,Y,Z,cmap=plt.cm.jet,norm = LogNorm())
ax.set_xlabel('detection time [s]', fontsize=19)
ax.set_ylabel('energy [MeV]', fontsize=19)
#ax.set_title('m = '+str(M)+' eV - ' + str(events) + ' events - D = '+str(D)+ ' Mpc \n'+str(det_type), fontsize=19)
ax.xaxis.set_tick_params(labelsize=19, width=2)
ax.yaxis.set_tick_params(labelsize=19, width=2)
ax.xaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
# defining custom minor tick locations:
ax.xaxis.set_minor_locator(plt.FixedLocator([50,500,2000]))
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.tick_params(axis='both',reset=False,which='both',length=8,width=2)
cbar = fig.colorbar(surf, shrink=1, aspect=20, fraction=.12,pad=.02)
cbar.set_label('probability',size=19)
# access to cbar tick labels:
cbar.ax.tick_params(labelsize=19)
#plt.xlim(0.0, 10.0)
plt.ylim(1, 39)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(top=0.87)
surf = ax.contourf(X,Y,Z, 8, cmap=plt.cm.jet)
#plt.plot(times,energies, 'ro')
#surf = ax.contourf(X,Y,Z,levels=[1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1],cmap=plt.cm.jet,norm = LogNorm())
#surf = ax.contourf(X,Y,Z,cmap=plt.cm.jet,norm = LogNorm())
ax.set_xlabel('detection time [s]', fontsize=19)
ax.set_ylabel('energy [MeV]', fontsize=19)
#ax.set_title('m = '+str(M)+' eV - ' + str(events) + ' events - D = '+str(D)+ ' Mpc \n'+str(det_type), fontsize=19)
#plt.plot(times,energies, 'ro')
ax.xaxis.set_tick_params(labelsize=19, width=2)
ax.yaxis.set_tick_params(labelsize=19, width=2)
ax.xaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
# defining custom minor tick locations:
ax.xaxis.set_minor_locator(plt.FixedLocator([50,500,2000]))
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.tick_params(axis='both',reset=False,which='both',length=8,width=2)
cbar = fig.colorbar(surf, shrink=1, aspect=20, fraction=.12,pad=.02)
cbar.set_label('probability',size=19)
# access to cbar tick labels:
cbar.ax.tick_params(labelsize=19)
#plt.xlim(0.0, 10.0)
plt.ylim(1, 39)


# Surface Plot of the arrival distribution - log scale
'''
fig = plt.figure()
ax = fig.add_subplot(111)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(top=0.87)
#surf = ax.contourf(X,Y,Z, 8, cmap=plt.cm.jet)
surf = ax.contourf(X,Y,Z,levels=[1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1],cmap=plt.cm.jet,norm = LogNorm())
#plt.plot(times,energies, 'ro')
#surf = ax.contourf(X,Y,Z,cmap=plt.cm.jet,norm = LogNorm())
ax.set_xlabel('time [s]', fontsize=19)
ax.set_ylabel('energy [MeV]', fontsize=19)
#ax.set_title('m = '+str(M)+' eV - ' + str(events) + ' events - D = '+str(D)+ ' Mpc \n'+str(det_type), fontsize=19)
ax.xaxis.set_tick_params(labelsize=19, width=2)
ax.yaxis.set_tick_params(labelsize=19, width=2)
ax.xaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
# defining custom minor tick locations:
ax.xaxis.set_minor_locator(plt.FixedLocator([50,500,2000]))
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.tick_params(axis='both',reset=False,which='both',length=8,width=2)
cbar = fig.colorbar(surf, shrink=1, aspect=20, fraction=.12,pad=.02)
cbar.set_label('# of events',size=19)
# access to cbar tick labels:
cbar.ax.tick_params(labelsize=19)
#plt.xlim(0.0, 10.0)
#plt.ylim(1, 39)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(top=0.87)
surf = ax.contourf(X,Y,Z, 8, cmap=plt.cm.jet)
#surf = ax.contourf(X,Y,Z,levels=[1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1],cmap=plt.cm.jet,norm = LogNorm())
#surf = ax.contourf(X,Y,Z,cmap=plt.cm.jet,norm = LogNorm())
ax.set_xlabel('time [s]', fontsize=19)
ax.set_ylabel('energy [MeV]', fontsize=19)
#ax.set_title('m = '+str(M)+' eV - ' + str(events) + ' events - D = '+str(D)+ ' Mpc \n'+str(det_type), fontsize=19)
#plt.plot(times,energies, 'ro')
ax.xaxis.set_tick_params(labelsize=19, width=2)
ax.yaxis.set_tick_params(labelsize=19, width=2)
ax.xaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
# defining custom minor tick locations:
ax.xaxis.set_minor_locator(plt.FixedLocator([50,500,2000]))
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.tick_params(axis='both',reset=False,which='both',length=8,width=2)
cbar = fig.colorbar(surf, shrink=1, aspect=20, fraction=.12,pad=.02)
cbar.set_label('# of events',size=19)
# access to cbar tick labels:
cbar.ax.tick_params(labelsize=19)
#plt.xlim(0.0, 10.0)
#plt.ylim(1, 39)
'''

'''
# plot original spectrum

def energy(E):
    return (E**3.8)*np.exp(-(1.0+3.8)*E/15.4)/4802.516160

def time(t, energy, mass, dist):
	sigma = 0.9389
	my = -1.0104
	a = 146.2284
	t_delta = dist*51.4635*(mass/energy)**2
	time = t - t_delta
	if time < 0.: return float(0.)
	value = float( (a/146.1884) * np.exp( - (np.log(time) - my)**2/(2*sigma*sigma) ) / (time*sigma*np.sqrt(2*np.pi)))
	return value

def spectrumLL(t,E,mass,dist):
    return energy(E)*time(t,E,mass,dist)


timeArray = np.arange(0.0001, TMAX, STEPT)
energyArray = np.arange(STEPE, EMAX-STEPE, STEPE)
myArray = [[spectrumLL(t,E,mass,distance) for t in timeArray] for E in energyArray]

X, Y = np.meshgrid(timeArray, energyArray)
Z = myArray

fig = plt.figure()
ax = fig.add_subplot(111)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(top=0.87)
surf = ax.contourf(X,Y,Z, 8, cmap=plt.cm.jet)
#plt.plot(times,energies, 'ro')
#surf = ax.contourf(X,Y,Z,levels=[1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1],cmap=plt.cm.jet,norm = LogNorm())
#surf = ax.contourf(X,Y,Z,cmap=plt.cm.jet,norm = LogNorm())
ax.set_xlabel('time after SN [s]', fontsize=19)
ax.set_ylabel('energy [MeV]', fontsize=19)
#ax.set_title('m = '+str(M)+' eV - ' + str(events) + ' events - D = '+str(D)+ ' Mpc \n'+str(det_type), fontsize=19)
#plt.plot(times,energies, 'ro')
ax.xaxis.set_tick_params(labelsize=19, width=2)
ax.yaxis.set_tick_params(labelsize=19, width=2)
ax.xaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
# defining custom minor tick locations:
ax.xaxis.set_minor_locator(plt.FixedLocator([50,500,2000]))
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.tick_params(axis='both',reset=False,which='both',length=8,width=2)
cbar = fig.colorbar(surf, shrink=1, aspect=20, fraction=.12,pad=.02)
cbar.set_label('probability',size=19)
# access to cbar tick labels:
cbar.ax.tick_params(labelsize=19)
#plt.xlim(0.0, 10.0)
plt.ylim(1, 39)
'''

plt.show()

