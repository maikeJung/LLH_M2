import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import scipy
from scipy.interpolate import interp1d

def f(x):
    if x < 600 and x > 500:
        return 1
    else:
        return 0

def g(x):
    if x < 10 and x > 0:
        return 1
    else:
        return 0

a = np.arange(0,1000,1)
b = []
c = []
for i in a:
    b.append(f(i))
    c.append(g(i))
print len(a),len(b)
#print a,b

conv = []
for n in a:
    tosum = 0
    for m in a:
        if n-m >= 0 and n-m < len(a):
            tosum += b[m]*c[n-m]
    conv.append(tosum)

#print conv

test = np.logspace(-5.0,1.0,num=1000)

aa = np.logspace(0.0,3.0, num= 1000)
print test
print len(test)
bb = []
cc = []
for i in aa:
    bb.append(f(i))
    cc.append(g(i))

#test log conv
conv_log = []
for t in range(len(aa)):
    tosum = 0
    for ac in range(len(aa)):
        if ac == 0:
            tosum += 0
        elif t%ac == 0 and ac > 0:
            tosum += bb[int(t/ac)]*cc[ac]*( (aa[ac]-aa[ac-1]) /aa[ac])
    conv_log.append(tosum)

print len(aa), len(bb)

# PLOT
fig=plt.figure()
ax = fig.add_subplot(111)
plt.gcf().subplots_adjust(top=0.85) # to be able to read the lables
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
axes = plt.gca()
axes.xaxis.set_tick_params(labelsize=17, width=1)
axes.yaxis.set_tick_params(labelsize=17, width=1)

#plt.plot(a,b, color="blue", label='original', alpha=0.6)
#plt.plot(a,c, color="green", label='g')
plt.loglog(a,conv, 'r-', label='conv')

#plt.plot(aa,bb,'b--', label='test log')
plt.plot(aa,conv_log, 'm-', label='conv test log')

#plt.xlim(0.6,1.4)
plt.xlabel('mass [eV]', fontsize=18)
plt.ylabel('# of pseudo experiments', fontsize=18)
plt.legend(loc='upper right', fontsize=14)

plt.show()
