"""Graphene random trajectories.

This python script is designed to determine the number of states arriving from a source, 
to a detector via two pinhole filters. These filters spatially restrict the trajectories,
so only those satisfying a certain condition can pass through. The method used to 
calculate the trajectories will start from the graphene band-structure, and determine the 
Fermi contour at a given energy.  

Following this, the equation of motion is solved, with a random additional offset added. 
This offset is the guiding center co-ordinate, and is chosen in a region. The resulting data
set is then checked to see whether it passes through to the collector. If it does, 
a count is then recorded. After a given numer of runs, the field is increased, and the 
process repeats. 

As a result, there are two parts. In part 1, the Fermi contour is determined, and a set 
created, (x, y) of all point on the Fermi arc. In part 2, there are two processes. Firstly, 
a magnetic field is chosen, which generates the scale (or size) of the Fermi contour
in real space. Then the guiding centre co-ordinates are chosen at random. Finally, it is 
checked as to whether the given trajectory can pass through all the pinhole filters. 
"""

#######

#Code for Graphene paper on Valley Splitting. Does the random focusing calculations

from scipy.integrate import dblquad
import scipy.special as sc
import scipy.optimize as so
import numpy as np
import math
import cmath
from numpy import loadtxt
from pylab import figure, plot, xlabel, ylabel, xticks, yticks, grid, legend, title, savefig, ylim, xlim, subplots_adjust
from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib.font_manager



#Introduce the dictionary

d = {} 

jmax = 5
#Call the filenames I will use for the two valley states. 
filename1 = "kxkyvalley1"
filename2 = "kxkyvalley2"

#And the filenames for the figures.
#and the series number, q
q=5
figfile1 = "Valleysplitting"+str(q)
figfile2 = "Trajectories"+str(q)

#Physical Parameters for the system
h = 1.054*10**(-34)
e = 1.609*10**(-19)

#Graphene hopping parameters.
t = 2.8*1.6*10**(-19)
tprime = 0.02*t
a = 1.42*10**(-10)

#Graphene Fermi velocity
vF = 3*t*a/(2*h)

#density 
n = 2*10**15

#fermi momentum 
k = np.sqrt(np.pi*n)
print(k)

#Define the location of K and K' gamma points
Kx = 2*np.pi/(3*a)
Ky = 2*np.pi/(np.sqrt(3)*3*a)

#determining the spacing of the points that I calculate the dispersion at
thetaspacing = np.pi/100.

#predefine the densities
narray = np.array([0.2, 0.6, 1.5, 3, 6])*10**(17)
narraylegend = [r'$n = 2 \times 10^{12}$cm$^{-2}$', r'$n = 6 \times 10^{12}$cm$^{-2}$', r'$n = 1.5 \times 10^{13}$cm$^{-2}$', r'$n = 3 \times 10^{13}$cm$^{-2}$', r'$n = 6 \times 10^{13}$cm$^{-2}$']

e0array = vF*h*np.sqrt(np.pi*narray)

def Dispersion(k, args):
	"""Dispersion function

	This function is used in combination with Bisection in the first loop to determine the 
	Fermi contour. 

	It does not include the contribution from t', as this is very small. 
	"""
	a, t, vF, s, theta, e0 = args 
	
	kx = k*np.cos(theta)
	ky = k*np.sin(theta)
	fkk = 2*np.cos(np.sqrt(3)*(-s*ky*a) + 2*np.pi/3) + 4*np.cos(np.sqrt(3)*(-s*ky*a/2)+ np.pi/3)*np.cos(3*(kx)*a/2 + np.pi)
	Dis = t*np.sqrt(3+fkk) - e0
	
	return Dis

def Arraydiscarder(arrayx, arrayy, args):
	"""Array analysis

	This function takes in two n x 1 arrays, which defines the discretised trajectory
	of the electron. If the electron doesn't pass through a specific region, the 
	function returns 0. If it does, it returns 1. 

	The function takes arguments, L, the length of the collimator, l, the focusing length, 
	and w, the width of the aperture. Finally there is delta, which is necessary due to the 
	discrete nature of the problem. 
	"""
	L, l, w, delta = args
	j=0
	for j in range(len(arrayx)):
		x = arrayx[j]
		y = arrayy[j]
		if y > -L/2 - delta and y < -L/2 + delta:
			if x > l - w/2 and x < l + w/2: 		
				u = 1
			else:
				u = 0
				break
		if y > L/2 - delta and y < L/2 + delta:
			if x > l - w/2 and x < l + w/2: 
				u = 1
			else:
				u = 0
				break	
		if x > 0 - delta and x < 0 + delta:
			if y > l - w/2 and y < l + w/2: 
				u = 1
			else:
				u = 0
				break
#		if x > L/2 - delta and x < L/2 + delta:
#			if y > l - w/2 and y < l + w/2: 
#				u = 1
#			else:
#				u = 0
#
#				break
#		if x > -L/2 - delta and x < -L/2 + delta:
#			if y > l - 2*w/3 and y < l + w/3: 
#				u = 1
#			else:
#				u = 0
#				break
		j= j+1
	return u, j	

def running_mean(x, N):
	"""Running mean

	This smooths out the stochastic noise of the random focusing output. Input is a numpy array, 
	x, and a scalar, N. 
	"""
	cumsum = np.cumsum(np.insert(x, 0, 0))
	return (cumsum[N:] - cumsum[:-N]) / float(N)
		


j=1
while j <= jmax:
	"""
	First loop Fermi energy
	
	This is the main loop for determining the Fermi contour in a particular valley, for a given 
	Fermi energy. 

	As the Fermi energy increases, the Fermi contour changes.
	"""	

	file1 = open(filename1+str(j)+".txt", "w")
	file2 = open(filename2+str(j)+".txt", "w")
	theta = -np.pi/3
	e0 = e0array[j-1]	
	
	while theta <= 2*np.pi/3:
	
		k = e0/(h*vF)
		#define the values of the limits for the root finding function.
		valuea = 0.3*k
		valueb = 1.5*k
		#next use bisect on the Rotation value function to find the k values
	
		#print(theta, Dispersion(valuea,args = ([a, t, vF, 1, theta, e0])), Dispersion(valueb,args = ([a, t, vF, -1, theta, e0]))) 
		kx = k*np.cos(theta)
		ky = k*np.sin(theta)
	
		#We use bisect to determine the values of kplus around the fermi surface
		kplus = so.bisect(Dispersion, valuea, valueb, args = ([a, t, vF, 1, theta, e0]))
		#we then repeat the procedure for the negative spin values	
		kminus = so.bisect(Dispersion, valuea, valueb, args = ([a, t, vF, -1, theta, e0]))
		###
		#Next we determine what the x and y positions are. 
		kxp1, kyp1, kxm1, kym1 = kplus*np.cos(theta)/k, kplus*np.sin(theta)/k, kminus*np.cos(theta)/k, kminus*np.sin(theta)/k
	
		#then we print out the particular values of kplus and kminus to a file. 
		ustring = theta, kminus, kplus
		file1.write(' '.join(map(str, ustring)))
		file1.write('\n')
	
		#And following this, we create a string for the kx and ky values. 
		wstring = kxp1, kyp1, kxm1, kym1
		file2.write(' '.join(map(str, wstring)))
		file2.write('\n')	 
	
	
		theta = theta + thetaspacing
	file1.close()
	file2.close()
	print(j)
	j=j +1
	

#########
#########
#########
#So in this section, I'll recalculate the plot, but with the condition on each data set. 

#First, load up the dataset.
j=1
while j <= jmax:
	d["kpx"+str(j)], d["kpy"+str(j)], d["kmx"+str(j)], d["kmy"+str(j)] = loadtxt(filename2+str(j)+".txt", unpack=True)
	d["theta"+str(j)], d["kminus"+str(j)], d["kplus"+str(j)] = loadtxt(filename1+str(j)+".txt", unpack=True)
	j= j + 1

"""
Start of part2

This is the end of part 1 of the code. I now have a set of arrays, with different values. 
In the next section, will construct a random set with different values of guiding center 
co-ordinate.

These are then formed into tuples, and copied into an array of the same length as that of the 
the set of trajectory tuples. Then all() is used to select out the elements such that
"""
j=0
i=0

#define the fixed parameters, focusing length, and size of the pinhole collimator.
l = 2000*10**(-9)
L = 1000*10**(-9)
w = 200*10**(-9)



Barray = np.array(range(-500, 500, 1))*(h*k/(e*l))*(1/800) + (h*k)/(e*l)
emptyp = []
emptym = []
emptyjp = []
emptyjm = []

for i in range(len(Barray)):
#here I define the array of magnetic field values I will check. 
	rescaledtrajectoryxp = np.array(d["kpx"+str(q)])*(h*k)/(e*Barray[i])
	rescaledtrajectoryyp = np.array(d["kpy"+str(q)])*(h*k)/(e*Barray[i])
	rescaledtrajectoryxm = np.array(d["kmx"+str(q)])*(h*k)/(e*Barray[i])
	rescaledtrajectoryym = np.array(d["kmy"+str(q)])*(h*k)/(e*Barray[i])
	delta = (h*k)/(e*Barray[i])*thetaspacing



	scale = 2*l/4
	gccoordx = np.array(np.random.random_sample(500)*scale) - l/4
	gccoordy = np.array(np.random.random_sample(500)*scale) - l/4

	if i == 500:
		randomstartx = gccoordx
		randomstarty = gccoordy
	a = 0 
	vp = 0
	vm = 0 
	j=0
	for j in range(len(gccoordx)):

		gcx = gccoordx[j]
		gcy = gccoordy[j]

		randomfocusingxp = rescaledtrajectoryxp + gcx
		randomfocusingyp = rescaledtrajectoryyp + gcy
		randomfocusingxm = rescaledtrajectoryxm + gcx
		randomfocusingym = rescaledtrajectoryym + gcy

		vpp, jp = Arraydiscarder(randomfocusingxp, randomfocusingyp, args=([L, l, w, delta]))
		vp  = vp + vpp
		vmm, jm = Arraydiscarder(randomfocusingxm, randomfocusingym, args=([L, l, w, delta]))
		vm  = vm + vmm
		if i == 500:
			emptyjp.append(jp)
			emptyjm.append(jm)
		j = j + 1
	emptyp.append(vp)
	emptym.append(vm)
	
	i = i+ 1
	
	
emptypsmooth = running_mean(emptyp, 61)
emptymsmooth = running_mean(emptym, 61)

print(len(emptyp), len(emptypsmooth))	

Barrayplot = Barray[0+30:len(Barray) - 30]
Barraylabel = np.array(range(-1, 2, 1))*(h*k/(e*l))*(1/3) + (h*k)/(e*l)
Barraylabel1 = np.array(range(-1, 2, 1))*(h*np.sqrt(np.pi*narray[q-1])/(e*l))*(1/3) + (h*np.sqrt(np.pi*narray[q-1]))/(e*l)


#here I'll define the sum of the two smoothed values, and extract the maximum. 
sumarray = emptypsmooth + emptymsmooth
maxvalue = np.max(sumarray)
print(maxvalue)

print(len(Barrayplot))
figure(1, figsize=(5, 5))


plt.plot(Barrayplot, emptypsmooth, lw=1, color="red", dashes = [2, 2])
plt.plot(Barrayplot, emptymsmooth, lw=1, color="blue", dashes = [2, 2])
plt.plot(Barrayplot, emptypsmooth + emptymsmooth, lw=1, color="green")


Bfieldvalue = np.round(Barraylabel1, 2)
print(Bfieldvalue)


y = []
labely = y
subplots_adjust(left=0.22, bottom=0.22, right=None, top=None, wspace=None, hspace=None)
yticks(y, labely, fontname="Times New Roman", fontsize=20)
grid(True)
x = Barraylabel
label = Bfieldvalue
xticks(x, label, fontname="Times New Roman", fontsize=20)
ylabel('Intensity', fontname="Times New Roman", fontsize=30)
xlabel('B(T)', fontname="Times New Roman", fontsize=30)
#ylim([-1.5, 1.5])
plt.text(Barrayplot[20], maxvalue, narraylegend[q-1], fontsize=15)
xlim([Barrayplot[0], Barrayplot[len(Barrayplot) - 1]])
#legend((r'$$', r'$x_2$'), prop=FontProperties(size=16))
savefig(figfile1+'.pdf', dpi=100)


###
###
###
#Here I print a plot of the various trajectories, which will be normalised by the scale. 

figure(2, figsize=(5, 5))

j=0
for j in range(len(randomstartx)):
	rescaledtrajectoryxp = (np.array(d["kpx"+str(q)])*(h*k)/(e*Barray[500]) + randomstartx[j])/l
	rescaledtrajectoryyp = (np.array(d["kpy"+str(q)])*(h*k)/(e*Barray[500]) + randomstarty[j])/l
	rescaledtrajectoryxm = (np.array(d["kmx"+str(q)])*(h*k)/(e*Barray[500]) + randomstartx[j])/l
	rescaledtrajectoryym = (np.array(d["kmy"+str(q)])*(h*k)/(e*Barray[500]) + randomstarty[j])/l
	
	xp = rescaledtrajectoryxp[0:emptyjp[j]]
	yp = rescaledtrajectoryyp[0:emptyjp[j]]
	xm = rescaledtrajectoryxm[0:emptyjm[j]]
	ym = rescaledtrajectoryym[0:emptyjm[j]]
	
	
	plt.plot(xp, yp, lw=0.5, color="red", alpha=0.3)
	plt.plot(xm, ym, lw=0.5, color="blue", alpha=0.3)
	
	j = j+1
	
x1, y1 = [0.5, 1-0.125/2], [-0.25 - np.pi/100, -0.25 - np.pi/100]
x2, y2 = [1+0.125/2, 1.3], [-0.25 - np.pi/100, -0.25 - np.pi/100]
x3, y3 = [0.8, 1-0.125/2], [0.25 - np.pi/100, 0.25 - np.pi/100]
x4, y4 = [1+0.125/2, 1.2], [0.25 - np.pi/100, 0.25 - np.pi/100]
x5, y5 = [0.8, 0.8], [-0.25 - np.pi/100, +0.25 - np.pi/100]
x6, y6 = [1.2, 1.2], [-0.25 - np.pi/100, +0.25 - np.pi/100]
x7, y7 = [ np.pi/100, np.pi/100], [0.7, 1-0.125/2]
x8, y8 = [ np.pi/100, np.pi/100], [1+0.125/2, 1.3]

plt.plot(x1, y1, 'k', marker = 'o')
plt.plot(x2, y2, 'k', marker = 'o')
plt.plot(x3, y3, 'k', marker = 'o')
plt.plot(x4, y4, 'k', marker = 'o')
plt.plot(x5, y5, 'k', marker = 'o')
plt.plot(x6, y6, 'k', marker = 'o')
plt.plot(x7, y7, 'k', marker = 'o')
plt.plot(x8, y8, 'k', marker = 'o')

y = [-1, 0, 1]
labely = [-1, " ", "l"]
subplots_adjust(left=0.16, bottom=0.16, right=None, top=None, wspace=None, hspace=None)
yticks(y, labely, fontname="Times New Roman", fontsize=20, fontstyle="italic")
grid(True)
x = [-1, 0, 1]
label = [-1, " ", "l"]
xticks(x, label, fontname="Times New Roman", fontsize=20, fontstyle="italic")
ylabel('y', fontname="Times New Roman", fontsize=30, fontstyle="italic", rotation=90)
xlabel('x', fontname="Times New Roman", fontsize=30, fontstyle="italic")
ylim([-0.4, 1.3])
xlim([0,1.3])
plt.text(0.6, 1.2, narraylegend[q-1], fontsize=15)

#legend((r'$$', r'$x_2$'), prop=FontProperties(size=16))
savefig(figfile2+'.pdf', dpi=100)


