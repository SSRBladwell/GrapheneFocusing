#Code for Paper 3, to finalise a few of the plots. This will calculate 
#trajectories, and the spin orientations for electrons and hole systems.

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

#Introduce the dictionary Bitches!!!

d = {} 

jmax = 5
#Call the filenames I will use for the two valley states. 
filename1 = "kandtheta"
filename2 = "positionsvalues"
filename3 = "Kspace"
filename4 = "Valley1"
filename5 = "startingpositions"
filename6 = "Valley2"

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
k = np.sqrt(2*np.pi*n)
print(k)

#Define the location of K and K' gamma points
Kx = 2*np.pi/(3*a)
Ky = 2*np.pi/(np.sqrt(3)*3*a)

#determining the spacing of the points that I calculate the dispersion at
thetaspacing = np.pi/100.

def Dispersion(k, args):
	a, t, vF, s, theta, e0 = args 
	
	kx = k*np.cos(theta)
	ky = k*np.sin(theta)
	fkk = 2*np.cos(np.sqrt(3)*(-s*ky*a) + 2*np.pi/3) + 4*np.cos(np.sqrt(3)*(-s*ky*a/2)+ np.pi/3)*np.cos(3*(kx)*a/2 + np.pi)
	Dis = t*np.sqrt(3+fkk) - e0
	
	return Dis

j=1
while j <= jmax:

	
	file1 = open(filename1+str(j)+".txt", "w")
	file2 = open(filename2+str(j)+".txt", "w")
	theta = -np.pi/3
	e0 = 0.5*(j)*e - 0.4*e

	while theta <= 2.0000001*np.pi:
	
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
	d["kxp"+str(j)], d["kyp"+str(j)], d["kxm"+str(j)], d["kym"+str(j)] = loadtxt(filename2+str(j)+".txt", unpack=True)
	d["theta"+str(j)], d["kminus"+str(j)], d["kplus"+str(j)] = loadtxt(filename1+str(j)+".txt", unpack=True)
	j= j + 1

#####Next, we need to determine the values for each curve, such that dkx = 0. 

##First, define what the total number of elements of the array is
noofelements = int((np.pi + 2*np.pi)/thetaspacing)

##next, we determine the i value, such that kx changes sign.

file = open(filename1+str(j)+".txt", "w")

j=1
while j <= jmax:
	#define the relevant strings. 
	print(j)
	kxp = d["kxp"+str(j)]
	kxm = d["kxm"+str(j)]
	kyp = d["kyp"+str(j)]
	kym = d["kym"+str(j)]	
	

	for i in range(0, 100):
		if i == 0:
			sign = np.sign(kxp[i+1] - kxp[i])
			print(sign)
		diff1 = sign*(kxp[i+1] - kxp[i])
		if diff1 <= 0:
			ip = i
			break

	for i in range(0, 100):
		if i == 0:
			sign = np.sign(kxm[i+1] - kxm[i])
			print(sign)
		diff1 = sign*(kxm[i+1] - kxm[i])
		if diff1 <= 0:
			im = i
			break


	d["kxpinj"+str(j)] = kxp - kxp[ip]
	d["kypinj"+str(j)] = kyp - kyp[ip]
	d["kxminj"+str(j)] = kxm - kxm[im]
	d["kyminj"+str(j)] = kym - kym[im]
	
	
	
	
	
	
	for i in range(ip + 3, 200):
		u = d["kxpinj"+str(j)][i] + d["kxp"+str(j)][ip] 
		if u <= 0.:
			icutp = i
			break
			
	for i in range(im + 3, 200):
		u = d["kxminj"+str(j)][i] + d["kxm"+str(j)][im] 
		if u <= 0.:
			icutm = i
			break
	print(icutm)
	
	print(j)
	d["kxpcut"+str(j)] = d["kxpinj"+str(j)][ip:icutp]
	d["kxmcut"+str(j)] = d["kxminj"+str(j)][im:icutm]
	d["kypcut"+str(j)] = d["kypinj"+str(j)][ip:icutp]
	d["kymcut"+str(j)] = d["kyminj"+str(j)][im:icutm]
	j = j + 1
j =1
	



########
########
########


#Next I form the distribution of injection angles. 

#First, I define the injection angle, in terms of the polar angle
#theta0p = theta[ip1]
#theta0m = theta[im1]

#print(theta0p, theta0m)

#Next, let us determine the number of counts for a particular value of theta

#COUNTARRAYp = list(range(40))
#COUNTARRAYm = list(range(40))

#thetaspread = np.pi/20


#for j in range(0, 40):
#	COUNTARRAYp[j] = np.floor(20*np.exp(-(theta1[j - 20 + ip1])**2/thetaspread))
#	COUNTARRAYm[j] = np.floor(20*np.exp(-(theta1[j - 20 + im1])**2/thetaspread))

#next up I create a sum of the trajectories, but counting the number of states that
#arrive at the collector. 

#print(kyp1[ip1], kxp1[ip1])

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#First, let's plot the Electron Trajectories. 
figure(2, figsize=(5, 5))

#plt.quiver(kxp, kyp, sigmaxp, sigmayp, lw=0.1, color='green')
#plt.quiver(kxm, kym, sigmaxm, sigmaym, lw=0.1, color='green')

plt.plot(d["kxp1"], d["kypinj"+str(1)], lw=1, color="red")
plt.plot(d["kxm1"], d["kyminj"+str(1)], lw=1, color="blue")
	
y = [-1.5, 0, 1.5]
labely = y
subplots_adjust(left=0.22, bottom=0.22, right=None, top=None, wspace=None, hspace=None)
yticks(y, labely, fontname="Times New Roman", fontsize=20)
grid(True)
x = [-1.5, 0, 1.5]
label = x
xticks(x, label, fontname="Times New Roman", fontsize=20)
xlabel('$y/r_c$', fontname="Times New Roman", fontsize=30)
ylabel('$x/r_c$', fontname="Times New Roman", fontsize=30, rotation=90)
ylim([-1.5, 1.5])
xlim([-1.5, 1.5])
#legend((r'$$', r'$x_2$'), prop=FontProperties(size=16))
savefig(filename3+'.pdf', dpi=100)

figure(3, figsize=(5, 5))


plt.plot(d["kxp4"], d["kypinj"+str(4)], lw=1, color="red")
plt.plot(d["kxm4"], d["kyminj"+str(4)], lw=1, color="blue")

y = [-1.5, 0, 1.5]
labely = y
subplots_adjust(left=0.22, bottom=0.22, right=None, top=None, wspace=None, hspace=None)
yticks(y, labely, fontname="Times New Roman", fontsize=20)
grid(True)
x = [-1.5, 0, 1.5]
label = x
xticks(x, label, fontname="Times New Roman", fontsize=20)
xlabel('$y/r_c$', fontname="Times New Roman", fontsize=30)
ylabel('$x/r_c$', fontname="Times New Roman", fontsize=30, rotation=90)
ylim([-1.5, 1.5])
xlim([-1.5, 1.5])
#legend((r'$$', r'$x_2$'), prop=FontProperties(size=16))
savefig(filename3+'1.pdf', dpi=100)



#First, let's plot the Electron Trajectories. 


j=1
while j < jmax:
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	
	figure(4 + j, figsize=(5, 5))

	plt.plot(d["kxp"+str(j)], d["kyp"+str(j)], lw=1, color="red")

	y = [-1, 0, 1]
	labely = y
	subplots_adjust(left=0.25, bottom=0.25, right=None, top=None, wspace=None, hspace=None)
	yticks(y, labely, fontname="Times New Roman", fontsize=40)
	grid(True)
	x = [-1, 0, 1]
	label = ["$-1$", "$0$","$1$"]
	xticks(x, label, fontname="Times New Roman", fontsize=40)
	xlabel('$k_x/k$', fontname="Times New Roman", fontsize=40)
	ylabel('$k_y/k$', fontname="Times New Roman", fontsize=50, rotation=90)
	ylim([-1.5, 1.5])
	xlim([-1.5, 1.5])
	#legend((r'$$', r'$x_2$'), prop=FontProperties(size=16))
	savefig(filename4+str(j)+'.pdf', dpi=100)

	
	
	figure(9 + j, figsize=(5, 5))

	plt.plot(d["kxm"+str(j)], d["kym"+str(j)], lw=1, color="blue")

	y = [-1, 0, 1]
	labely = y
	subplots_adjust(left=0.25, bottom=0.25, right=None, top=None, wspace=None, hspace=None)
	yticks(y, labely, fontname="Times New Roman", fontsize=40)
	grid(True)
	x = [-1, 0, 1]
	label = x
	xticks(x, label, fontname="Times New Roman", fontsize=40)
	xlabel("$k_x/k$",  fontname="Times New Roman", fontsize=50)
	ylabel("$k_y/k$",  family='Times New Roman', fontsize=50, rotation=90)
	ylim([-1.5, 1.5])
	xlim([-1.5, 1.5])
	#legend((r'$$', r'$x_2$'), prop=FontProperties(size=16))
	savefig(filename6+str(j)+'.pdf', dpi=100)
	j= j+ 1






#Ok, next I want to do a spread plot, for all the possible trajectories. 
#Spacing of the spread of the beam. 




figure(20, figsize=(4, 5))


plt.plot(d["kxpinj4"], d["kypinj4"], lw=1, color="red")
plt.plot(d["kxminj4"], d["kyminj4"], lw=1, color="blue")

y = [-1.5, 0, 1.5]
labely = y
subplots_adjust(left=0.22, bottom=0.22, right=None, top=None, wspace=None, hspace=None)
yticks(y, labely, fontname="Times New Roman", fontsize=20)
grid(True)
x = [-1.5, 0, 1.5]
label = x
xticks(x, label, fontname="Times New Roman", fontsize=20)
xlabel('$y/r_c$', fontname="Times New Roman", fontsize=30)
ylabel('$x/r_c$', fontname="Times New Roman", fontsize=30, rotation=90)
ylim([0, 1.5])
xlim([-1, 0.2])
#legend((r'$$', r'$x_2$'), prop=FontProperties(size=16))
savefig(filename3+'4.pdf', dpi=100)












