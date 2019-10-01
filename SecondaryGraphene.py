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

jmax = 35
#Call the filenames I will use for the two valley states. 
filename1 = "kandtheta"
filename2 = "positionsvalues"
filename3 = "Splitting"
filename4 = "SplittingPicture"

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
thetaspacing = np.pi/400.

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
	e0 = 0.025*(j)*e

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
		ustring = theta, kminus, kplus, k
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
	j= j + 1

#####Next, we need to determine the values for each curve, such that dkx = 0. 

##First, define what the total number of elements of the array is
noofelements = int((np.pi + 2*np.pi)/thetaspacing)

##next, we determine the i value, such that kx changes sign.

file3 = open(filename3+".txt", "w")

j=1
while j <= jmax:
	#define the relevant strings. 
	print(j)
	kxp = d["kxp"+str(j)]
	kxm = d["kxm"+str(j)]
	kyp = d["kyp"+str(j)]
	kym = d["kym"+str(j)]	
	n = (0.024*(j)*e/(h*vF))**2/np.pi
	

	for i in range(0, 1000):
		if i == 0:
			sign = np.sign(kxp[i+1] - kxp[i])
			print(sign)
		diff1 = sign*(kxp[i+1] - kxp[i])
		if diff1 <= 0:
			ip = i
			break
	
	for i in range(0, 1000):
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
	
	
	
	for i in range(ip + 3, 2000):
		if d["kxpinj"+str(j)][i] + 1 <= 0:
			icutp = i
			break
			
	for i in range(im + 3, 2000):
		u = d["kxminj"+str(j)][i] + 1
		if u <= 0.:
			icutm = i
			break

	
	
	print(j)
	y1 = d["kypinj"+str(j)][icutp]
	y2 = d["kyminj"+str(j)][icutm]
	splity = y1 - y2
	stringy = [splity, n]
	file3.write(" ".join(map(str,stringy)))
	file3.write('\n')	
	j = j + 1
j =1
file3.close()



########
########
########
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#Next, I load the text file with the results
splitting, n = loadtxt(filename3+".txt", unpack=True)

figure(20, figsize=(8, 6))


plt.plot(n, splitting, lw=1.5, color="k")
plt.axhline(0.3/1.25, lw=0.5, color='g')
plt.axvline(2.7*10**(16), lw =0.5, color='b')



y = [0, 0.25, 0.5]
labely = y
subplots_adjust(left=0.22, bottom=0.22, right=None, top=None, wspace=None, hspace=None)
yticks(y, labely, fontname="Times New Roman", fontsize=20)
grid(True)
x = [0, 2*10**(17), 4*10**(17)]
label = ["$0$", "$2 \\times 10^{13}$", "$4 \\times 10^{13}$" ]
xticks(x, label, fontname="Times New Roman", fontsize=20)
xlabel('$n$ (cm$^{-2}$)', fontname="Times New Roman", fontsize=30)
ylabel('$\Delta x/r_c$', fontname="Times New Roman", fontsize=30, rotation=90)
#legend((r'$$', r'$x_2$'), prop=FontProperties(size=16))
xlim([0, 6*10**(17)])
savefig(filename3+'4.pdf', dpi=100)

