

'''
==================================================================================================
Script to write to file the complex refractive index for 30 wavelengths within the visible range
having in input the dielectric function
==================================================================================================

sys.argv[1] = compound name
sys.argv[2] = file name with the dielectric function
'''



import os , json , sys
import numpy as np
import colour
import matplotlib.pyplot as plt
import matplotlib


# Input from command-line
compound_name = sys.argv[1]
eps_filename = sys.argv[2]

## Here I take the dielectric function from the outputfile of Yambo
data = np.genfromtxt(eps_filename)
energies  = data[:,0]
eps_im = data[:,1]
eps_re = data[:,2]

# n(E) and k(E)
norm_epsilon = np.sqrt(eps_re**2 + eps_im**2)
refractive_index = np.sqrt( ( eps_re + norm_epsilon  ) / 2. )
ext_coeff = np.sqrt( ( -eps_re + norm_epsilon  ) / 2. )

# extract the visible part of the complex refractive index
wavelenght_nm = 1239.8 / energies
wavelenght_visible = []
refractive_index_visible = []
ext_coeff_visible = []
for i in xrange(len(wavelenght_nm)):
	if  370.0 <= wavelenght_nm[i] <= 820.0:
            wavelenght_visible.append(wavelenght_nm[i])            
            refractive_index_visible.append(refractive_index[i])
            ext_coeff_visible.append(ext_coeff[i])
        
# Polynomial fit of n(lambda) and k(lambda)
fit_params_refractive = np.polyfit(wavelenght_visible, refractive_index_visible, 11)
polynomial_refractive = np.poly1d(fit_params_refractive) 

fit_params_ext = np.polyfit(wavelenght_visible, ext_coeff_visible, 11)
polynomial_ext = np.poly1d(fit_params_ext)

# Plot fit  
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.figure(figsize=(15,10))
plt.plot(wavelenght_visible,refractive_index_visible, color='black', linewidth=1.0, label='Data')
plt.plot(wavelenght_visible, polynomial_refractive(wavelenght_visible), color='blue', linewidth=1.0, label='Fit')
plt.xlabel(r'Wavelenght [nm]',size=14)
plt.ylabel(r'n($\boldsymbol{\lambda}$)',size=14)
#plt.xlim(380.0,780.0)
plt.title(r'Fit of the refractive index', size=19)
plt.legend()
plt.show()    

matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.figure(figsize=(15,10))
plt.plot(wavelenght_visible,ext_coeff_visible, color='black', linewidth=1.0, label='Data')
plt.plot(wavelenght_visible, polynomial_ext(wavelenght_visible), color='blue', linewidth=1.0, label='Fit')
plt.xlabel(r'Wavelenght [nm]',size=14)
plt.ylabel(r'k($\boldsymbol{\lambda}$)',size=14)
#plt.xlim(380.0,780.0)
plt.title(r'Fit of the extinction coefficient', size=19)
plt.legend()
plt.show()    

## Write n(lambda) and k(lambda) to file (30 data points)
with open('{}.eta.dat'.format(compound_name),'w') as o: 
	for i in xrange(370,820,15):	
		o.write(str(i))
		o.write(' '+str(polynomial_refractive(i)))
		o.write('\n')
with open('{}.k.dat'.format(compound_name),'w') as o: 
	for i in xrange(370,820,15):
		o.write(str(i))
		o.write(' '+str(polynomial_ext(i)))
		o.write('\n')




