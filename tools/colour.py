# -*- coding: utf-8 -*-

import numpy as np
import os , math
import matplotlib.pyplot as plt
import matplotlib

'''
Python module to calculate reflectivity and CIE colours from the knowledge of the complex dielectric function
'''

def calcReflectivity(epsilon):
    '''
    Calculate reflectivity (at normal incidence) from the dielectric function.
    
    :param epsilon: complex numpy array with the dielectric function (real_part + imaginary_part*j)
    :return reflectivity: numpy array with the reflectivity
    '''
    
    epsilon_re = epsilon.real
    epsilon_im = epsilon.imag
    norm_epsilon = np.sqrt(epsilon_re**2 + epsilon_im**2)

    refractive_index = np.sqrt( ( epsilon_re + norm_epsilon  ) / 2. )
    abs_coeff = np.sqrt( ( -epsilon_re + norm_epsilon  ) / 2. )
    reflectivity = ( (refractive_index - 1)**2 + abs_coeff**2 ) / ( (refractive_index + 1)**2 + abs_coeff**2 )
    
    return reflectivity
    


def calcColour(energy_ev, reflectivity, file_d65illuminant, file_cmf, wavelenght_steps = '1nm', do_plot=True):
    '''
    Calculate colour coordinates from reflectivity.
    
    :param energy_ev: array with energy steps (must be in eV)
    :param reflectivity: array with reflectivity values corresponding to 'energy_ev'
    :param file_d65illuminant: file with CIE standard illuminant D65. Must contain two columns: 
        wavelenghts + D65 data
    :param file_cmf: file with Color Matching Functions (CMFs). Must contain four columns:
        wavelenghts + CMF_x + CMF_y + CMF_z
    :param wavelenght_steps: can be '1nm' (default) or '5nm'. These are the supported wavelenght steps
        for data in 'file_d65illuminant' and 'file_cmf'.
    :return a dictionary of the form
        {'Tristimulus values': CIE-XYZ, 
        'Chromaticity coordinates': CIE-xyY,
        'CIELAB': CIE-L*a*b*,
        'Yellowness index': D1925 yellowness index,
        'sRGB': standard RGB,
        'HEX': hexadecimals,
        }         
    '''
    
    data = np.genfromtxt(file_d65illuminant)
    d65_illuminant  = data[:,1]    
    data = np.genfromtxt(file_cmf)
    cmf_x  = data[:,1]
    cmf_y  = data[:,2]
    cmf_z  = data[:,3]
    
    if len(d65_illuminant) != len(cmf_x):
        raise Exception("D65 and CMFs data do not have the same lenght")
    
    if len(energy_ev) != len(reflectivity):
        raise Exception("'energy_ev' and 'reflectivity' do not have the same lenght")
    
    
    wavelenght_nm = 1239.8 / energy_ev
    
    # extract the visible part of the reflectivity curve
    wavelenght_visible = []
    reflectivity_visible = []
    for i in xrange(len(wavelenght_nm)):
        if  380.0 <= wavelenght_nm[i] <= 780.0:
            wavelenght_visible.append(wavelenght_nm[i])            
            reflectivity_visible.append(reflectivity[i])
        
    
    # Polynomial fit of reflectivity_visible
    fit_params, fit_residuals, _, __, ___ = np.polyfit(wavelenght_visible, reflectivity_visible, 11,full=True)
    polynomial = np.poly1d(fit_params) 
    
    if do_plot:
      matplotlib.rc('text', usetex=True)
      matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
      plt.figure(figsize=(15,10))
      plt.plot(wavelenght_visible,reflectivity_visible, color='black', linewidth=1.0, label='Data')
      plt.plot(wavelenght_visible, polynomial(wavelenght_visible), color='blue', linewidth=1.0, label='Fit')
      plt.xlabel(r'Wavelenght [nm]',size=14)
      plt.ylabel(r'R($\boldsymbol{\lambda}$)',size=14)
      #plt.xlim(380.0,780.0)
      plt.title(r'Fit of the reflectivity', size=19)
      plt.legend()
      plt.show()    
    
    wavelenght_uniform=[]
    reflectivity_fit=[]
    if wavelenght_steps == '1nm':
        for i in xrange(380,781):
            wavelenght_uniform.append(i)
            reflectivity_fit.append(polynomial(i))
    elif wavelenght_steps == '5nm':
        for i in xrange(380,781,5):
            wavelenght_uniform.append(i)
            reflectivity_fit.append(polynomial(i))
    else:
        raise Exception("Supported 'wavelenght_steps' are only '1nm' and '5nm'")    
        
    wavelenght_uniform=np.array(wavelenght_uniform)
    reflectivity_fit=np.array(reflectivity_fit)
   
    
    if len(reflectivity_fit) != len(d65_illuminant):
       raise Exception("Fit of reflectivity in the visible range do not have the same lenght of D65 and CMFs data")  
       
    # Calculate renormalization constant k
    k = 100.0 / sum(cmf_y*d65_illuminant)
    
    # Calculate Xn, Yn and Zn  (Yn=100)
    x_n = k*sum(cmf_x*d65_illuminant)
    y_n = k*sum(cmf_y*d65_illuminant)
    z_n = k*sum(cmf_z*d65_illuminant)
        
    # Tristimulus values
    x = k * sum(cmf_x*reflectivity_fit*d65_illuminant)
    y = k * sum(cmf_y*reflectivity_fit*d65_illuminant)
    z = k * sum(cmf_z*reflectivity_fit*d65_illuminant)
     
    # D1925 yellowness index (for D65 illuminant and 1931 CIE observer)
    d1925_index = 100.0*( 1.2985*x-1.1335*z )/y 
    
    # Chromaticity coordinates
    x_chrom = x / (x+y+z)
    y_chrom = y / (x+y+z)
    z_chrom = z / (x+y+z) 
    
    # Calculate L,a,b of CIELAB
    def f(x):
        if x > (24./116.)**3:
           return x**(1./3.)
        elif x <= (24./116.)**3:
           return (841./108.)*x + (16./116.)   
        
    l = 116.*f(y/y_n) - 16.
    a = 500.*(f(x/x_n) - f(y/y_n)  )
    b = 200.*(f(y/y_n) - f(z/z_n)  )
        
    # Calculate chroma and hue
    chroma = math.sqrt(a**2 + b**2)
    hue = math.degrees(math.atan(b/a))


    # Standard RGB (sRGB)
    var_X = x / 100.
    var_Y = y / 100.
    var_Z = z / 100.

    var_R = var_X *  3.2406 + var_Y *(-1.5372) + var_Z * (-0.4986)
    var_G = var_X * -0.9689 + var_Y *  1.8758 + var_Z *  0.0415
    var_B = var_X *  0.0557 + var_Y *(-0.2040) + var_Z *  1.0570

    if  var_R > 0.0031308:
        var_R = 1.055 * ( var_R**( 1. / 2.4 ) ) - 0.055
    else:                     
        var_R = 12.92 * var_R
    if var_G > 0.0031308: 
        var_G = 1.055 * ( var_G**( 1. / 2.4 ) ) - 0.055
    else:                     
        var_G = 12.92 * var_G
    if var_B > 0.0031308: 
        var_B = 1.055 * ( var_B**( 1. / 2.4 ) ) - 0.055
    else:                     
        var_B = 12.92 * var_B

    sR = var_R * 255.
    sG = var_G * 255.
    sB = var_B * 255.

    # Hexadecimal colour 
    def clamp(x): 
        return max(0, min(x, 255))  # to ensure that  0 < sR,sG,sB < 255

    hexadecimal = '#%02x%02x%02x' % (clamp(sR), clamp(sG), clamp(sB))

    return {'Colours':
                {'Tristimulus values': {'X': x, 'Y': y, 'Z': z}, 
                'Chromaticity coordinates': {'x':x_chrom, 'y':y_chrom, 'z':z_chrom,},
                'CIELAB': {'L':l, 'a': a, 'b': b, 'Chroma': chroma, 'Hue': hue},
                'Yellowness index': d1925_index,
                'sRGB': {'R':sR, 'G':sG, 'B':sB},
                'HEX': hexadecimal,
                },
            'Fit_residuals': fit_residuals,
            }

