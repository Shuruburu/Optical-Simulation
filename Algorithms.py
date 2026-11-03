#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 16:52:07 2026

@author: shurubura
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sci 
c = 3*1e8

hbar = 1
from pathlib import Path

base_dir = Path(__file__).resolve().parent
save_dir  = base_dir / "plots"
n_quartz = 1.55 

"""
 the funny n letter tellus us , but we have internal losses 
 
"""

from matplotlib.lines import Line2D  


def find_collimation(z_vals, q0, g):
    for z in z_vals:
        qz = propagate_q_grin(q0, z, g)
        if np.isclose((1/qz).real, 0, atol=1e-6):
            return z, qz
    return None, None


def Calculate_the_length_of_the_mirror(x , c):
    y1 = x-np.sqrt(x**2 - c**2)
    y2 = x+np.sqrt(x**2 - c**2)
    return  y1 , y2











def grin_abcd(z, g):
    """
    Returns ABCD matrix elements for a GRIN medium.
    
    Parameters:
        z : propagation distance
        g : gradient parameter
        
    Returns:
        A, B, C, D
    """
    A = np.cos(g * z)
    B = (1 / g) * np.sin(g * z)
    C = -g * np.sin(g * z)
    D = np.cos(g * z)
    return A, B, C, D


def propagate_q_grin(q0, z, g):
    """
    Propagate Gaussian beam q-parameter through GRIN medium.
    
    Parameters:
        q0 : complex initial q parameter
        z  : propagation distance
        g  : gradient parameter
        
    Returns:
        q(z)
    """
    A, B, C, D = grin_abcd(z, g)
    return (A * q0 + B) / (C * q0 + D)


def beam_parameters(q, wavelength):
    """
    Extract beam radius w and curvature R from q.
    
    Parameters:
        q : complex beam parameter
        wavelength : wavelength
        
    Returns:
        w : beam radius
        R : wavefront curvature (np.inf if collimated)
    """
    inv_q = 1 / q
    
    # Curvature
    if np.isclose(inv_q.real, 0):
        R = np.inf
    else:
        R = 1 / inv_q.real
    
    # Beam radius
    w = np.sqrt(-wavelength / (np.pi * inv_q.imag))
    
    return w, R

def Visualize_the_beam():
    # Parameters
    wavelength = 1.55e-6   # 1550 nm
    w0 = 50e-6            # 50 microns
    R0 = -0.01            # diverging beam (negative curvature)
    g = 100               # GRIN parameter (1/m)

    # Initial q
    q0 = 1 / (1/R0 - 1j * wavelength / (np.pi * w0**2))

    # Propagate
    z_vals = np.linspace(0, 0.05, 500)

    w_vals = []
    R_vals = []

    for z in z_vals:
        qz = propagate_q_grin(q0, z, g)
        w, R = beam_parameters(qz, wavelength)
        w_vals.append(float(w))
        R_vals.append(float(R))
    


def Fixing_the_position(Ls ,lam , z2 ,R1 ):
    Cavity_Modes = ((Ls * 2) / lam).astype(int)
    R2 = np.empty(R1.shape)
    waist = np.empty(R1.shape)
    dict_waists = {}
    dict_R2 = {}
    for L in Ls:
        for  n , R  in enumerate(R1):
             R2[n] , waist[n] = cavity_fix_z2(z2 , R , L , Cavity_Modes )
             
        print(min(R2))
        print(min(waist))
        dict_waists[L] =  waist*1e6
        dict_R2[L]  = R2*1e3
    Plot_dict(dict_waists, R1*1e3 ,   "Waist[micrometers]" , "R1(Curvuture of the first mirror)[mm]" , title  = "Curvuture pf t")
    Plot_dict(dict_R2 , R1*1e3 , "R2[mili]" ,"R1(Curvuture of the first mirror)[mm]" ,  title  = "When while z = -3  to ")
    




def Plot_2D_single(R1, R2, z , zlabel , ):
    """
    R1, R2 : 2D meshgrid arrays
    waist  : single 2D array
    """

    plt.figure()

    plt.imshow(
        z,
        extent=[R1.min()*1e3, R1.max()*1e3, R2.min()*1e3, R2.max()*1e3],
        origin='lower',
        aspect='auto'
    )

    plt.colorbar(label= zlabel)
    plt.title("2D Waist Map")
    plt.xlabel("R1 [mm]")
    plt.ylabel("R2 [mm]")

    plt.tight_layout()
    plt.show()






def Average_waist(R1, R2, L, width, z0, lam):
    """
    z0: position of the HBAR
    width: width of HBAR
    R1 : Curvature of the first mirror
    R2 : Curuvuture of the 
    """
    cavity_mode = ((L * 2) / lam).astype(int)
    zs = np.linspace(z0 - width/2, z0 + width/2, 500)

    waist, z1, z2, zr = Cavity(R1, R2, L, cavity_mode, n=1, Ray=True)
    
    waists_av = Waist_z(waist, zr, zs)

    # Take the average over the sampled region
    average_waist = np.mean(waists_av)

    return average_waist








def Plot_surface(R1, R2, waists, z_axis, mode="surface" , title  =  None):
    color_dict = {
        1: "blue", 2: "red", 3: "green", 4: "orange", 5: "purple",
        6: "black", 7: "yellow", 8: "cyan", 9: "magenta", 10: "pink",
        11: "brown", 12: "gray", 13: "lime", 14: "teal", 15: "navy",
        16: "maroon", 17: "olive", 18: "violet", 19: "indigo", 20: "gold"
    }

    colors = list(color_dict.values())

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    # Keep track of colors and labels for legend
    legend_elements = []

    if mode == "surface":
        for n, waist in enumerate(waists):
            color = colors[n % len(colors)]
            ax.plot_surface(R1*1e3, R2*1e3, waist*1e6, color=color, alpha=0.7)
            # Create a proxy artist for legend
            legend_elements.append(Line2D([0], [0], color=color, lw=4, label=f"mode {n+1}"))
    elif mode == "scatter":
        for n, waist in enumerate(waists):
            color = colors[n % len(colors)]
            ax.scatter(R1*1e3, R2*1e3, waist*1e6, color=color)
            legend_elements.append(Line2D([0], [0], marker='o', color=color, label=f"Waist {n+1}", linestyle=''))

    ax.set_xlabel('R1 [millimeters]')
    ax.set_ylabel('R2 [millimeters]')
    ax.set_zlabel('waist [micrometers]')

    # Add legend using the proxy elements
    ax.legend(handles=legend_elements, loc='upper right')

    fig.tight_layout()
    if title != None:
        ax.set_title(title)
        plt.savefig(save_dir /  f"{title} - 3D surface plot.jpg")
    plt.show()


def Position_of_sample_variation(L, Lam, R1, R2, z0_arr , width ):  
    # Position of the sample 
    arr_av_waist = []

    # Define R ranges

    # Create grid
    R1_grid, R2_grid = np.meshgrid(R1, R2)
    Rs = np.column_stack((R1_grid.ravel(), R2_grid.ravel()))

    for z0 in z0_arr:
        av_waist = {}  # reset for each L

        for R in Rs:
            R_key = (float(R[0]), float(R[1]))  # make hashable
            temp = Average_waist(R[0], R[1], L, width, z0, Lam)
            av_waist[R_key] = float(temp)
        arr_av_waist.append(av_waist)

    return arr_av_waist


# Need to add just part of the code where I dont simulate all of the R I just go over one R1, inf 

def Simulate_different_R1_R2(Ls, lam  , R1 , R2 , Plot = True):
    
    Cavity_Modes = ((Ls * 2) / lam).astype(int)
    minLs = []
    arr_waists = []
    arr_z2 = []
    Min_geometry =  []
    R1, R2 = np.meshgrid(R1, R2)
    Rs = np.column_stack((R1.ravel(), R2.ravel()))

    for idx, Cavity_mode in enumerate(Cavity_Modes):
        waists = np.empty(len(Rs))
        z2 = np.empty(len(Rs))

        for i, R in enumerate(Rs):
            w, z1, z2_val = Cavity(R[0], R[1], Ls[idx], Cavity_mode, n=1)
            waists[i] = w
            z2[i] = z2_val

        waists = waists.reshape(R1.shape)
      #  print(waists)
        z2 = z2.reshape(R1.shape)
        arr_waists.append(waists)
        arr_z2.append(z2)
        min1 , min_indx = Find_the_min(waists)
        Min_geometry.append(Min_geometry)
        minLs.append(min1)
        print("The smallest possible waist for this lenfgh of the mirror is {}  microns " .format(min1*1e6)) # showing the unities in microns
        print("The geomeretry of the cavity was {}". format(Rs[min_indx]))
        if Plot == True:
            Plot_2D_single(R1, R2, waists*1e6 , 'Waist [µm] ')
            Plot_2D_single(R1, R2,  z2*1e3, 'Position of the backmirror[mm]')
        
    arr_waists = arr_waists
   
    #arr_z2 = np.array(arr_z2)

    #Plot_surface(R1, R2, arr_z2[2], 'z2 [micrometers]', mode="scatter" , )
    if Plot == True:
        Plot_surface(R1, R2, arr_waists , 'waist [micrometers]' , title  = " Surface plot title x")
    
    return arr_waists, arr_z2



def Find_the_min(arr):
    # the function return the lowest value in an array and the index of it in the array
    min1 = float('inf')  # start with very large number
    min_idx =  -1
    for n ,row in enumerate(arr):
        for element in row:
            if element < min1:
                min1 = element
                min_idx = n
    return min1 , n 

def Match_the_modes(acoutic_modes , Optical_modes ):
    """
    Optical modes in radians like the one with two pi factor
    Acoustic mode in radians 
    With this function we want to match the acosutic mode to optical mode , this is to ensure the high possible coupling in terms of lognitudal degree
    The transversal coupling depands on the transversal mode (not that easy)
    """
    matching_modes = {}
    for Acoustic in acoutic_modes:
        for n , optical in enumerate(Optical_modes):
                if abs(Bernouli_frequency(Optical_modes) - optical) < 0.3 :
                    matching_modes[n] = optical
    return matching_modes 


def Radius_(R):
    """
    How does the Hi form DASE, which we do know this we need to improve the coupling to the quibit , what the fuck are the ones
    from the antenna, how do they coupled , 
    """    
    return
def Simulate_different_coupling_rate(photo_tensor,  omega_broilin):
    R_in  = 2    
    return

def Omega_plot(Ls, m):
    Omega_Bs = np.empty((Ls.size -1, 20 ))
    for i in range(0,Ls.size -1 ):
        Omega_Bs[i] =Simulate_acoustic_modes(Ls[i+1],  m  ,n = 1 )

    Plot(Ls[1:]*1e3, Omega_Bs.T*1e-9, "L[milimeter]", 'Omega[GHz]')



def Simulate_acoustic_modes(L,  m  ,vl = 6000 , n =1 ):
    
    return  2*np.pi * vl/(2*L)





def Bernouli_frequency( omega , vl = 6000  ,n =1 ):
    return 2*omega*vl* n / c 



def plot_mode_differences(Ls, c, ms):
    prev_modes = None
    
    for m in ms:
        modes = Optical_modes(Ls[1:-1], c, m)

        if prev_modes is not None:
            diff = modes - prev_modes
            Plot(diff * 1e-9, Ls[1:-1] * 1e3,
                    "ΔOmega [GHz]",
                    "Distance [millimeter]",
                    1)

        prev_modes = modes

def Plot_dict(dictionary, omega, ylabel, xlabel, title=None):
    color_dict = {
        1: "blue",
        2: "red",
        3: "green",
        4: "orange",
        5: "purple",
        6: "black",
        7: "yellow",
        8: "cyan",
        9: "magenta",
        10: "pink",
        11: "brown",
        12: "gray",
        13: "lime",
        14: "teal",
        15: "navy",
        16: "maroon",
        17: "olive",
        18: "violet",
        19: "indigo",
        20: "gold"
    }

    fig, ax = plt.subplots(figsize=(10, 6))

    for n, key in enumerate(dictionary.keys()):
        color_index = (n % len(color_dict)) + 1  # wrap around if more than 20
        ax.plot(omega, dictionary[key], color=color_dict[color_index], label=str(key))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
        plt.savefig(save_dir /  f" {title}.jpg")
    ax.legend()
    ax.grid(True)
    plt.show()
    
def Simulate_coupling(L, loss):
    
    omega0 = Optical_modes(L, c, 1)

    sweep_omega = np.linspace(omega0 - omega0*0.01, omega0 + omega0*0.01, 10000)

    Ref1 = np.linspace(0.98, 0.999, 10)

    results = {}

    for ref in Ref1:
        results[float(ref)] = Coupling_constant(L, ref, 0.9999, loss, sweep_omega, omega0)

    Plot_dict(results ,sweep_omega*1e-9 , "Reflexivity" , "Omega[GHz]")
    
def Coupling_constant(L , Ref1 , Ref2 ,  loss , omega , omega0):
    """
        
    """
    FSR=Optical_modes(L  , c ,1 )
    kappa1 = (1-Ref1) *FSR
    kapppa2 = (1-Ref2) * FSR
    kapppa_int = loss *FSR
    kappa_tot = kappa1 + kapppa2 + kapppa_int
    detuning = omega - omega0
    
    Rc = np.abs(1 - (kappa1 / ( - 1j*detuning  + kappa_tot/2)))**2
    return Rc


def Plot(x, ys, xlabel, ylabel, vertical=None):

    color_dict = {
        1: "blue", 2: "red", 3: "green", 4: "orange", 5: "purple",
        6: "black", 7: "yellow", 8: "cyan", 9: "magenta", 10: "pink",
        11: "brown", 12: "gray", 13: "lime", 14: "teal", 15: "navy",
        16: "maroon", 17: "olive", 18: "violet", 19: "indigo", 20: "gold"
    }

    fig, ax = plt.subplots(figsize=(10, 6))

    for i, y in enumerate(ys):
        ax.scatter(x, y, color=color_dict[i % 20 + 1], label=str(i+1))
        

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if vertical is not None:
        ax.plot(x, np.ones(x.size) * vertical, color=color_dict[2], label="vertical")

    ax.legend()
    ax.grid(True)

    plt.show()


def Raylaigh_range_1(waist, Lam):
    return  waist**2 *  np.pi / Lam


def Radius_curvature(z , zr ):
    return z*(1+  (zr/z)**2)

def Modes(L, m , n =1):
    return  c* n  * m/ (2 * L)

def Optical_modes(L, cn, m):
    return  cn / (2 * L)

def Position_z2(R1, R2,  L ):
    return L*(R1 - L) / (R2 + R1  -  2*L)
 
def Position_z1(z2  , L): 
    
    return z2 + L


def Rayleigh(Lam , waist):
    return Lam / (np.pi * waist**2) 

def Rayleigh_range(z, R ):
    """
    The function returns the Rayleigh range for the specific curvature
    """
    return np.sqrt(R / z - 1) * z

def Waist_of_the_beam(Lam, Ray):
    return np.sqrt(Ray * Lam / (np.pi))


def cavity_fix_z2(z2, R1, L, m, n=1):
    """Fix cavity geometry assuming R1 < R2 and beam waist at z=0."""

    modes = Modes(L, m, n=n)
    if modes == 0:
        raise ValueError("Modes returned zero → invalid frequency")

    Lam = c / modes
    Lam_eff = Lam / n

    zr = Rayleigh_range(z2 + L, R1)
    w0 = Waist_of_the_beam(Lam_eff, zr)

    R2 = Radius_curvature(-z2, zr)  # fixed typo

    return R2, w0

def Waist_z(w0 , zr, z ):
    return w0 *np.sqrt( 1 + (z/zr)**2)


def Cavity(R1, R2, L, m, n=1, Ray = False):
    """
    Calculate Gaussian beam waist and mirror positions for a cavity.
    R1, R2: mirror radii of curvature
    L: cavity length
    m: mode index for frequency
    n: refractive index (default 1)
    """


    # Calculate wavelength in the medium
    modes = Modes(L, m, n=n)  # Make sure this returns frequency
    Lam = c / modes  # Wavelength in vacuum
    Lam_eff = Lam / n  # Wavelength inside medium if refractive index != 1

    if R1 < R2:
        z2 = Position_z2(R1, R2, L)
        z2 = -z2
        z1 = Position_z1(z2, L)  # z1 = z2 - L

        zr = Rayleigh_range(z1, R1)  # z_R
        w0 = Waist_of_the_beam(Lam_eff, zr)  # w0
    else:
        z1 = Position_z2(R2, R1, L)
        z1 = -z1 # making sure the translation works right
        z2 = Position_z1(z1, L)  # z2 = z1 + L

        zr = Rayleigh_range(z2, R2)
        w0 = Waist_of_the_beam(Lam_eff, zr)
        if Ray == True:
            return w0 ,z1, z2, zr
    return w0, z1, z2
def Calculate_the_the_order(distance, wavelegnth):
    return  (distance/wavelegnth - np.pi/16)/ (np.pi/2)
def Calculate_the_wavelength(distance , order):
    return distance/(np.pi/16 + np.pi/2*order)


def Function_modulation(z , wave_vector , distance , n  ):
    return np.sin( 2*wave_vector *  z * n ) * np.cos( 2* wave_vector *z* n  +  2* wave_vector * distance )



def Select_replace_nans(arr):
    for i, element in enumerate(arr):
        if np.isnan(element):
           arr[i]= 0
    return arr


def GRIN_lenses_waist(g, L , lam, waist ,z  ):
    return  np.array( waist / ( np.sqrt( (np.cos(g*L) - g*z *np.sin(g*L)  )**2  +  (g*Raylaigh_range_1(waist, lam ) *np.sin(g*L) ) **2)) )

def Pitch_to_length(g ):
    return  np.pi *2 * 0.25 / g