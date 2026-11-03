#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 10:24:16 2026

@author: shurubura
"""

import numpy as np
import Algorithms as al
import Gaussian_beam as Gaus
import logging 
import tracks
from pathlib import Path
import copy
import json
import xarray as xr
import matplotlib.pyplot as plt
from itertools import product
from multiprocessing import Pool
c = 3*1e8

hbar = 1.054*1e-34
n_quartz = 1.55 

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)



"""
 the funny n letter tellus us , but we have internal losses 
 
"""

from matplotlib.lines import Line2D  


def Different_mirror_type(L , lam ):
   
    #Ls = np.linspace(9e-3, 24e-3,  10)
    L = np.array([11.84e-3])
    print("Concave concave mirror")
    R1 = np.linspace(L[0] +0.1e-3, 1e-1, 200)
    R2 = np.linspace(L[0] + 0.1e-3, 1e-1, 200)
    al.Simulate_different_R1_R2(L  , lam , R1 , R2 )
    print("Planar concave mirror")
    R1 = np.array([np.inf])
    R2 = np.linspace(L[0]+ 0.1e-3, 1e-1, 200)
    
    arr_waists , arr_z2  =al.Simulate_different_R1_R2(L , lam, R1 , R2 , Plot  = False)
    #print(arr_waists.shape) 
    al.Plot(R2*1e3 , arr_waists, 'R1 [millimeters]' ,  'waist [micrometers]' )
    z2   = -3e-3
   # Fixing_the_position(L  , lam  , z2, R2)
    Cavity_Mode = ((L[0] * 2) / lam)
    waist, z1,  z2 , zr =al.Cavity(np.inf, 15e-3, L[0], Cavity_Mode, n=1 , Ray = True)

    print("Our waist with te change length should be roughly equal to {}".format(waist))

def Optmiaztion_waist(L, Lam):
    width: float = 1e-3 # width of the HBAR
    z0_arr= np.array([3.5e-3]) # postion of the sample  : In the old set up we want to have the highest average width gie
    
    R1 = np.array([np.inf])
    R2 = np.linspace(L[0]+ 0.1e-3, 1e-1, 200)
    arr_average_waists=al.Position_of_sample_variation(L, Lam , R1 ,R2, z0_arr ,  width)
    min1 = np.inf
    mirror = 0
    for average_waists in arr_average_waists:
        for mirrors in average_waists:
            
            if average_waists[mirrors] < min1:
                min1 = average_waists[mirrors] 
                mirror = mirrors
        print(f"Minmal average waist has been acheived for {mirror}  and it value was {min1}")

def Visualization_of_GRIN(Lam , waist ):
    
    #NEED to be improved
    dict  =  {"g": 357.0, "pitch": 0.25263157894736843, "d1": 0.00044719063545150503, "d": 0.009690987124463517}
    gs =  g
    L = al.Pitch_to_length(gs)
    L  = 5e-3
    z = 0
    waists = al.GRIN_lenses_waist(gs , L , Lam, waist, z)
    
    al.Plot( gs   ,[waists*1e6], "GRIN Index[1/m]"  ," waist [microns]" )
    g  = 2000
    Ls = np.linspace(4.26e-3 , 7e-3 , 400)    
    waists = al.GRIN_lenses_waist(g , Ls , Lam, waist, z)
    
    al.Plot(  Ls*1e3   ,[waists*1e6], "Length[mm]"  ," waist [microns]" )
    waist = np.linspace((10.4/2)*1e-6 ,  40e-6 , 200)
    L  =   al.Pitch_to_length(g)
    waists = al.GRIN_lenses_waist(g , L , Lam, waist, z)

    al.Plot(waist*1e6   ,[waists*1e6], "waist[microns]"  ," waist [microns]" )

def Cavity_waist( Track , save_track , parameters = {}, Prestent_the_track = False, Show = False):
    """
    Parameters: Expected a dictionary with the name of the parameter as inside of the Track library
    """
    My_track = Track(**parameters)
    if Prestent_the_track:
        tracks.Track_presenter(My_track)
    tracks.Track_save(My_track , Path(__file__).parent / "save" / save_track)
    Simulation=Gaus.Beam_evolution(My_track, 1550e-9)
    R1=tracks.Track_front_mirror(My_track)
    R2 = tracks.Track_back_mirror(My_track)
    Gaus.Find_the_waist_analyticaly(My_track)
    logger.debug(" ###################################################################")
    waist_at_the_front_mirror ,waist = Simulation.Find_the_waist()
    waist_at_the_front_mirror  , waist ,=Simulation.Adjust_for_the_radius(waist)
    logger.debug("##################################################################")
    
    if Show == True:
        points = 6000
        positions = np.linspace(0, tracks.Find_length_track(My_track) ,  points)
        wavelength = 1550e-9
        wavelengths = wavelength*np.ones(points)
        beam_lists = []
        direction = "right"
        zr =  Gaus.Raylaigh_range_1(waist_at_the_front_mirror , wavelength)
        beam =Gaus.Gaussian_Beam(zr ,  R1 , Position= False )
        beam_lists.append(copy.deepcopy(beam))
        wavelengths = Gaus.Propagate(beam_lists ,positions, direction , My_track, wavelengths, Refraction_curve= True , )
        Plot = Gaus.Plot_the_beam(beam_lists , positions , wavelengths)
        Plot.Plot_the_track(f"Concave-Concave R1 = {R1*1e3}[mm] and R2 = {R2*1e3}[mm]")
        Plot.Plot_the_track_radius(f"radius of curvture for gaussian beam inside R1 = {R1*1e3}[mm] and R2 = {R2*1e3}[mm]")

def Parameter_sweep(Track):
    frequency = tracks.Calculate_the_briouilin_frequnecy(1.527)
    L = 3 * 1e8 / (2 * frequency)
    logger.info("Targeted brillouin frequency: %e", L)
    
    gap_position = np.linspace(1e-6, 10e-6, 100)
    distance_to_the_front_mirror = np.linspace(1e-3, 6e-3, 1000)
    radius_of_curvature = np.linspace(50e-3, 1.7, 500)
    # --- Sweep 1: gap_position ---
    waists_gap = []
    for gap in gap_position:
        distance = tracks.Track_length_adjuster(Track(gap_length=gap))
        
        distance_to_the_mirror = L - distance
        my_track = Track(gap_length=gap, d4=distance_to_the_mirror)
        waist_at_the_front_mirror ,waist=Gaus.Find_the_waist_analyticaly(my_track)
        waists_gap.append(waist)
        waist0 = waist
    Gaus.Plot_the_sweep(waists_gap,      gap_position * 1e6,                'gap_sweep',      'Gap length (µm)',        'Beam waist (µm)', 'Waist vs Gap Position',  
                        mark_x=  3 , mark_y = waist*1e6)
    # --- Sweep 2: distance_to_the_front_mirror ---
    
    waists_distance = []
    for d_mirror in distance_to_the_front_mirror:
        distance = tracks.Track_length_adjuster(Track(distance_to_the_sample=d_mirror))
        
        distance_to_the_mirror = L - distance
        my_track = Track(distance_to_the_sample =d_mirror , d4=distance_to_the_mirror)
        waist_at_the_front_mirror ,waist=Gaus.Find_the_waist_analyticaly(my_track)
        waists_distance.append(waist)
    Gaus.Plot_the_sweep(waists_distance, distance_to_the_front_mirror * 1e3, 'distance_sweep', 'Distance to sample (mm)', 'Beam waist (µm)', 'Waist vs Distance to Sample' ,
                         mark_x = 4 , mark_y= waist0 *1e6)
    # --- Sweep 3: radius_of_curvature ---
    waists_roc = []
    for roc in radius_of_curvature:
        distance = tracks.Track_length_adjuster(Track(R=roc))
        distance_to_the_mirror = L - distance
        my_track = Track(R=roc, d4=distance_to_the_mirror)
        
        waist_at_the_front_mirror ,waist=Gaus.Find_the_waist_analyticaly(my_track)
        waists_roc.append(waist)
    Gaus.Plot_the_sweep(waists_roc,      radius_of_curvature * 1e3,          'roc_sweep',      'Radius of curvature (mm)','Beam waist (µm)', 'Waist vs Radius of Curvature',
                        mark_x = 1300 , mark_y= waist0*1e6 )
    """
    waists_list = []
    waists = []
    
    for roc in radius_of_curvature:
        waists = []
        for d_mirror in distance_to_the_front_mirror:
            distance = tracks.Track_length_adjuster(Track(distance_to_the_sample=d_mirror , R=roc))
            distance_to_the_mirror = L - distance
            my_track = Track(distance_to_the_sample =d_mirror , d4=distance_to_the_mirror)
            waist=Gaus.Find_the_waist_analyticaly(my_track)
            waists.append(waist)
        waists_list.append(waists)
    Gaus.plot_waist_colormap(waists_list , radius_of_curvature , distance_to_the_front_mirror  , " Color map of the Radius of curvuture and the postion of the sample")
    """   
    return (
        waists_gap,
        waists_distance,
        waists_roc,
    )

def Grin_lens(Cavity_track ,Lens_Track, parameters =  None):
    path = Path(__file__).parent  / "save" / " GRIN lenses paramater_in_the_lab.json"
    Simulation=Gaus.Beam_evolution(Lens_Track, 1550e-9)
    Simulation.GRIN_lens_distance_mirror(path , parameters = parameters )
    
    with open( path , "r") as f:
            solution = json.load(f)
    
    
    Simulation1 = Gaus.Beam_evolution(Cavity_track() , 1550e-9)
    
    Simulation1.Match_lens_and_cavity(solution , Lens_Track )

def Waist_Radius_2d_depadance(Track, start  = 12e-3, stop= 40e-3 ):
    R2s = np.linspace(start ,stop ,200)
    R1s = -R2s
    best_waist_list  = []
    for R1 in R1s:
        list1 = []
        for R2 in R2s:
            front_waist, best_waist = Gaus.Find_the_waist_analyticaly(Track(front_mirror = R1 , back_mirror = R2))
             # here we try to vary the plotting 
            list1.append(best_waist*1e6)
        best_waist_list.append(list1)
    
    Gaus.plot_waist_colormap(best_waist_list ,  -R1s*1e3 , R2s*1e3 , "Radius of curvuture vs and waists[microns] confocal.png", 
                           xlabel= "R1 (mm)" , ylabel= "R2 (mm) ")
    
def Higher_order_modes(Track):
    wavelength = 1550e-9
    My_track = Track()
    R1_vals = np.linspace(-6e-3, -10e-3, 100)
    R2_vals = np.linspace(6e-3, 10e-3, 100)
    higer_order_modes = np.arange(1, 3)
    Simulation = Gaus.Beam_evolution(My_track, wavelength)
    waist_at_the_front_mirror, waist = Gaus.Find_the_waist_analyticaly(My_track)
    length_of_the_cavity = tracks.Find__optical_length_track(My_track)

    cm = 1 / 2.54  # cm to inch conversion

    difference = xr.DataArray(
        data=np.full((len(R1_vals), len(R2_vals)), np.nan),
        dims=["R1[m]", "R2[m]"],
        coords={
            "R1[m]": R1_vals,
            "R2[m]": R2_vals,
            "Title": "The frequency difference between the fundumental mode and the first order higher mode"
        },
        attrs={
            "wavelength": wavelength,
            "units": "Δv[GHz]",
            "description": "The difference between the 1 order mode and higer order modes_zoomed",
        }
    )
    for m in higer_order_modes:
        for r1 in R1_vals:
            for r2 in R2_vals:
                My_track = Track(front_mirror=r1, back_mirror=r2)
                waist_at_the_front_mirror, waist = Gaus.Find_the_waist_analyticaly(My_track)
                zr = Gaus.Raylaigh_range_1(waist, wavelength)
                delta_v = Gaus.Calucalate_Higher_differenace(m, length_of_the_cavity, zr) * 1e-9
                difference.loc[r1, r2] = delta_v

        plot_map = difference.assign_coords({
            "R1[m]": difference.coords["R1[m]"].values * (-1e3),
            "R2[m]": difference.coords["R2[m]"].values * 1e3,
        }).rename({"R1[m]": "R1[mm]", "R2[m]": "R2[mm]"})

        fig, ax = plt.subplots(figsize=(8 * cm, 6 * cm))
        im = plot_map.plot(ax=ax, add_colorbar=True, cbar_kwargs={'pad': 0.02})

        # Add the 12.6 GHz contour on top
        contour = ax.contour(
            R2_vals * 1e3,
            -R1_vals * 1e3,
            difference.values,
            levels=[12.6],
            colors='white',
            linewidths=1.0
        )
        ax.clabel(contour, fmt='12.6 GHz', fontsize=6, inline=True)

        ax.set_title("Higher order modes", fontsize=8)
        ax.tick_params(labelsize=6)
        ax.xaxis.label.set_size(7)
        ax.yaxis.label.set_size(7)

        fig.tight_layout()
        plt.savefig(Path(__file__).parent / "plots" / f"higher order mode {m}.png", dpi=150, bbox_inches="tight")
        plt.close()

def Length_of_mirrors(diamater_of_the_opitcs):
    radius_of_the_opitcs = diamater_of_the_opitcs/2 # transf
    x = np.linspace(6.28, 10, 100)
    y1, y2 = al.Calculate_the_length_of_the_mirror(x, radius_of_the_opitcs)
    a =al.Calculate_the_length_of_the_mirror(7 , radius_of_the_opitcs)
    print(a)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(x, y1, label='Mirror 1', color='steelblue', linewidth=1.5)
    #ax.plot(x, y2, label='Mirror 2', color='coral', linewidth=1.5, linestyle='--')
    ax.set_xlabel('Radius of curvature (mm)')
    ax.set_ylabel('Spatial extent (mm)')
    ax.set_title(f'Mirror spatial extent vs radius of curvature {diamater_of_the_opitcs/2.54 *1e-1} [in]')
    ax.legend()
    ax.grid(True, linewidth=0.4, alpha=0.5)
    plt.tight_layout()
    plt.savefig(Path(__file__).parent / "plots" /  f" Mirror spatial extend_{diamater_of_the_opitcs/2.54 *1e-1} in.jpg")
    plt.show()

def Coupling_approximation(Track, R1_vals, R2_vals, title_average, title_coupling):

    wavelength = 1550e-9
    points = 500

    # build all (r1, r2) combos
    combos = [
        (r1, r2, Track, wavelength, points)
        for r1, r2 in product(R1_vals, R2_vals)
    ]

    with Pool() as p:
        results = p.map(Gaus.run_single, combos)  # 40,000 runs in parallel

    # fill the xarray from results
    average_waist_map = xr.DataArray(
        data=np.full((len(R1_vals), len(R2_vals), 2), np.nan),
        dims=[
            "Radius of curvature of the back mirror (mm)",
            "Radius of curvature of the front mirror (mm)",
            "Title",
        ],
        coords={
            "Radius of curvature of the back mirror (mm)": R1_vals,
            "Radius of curvature of the front mirror (mm)": R2_vals,
            "Title": ["average_waist_um", "coupling_sum"],
        },
        attrs={
            "wavelength": wavelength,
            "units": "micrometers",
            "description": "Average beam waist and coupling sum inside HBAR crystal",
        }
    )

    for r1, r2, av, coupling_sum in results:
        average_waist_map.loc[r1, r2, "average_waist_um"] = av
        average_waist_map.loc[r1, r2, "coupling_sum"] = coupling_sum

    # Plotting of the arrays
    cm = 1 / 2.54  # cm to inch conversion

    fig, ax = plt.subplots(figsize=(10 * cm, 6 * cm))
    plot_map = average_waist_map.assign_coords({
        "Radius of curvature of the back mirror (mm)":
            average_waist_map.coords["Radius of curvature of the back mirror (mm)"].values * 1e3,
        "Radius of curvature of the front mirror (mm)":
            average_waist_map.coords["Radius of curvature of the front mirror (mm)"].values * 1e3,
    })

    # Plot each metric separately
    waist_plot = (plot_map.sel(Title="average_waist_um") * 1e6)
    waist_plot.attrs["long_name"] = "Average waist"
    waist_plot.attrs["units"] = "µm"
    waist_plot.plot(ax=ax)
    ax.set_title("Average waist inside of the Hbar")
    plt.savefig(Path(__file__).parent / "plots" /  title_average, dpi=300, bbox_inches="tight")
    plt.close()

    fig, ax = plt.subplots(figsize=(10 * cm, 6 * cm))
    plot_map = average_waist_map.assign_coords({
        "Radius of curvature of the back mirror (mm)":
            -average_waist_map.coords["Radius of curvature of the back mirror (mm)"].values * 1e3,
        "Radius of curvature of the front mirror (mm)":
            average_waist_map.coords["Radius of curvature of the front mirror (mm)"].values * 1e3,
    })
    coupling_plot = plot_map.sel(Title="coupling_sum")
    coupling_plot.attrs["long_name"] = "Coupling strength"
    coupling_plot.attrs["units"] = "Hz"
    coupling_plot.plot(ax=ax)
    ax.set_title("The optomechanical coupling to the HBAR")
    plt.savefig(Path(__file__).parent / "plots" / title_coupling, dpi=300, bbox_inches="tight")
    plt.close()

def Senstivity_check():
    base_dir = Path(__file__).resolve().parent
    save_dir  = base_dir / "plots"
    wavelength = 1550e-9

    distance = 4e-3
    n = 1.54
    order = al.Calculate_the_the_order(distance, wavelength)

    z = np.linspace(0, 1e-3, 100)

    distance_differance = np.linspace(1, 50, 100) * 1e-6

    order_diff = al.Calculate_the_the_order(distance, wavelength)

    delta_wavelength = al.Calculate_the_wavelength(distance_differance, order_diff)

    wavelength = np.linspace(1549e-9, 1551e-9, 1000)
    distance_variation = np.linspace(-10e-6, 10e-6, 100) + distance 
    K_vector = 2 * np.pi / wavelength

    integral_modulation = np.zeros((distance_variation.shape[0], wavelength.shape[0]))

    for j, dist in enumerate(distance_variation):
        for i, k in enumerate(K_vector):
            modulation_values = al.Function_modulation(z, k, dist , n )
            integral_modulation[j][i] = np.trapezoid(modulation_values, z)

    # Normalize so the data spans -1 to 1 (g / g_max)
    g_max = np.max(np.abs(integral_modulation))
    integral_modulation = integral_modulation / g_max

    # cm -> inches conversion (matplotlib figsize is in inches)
    cm = 1 / 2.54
    fig_width = 14 *cm 
    fig_height = 10 *cm

    # 1D plot: modulation vs wavelength, at the first distance_variation value
    plt.figure(figsize=(fig_width, fig_height))
    plt.plot(wavelength * 1e9, integral_modulation[0, :])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('g/g_max')
    plt.ylim(-1, 1)
    plt.title('Modulation vs Wavelength')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig( save_dir   /  "the phazematching modulation vs the wavelength.png")
    plt.show()

    # 2D plot: modulation vs wavelength AND distance_variation
    plt.figure(figsize=(fig_width, fig_height))
    plt.pcolormesh(wavelength * 1e9, distance_variation * 1e6, integral_modulation, shading='auto', vmin=-1, vmax=1)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Distance  (µm)')
    plt.title('Phaze matching')
    plt.colorbar(label='g/g_max')
    plt.tight_layout()
    plt.savefig( save_dir / "Distance_phaze.png", dpi = 300 )
    plt.show()
    return wavelength, order, delta_wavelength, integral_modulation 

def Coupling_run_all(Track):
    R1_vals = -np.linspace(12e-3, 35e-3, 300)
    R2_vals = np.linspace(12e-3, 35e-3, 300)   
    Coupling_approximation(Track, R1_vals , R2_vals , "Average_waist_confocal_old.png", "Coupling_confcocal_old.png")
    #R1_vals = -np.linspace(6e-3, 9e-3, 300)
    #R2_vals = np.linspace(6e-3, 9e-3, 300)   
    #Coupling_approximation(Track, R1_vals , R2_vals , "Average_waist_cocentric_old.png","Coupling_cocentric_not_adjusted.png")
    R1_vals = -np.linspace(12e-3, 20e-3, 300)
    R2_vals = np.linspace(12e-3, 20e-3, 300) 
    Coupling_approximation(Track, R1_vals , R2_vals , "Average_waist_confocal_zoomed_old.png","Coupling_cofocal_zoomned_not_adjusted.png")