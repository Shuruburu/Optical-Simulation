#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 16:51:37 2026

@author: shurubura
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
c = 3*1e8

import numpy as np

import Algorithms as al
import scripts
import tracks
import logging
import Gaussian_beam as Gaus  
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
if __name__ == "__main__":
    
        
    FSR = tracks.FSR(wavelength= 1550e-9)
    logger.info("Our FSR is equal to %e" , FSR)
    

    #Gaus.Visualize_the_beam()
    #scripts.Parameter_sweep(tracks.Track)
    scripts.Cavity_waist(tracks.Track , "inf" , Show = True)
    scripts.Cavity_waist(tracks.Track_concave_concave , "Concave_15[mm]-concave_15[mm]" , Show = True, parameters= {"front_mirror": -15e-3 , "back_mirror" :  15e-3 } )
    scripts.Cavity_waist(tracks.Track_concave_concave , "Concave_7[mm]-concave_7[mm]", Show  = True , parameters =  {"front_mirror": -7e-3 , "back_mirror" :  7e-3 } )

    scripts.Coupling_run_all(tracks.Track_concave_concave)
  
    