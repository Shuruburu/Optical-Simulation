
import numpy as np
import scipy as sci
import copy
import logging
import matplotlib.pyplot as plt
from pathlib import Path
import cmath
import tracks
import xarray as xr
import json
from matplotlib.colors import LogNorm
import math
from lmfit import Model
import os
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
# Need adjust for the track !
CM_TO_INCH = 1 / 2.54

h_bar =  1.054* 1e-34 # reduced plancs constant
speed_of_light = 3 * 1e8 # we are pretty fast
def Find_the_position_to_the_medium(Track, sample = "Quartz"):
    prev = None
    Flag = False
    for key in Track.keys():
        
        if key == sample:
            if prev == None:
                logger.error("U can't have a cavity without one of the mirror please make a correct track")    
            if Track[prev]["n"] == 1:
                return abs(Track[prev]["position_beg"] - Track[prev]["position_end"])
            else: 
                Flag  = True # Flaging that the program should end with the next iteration (this means that the gap between the medium and the mirror is )
        elif Flag == True:
            return abs(Track[key]["position_beg"] - Track[key]["position_end"])
        prev = key
    
def Concave_concave(Track, wavelength=1550e-9, n = 4000):
    L = tracks.Find_length_track(Track)
    R2 = tracks.Track_back_mirror(Track)
    R1s = np.linspace(R2, 50e-3, n)
    target_position = Find_the_position_to_the_medium(Track)

    waists = []
    tol = 0.01
    best_error = np.inf
    best_radius = None  # Fix 3: initialize before loop
    best_z1 = None
    best_w0 = None

    for R1 in R1s:
        my_cav = Cavity(L, R1, R2)
        w0 = my_cav.Waist(wavelength)
        waists.append(w0)
        z1 = my_cav.z1
        error = abs(target_position - abs(z1))

        if error < best_error:  
            best_error = error
            best_radius = R1   
            best_z1 = z1
            best_w0 = w0

    if best_error < target_position * tol:
        logger.info(
            "The correct position of the waist has been found for the %e radius of curvature and it was %e with waist %e",
            best_radius, best_z1, best_w0  #
        )
        return best_radius  
    else:
        logger.warning(
            "No radius found within tolerance. Best error: %e (tolerance: %e)",
            best_error, target_position * tol
        )
        return None  

    

    #Adjust the front mirror and the proced with finding the object 
    #then change the procceed with finding waist inside of the cavity that stablizes the the waist using the back propagation


def Raylaigh_range_1(waist, Lam):
    return  waist**2 *  np.pi / Lam


class Analytical:
    def __init__(self, Track ,R1 , R2, wavelength):
        self.Track = Track
        self.matrix = np.identity(2)
        self.solution = None
        self.Solution_1 = None
        self.R1 = R1
        self.R2 = R2
        self.wavelength = wavelength
    # ---- Individual matrices ----
    def propagation(self, d,  n  ):
        return np.array([[1, d/ n],
                         [0, 1]])

    def grin(self, g, L):
        return np.array([[np.cos(g * L), (1 / g) * np.sin(g * L)],
                         [-g * np.sin(g * L), np.cos(g * L)]])

    def lens(self, f):
        return np.array([[1, 0],
                         [-1 / f, 1]])

    def curved_reflection(self, R):
        return np.array([[1, 0],
                         [-2 / R, 1]])

    def curved_refraction(self, R, n1, n2):
        return np.array([[1, 0],
                         [-(n2 - n1) / (n2 * R), n1 / n2]])

    def flat_refraction(self, n1, n2):
        return np.array([[1, 0],
                         [0, n1 / n2]])
    def Round_trip(self, n_round_trips=20): # problem here is we don;t start from the and we need to know how to start form the specifc point # which complicates things
        keys = list(self.Track.keys())

        for i in range(n_round_trips):

            # Start at mirror 1 (flat)
           
            self.matrix = self.curved_reflection(-self.R1) @ self.matrix

            # -------- Forward pass (left → right) --------
            prev = None
            for key in keys:
                elem = self.Track[key]

                if prev is not None:
                    prev_elem = self.Track[prev]
                    R_interface = elem["surface_left"]
                    if np.isinf(R_interface):
                        M = self.flat_refraction(prev_elem['n'], elem['n'])
                    else:
                        M = self.curved_refraction(R_interface, prev_elem['n'], elem['n'])
                    self.matrix = M @ self.matrix

                distance = abs(elem["position_end"] - elem["position_beg"])
                self.matrix = self.propagation(distance, elem["n"]) @ self.matrix

                prev = key

            # Hit mirror 2 (curved)
            self.matrix = self.curved_reflection(self.R2) @ self.matrix

            # -------- Backward pass (right → left) --------
            prev = None
            for key in reversed(keys):
                elem = self.Track[key]

                if prev is not None:
                    prev_elem = self.Track[prev]
                    R_interface = elem["surface_right"]
                    if np.isinf(R_interface):
                        M = self.flat_refraction(prev_elem['n'], elem['n'])
                    else:
                        M = self.curved_refraction(R_interface, prev_elem['n'], elem['n'])
                    self.matrix = M @ self.matrix

                distance = abs(elem["position_end"] - elem["position_beg"])
                self.matrix = self.propagation(distance, elem["n"]) @ self.matrix

                prev = key
    def Solve(self):
        #equation cq^2 + (D-A)q - B =  0
        self.Round_trip()
        delta = cmath.sqrt( (self.matrix[1][1]- self.matrix[0][0])**2  + 4 * self.matrix[0][1]*self.matrix[1][0] )
        self.solution = (self.matrix[1][1]- self.matrix[0][0] - delta)/ (2*self.matrix[1][0])
        self.solution_1 = (self.matrix[1][1]- self.matrix[0][0] + delta)/ (2*self.matrix[1][0])

    def Get_the_waist(self):
        return np.sqrt((self.solution.imag* self.wavelength)/np.pi) ,  np.sqrt((self.solution_1.imag* self.wavelength)/np.pi)
    def Get_the_radius(self):
        q_inv = 1 / self.solution
        q_inv1 = 1 / self.solution_1

        real0 = np.real(q_inv)
        real1 = np.real(q_inv1)

        R0 = np.inf if np.isclose(real0, 0) else 1 / real0
        R1 = np.inf if np.isclose(real1, 0) else 1 / real1

        return R0, R1
    def Calculate_w(self , wavelength):
        
        q_inv = 1 / self.solution
        q_inv1 = 1 / self.solution_1
        imagine1 = np.imag(q_inv)
        imagine2 = np.imag(q_inv1)
        
        w_z1  = np.sqrt(wavelength/(imagine1* np.pi))
        
        w_z2 = np.sqrt(wavelength/(imagine2* np.pi))
        if np.isnan(w_z2):
            return w_z1 
        elif np.isnan(w_z1):
            return w_z2
        else:
            logger.warning("Something is weird inside of wZ calcualtion ")
            return None
        
class Plot_the_beam:
    def __init__(self, beams, z, wavelength):
        self.beams = beams
        self.z = np.asarray(z)
        self.wavelength = wavelength  # same shape as z
        self.relative_path = Path(__file__).parent / "plots"
        self.fig_size = (14 * CM_TO_INCH, 10 * CM_TO_INCH)  # width=14cm, height=10cm
        self.color_map = {
            0: "red", 1: "blue", 2: "green", 3: "yellow", 4: "purple",
            5: "orange", 6: "cyan", 7: "magenta", 8: "brown", 9: "pink"
        }

    def waist_from_q(self, q_array, wavelength_array):
        inv_q = 1 / q_array
        return np.sqrt(-wavelength_array / (np.pi * np.imag(inv_q)))

    def _highlight_mirror_gap(self, mirror_gap, color="green", alpha=0.15, label="GRIN Mirror distance (d)"):
        """mirror_gap = (L, d): shade the two 'd'-length free-space regions
        adjacent to the mirror, assuming symmetric layout centered at L/2."""
        if mirror_gap is None:
            return
        begin, d = mirror_gap
        begin_mm = begin * 1e3
        d_mm = d * 1e3
        plt.axvspan(begin_mm, begin_mm + d_mm, color=color, alpha=alpha, label=label)

    def Plot_the_track(self, name, mirror_gap=None, save = None):
        q_array = np.array([beam.q for beam in self.beams])
        w_values = self.waist_from_q(q_array, self.wavelength) * 1e6  # µm
        z_mm = self.z * 1e3  # mm

        plt.figure(figsize=self.fig_size)
        plt.plot(z_mm, w_values, 'b', label="Beam envelope")
        plt.plot(z_mm, -w_values, 'b')
        plt.fill_between(z_mm, -w_values, w_values, alpha=0.3)
        plt.axhline(0, linestyle='--', linewidth=1)

        self._highlight_mirror_gap(mirror_gap)

        plt.xlabel("Propagation distance z [mm]")
        plt.ylabel("Beam radius [µm]")
        plt.title(name)
        plt.grid(True)
        plt.legend(loc='upper right')
        if save == None:
            plt.savefig(self.relative_path / f"{name}.png" , dpi = 300)
        else:
            plt.savefig(self.relative_path / f"{save}.png" , dpi= 300)
        plt.show()

    def Plot_the_track_radius(self, name, mirror_gap=None , save = None):
        q_array = np.array([beam.q for beam in self.beams])
        q_inv = 1 / q_array
        a = np.real(q_inv)
        R_values = 1 / a * 1e3  # mm
        z_mm = self.z * 1e3  # mm

        plt.figure(figsize=self.fig_size)
        plt.plot(z_mm, R_values, 'b', label="Beam envelope")
        plt.axhline(0, linestyle='--', linewidth=1)

        self._highlight_mirror_gap(mirror_gap)

        plt.xlabel("Propagation distance z [mm]")
        plt.ylabel("Beam radius of curvature [mm]")
        plt.ylim(-100, 100)
        plt.title(name)
        plt.grid(True)
        plt.legend(loc='upper right')
        if save == None:
            plt.savefig(self.relative_path / f"{name}.png" , dpi = 300)
        else:
            plt.savefig(self.relative_path / f"{save}.png" , dpi =300)
        plt.show()

    def Plot_the_relation(self, name, ds):
        w_values_list = []
        for beams in self.beams:
            w_values = []
            for beam in beams:
                w = self.waist_from_q(beam.q) * 1e6  # µm
                w_values.append(w)
            w_values_list.append(np.array(w_values))

        z_m = self.z
        plt.figure(figsize=self.fig_size)
        for nr, w_values in enumerate(w_values_list):
            plt.plot(z_m, w_values, color=self.color_map[nr], label=f"distance to lens: {ds[nr]*1e3:.2f} mm")

        plt.grid(True)
        plt.xlabel("GRIN index [1/m]")
        plt.ylabel("Beam radius after GRIN lens [µm]")
        plt.title("Gaussian Beam Propagation (Envelope)")
        plt.legend(loc='upper right')
        plt.savefig(self.relative_path / f"{name}.png", dpi = 300)

    def Show_the_plot(self):
        plt.show()




#Class simulating the profile of the Gaussian Beam propagatin inside of the medium
class Gaussian_Beam:
    """
    Gaussian beam represented using the complex q-parameter:
        q = z + i * z_R

    The beam is transformed using ABCD matrices:
        q' = (Aq + B) / (Cq + D)
    """

    def __init__(self, Rayleigh_range: float, z: float,  Position = True):
        """
        Parameters:
        -----------
        Rayleigh_range : float
            Rayleigh range (z_R) [meters]
        z : float
            Beam position [meters]
        """
        if Position == True :
            self.q = complex(z, Rayleigh_range)
        else:
            a =1/z
            b =  - 1 /Rayleigh_range
            q_inv = complex(a ,  b  )
            self.q = 1/q_inv
            self.zr = None
    def _apply_matrix(self, M: np.ndarray):
        """
        Apply an ABCD matrix to the Gaussian beam.

        Parameters:
        -----------
        M : numpy.ndarray
            2x2 ABCD matrix [[A, B], [C, D]]
        """
        A, B = M[0]
        C, D = M[1]

        # ABCD law for Gaussian beams
        self.q = (A * self.q + B) / (C * self.q + D)

    def propagation(self, d: float , n):
        """Free-space propagation"""
        M = np.array([[1, d/n],
                      [0, 1]])
        self._apply_matrix(M)

    def grin(self, g: float, L: float):
        """GRIN medium"""
        M = np.array([[np.cos(g * L), (1 / g) * np.sin(g * L)],
                      [-g * np.sin(g * L), np.cos(g * L)]])
        self._apply_matrix(M)

    def lens(self, f: float):
        """Thin lens"""
        M = np.array([[1, 0],
                      [-1 / f, 1]])
        self._apply_matrix(M)

    def curved_reflection(self, R: float):
        """Curved mirror reflection"""
        M = np.array([[1, 0],
                      [- 2 / R, 1]])
        self._apply_matrix(M)

    def curved_refraction(self, R: float, n1: float, n2: float):
        """Curved interface refraction"""
        M = np.array([[1, 0],
                      [-(n2 - n1) / (n2 * R), n1 / n2]])
        self._apply_matrix(M)

    def flat_refraction(self, n1: float, n2: float):
        """Flat interface refraction"""
        M = np.array([[1, 0],
                      [0, n1 / n2]])
        self._apply_matrix(M)

    def get_q(self):
        """Return current q-parameter"""
        return self.q
    def radius_of_curvature(self):
        """
        Compute the radius of curvature R from the q-parameter.

        Returns:
        --------
     R : float
          Radius of curvature [meters]
        """
        inv_q = 1 / self.q
        if inv_q.real == 0:
            return np.inf
        return 1 / inv_q.real
    def colimate(self, g  ):
        M = np.array([[0, (1 / g)],
                      [-g , 0]])
        self._apply_matrix(M)
    def beam_waist(self, wavelength: float):
        """
        Compute the beam waist w from the q-parameter.

        Parameters:
        -----------
        wavelength : float
         Wavelength of the beam [meters]

        Returns:
         --------
        w : float
            Beam waist [meters]
        """
        inv_q = 1 / self.q
        w = np.sqrt(-wavelength / (np.pi * inv_q.imag))
        return w
    def Calucalate_w0(self , wavelength):


        q = self.q  # NOT 1/q

        z_R = np.imag(q)

        # waist from Rayleigh range
        w0 = np.sqrt(wavelength * z_R / np.pi)

        return w0

def Calucalate_Higher_differenace(  m ,  l,  z0, c = 3*1e8 , ):
    return    m*c/(2*l) *2* np.arctan(l/(2*z0) /np.pi)

def Calucalate_Higher(  z0 , long_mode ,  Cavity_length , transveral_mode, c = 3*1e8):
    return  long_mode +  transveral_mode*c/(2*Cavity_length)* np.arctan(Cavity_length/(2*z0) /np.pi)

class Beam_evolution:
    def __init__(self, Track , wavelength):
        self.wavelength = wavelength 
        self.Track = Track

    def Find_the_waist(self, n_round_trips=4, tol=0.5e-3):
        My_track = self.Track
        wavelength = self.wavelength

        L = tracks.Find__optical_length_track(My_track)
        R1 = tracks.Track_front_mirror(My_track)
        R2 = tracks.Track_back_mirror(My_track)

        logger.debug("The length of the cavity is %e", L)
        logger.debug("The curvature of the first mirror is %e", R1)
        logger.debug("The curvature of the second mirror is %e", R2)
        logger.debug("The wavelength of the laser is %e", wavelength)

        My_cav = Cavity(L, -R1, R2)
        My_cav.Waist(wavelength)

        sweep = np.linspace(40e-6 , 300e-6, 1000)
        logger.debug("The cavity waist %e", My_cav.waist)

        best_waist = None
        best_error = np.inf
        best_Radius = None
        best_flat = None

        for waist in sweep:
            zr = Raylaigh_range_1(waist, wavelength)
            if R1 == np.inf:
                beam = Gaussian_Beam(zr, 0)
            else:
                beam = Gaussian_Beam(zr , R1, Position= False)
            q_int = beam.q
           
            for _ in range(n_round_trips):
                self._forward_pass( beam , My_track)
                
                beam.curved_reflection(R2)
                self._backward_pass( beam , My_track)
                
                beam.curved_reflection(-R1 ) #we add minus beacuse of change of the vectors direction
                Radius = beam.radius_of_curvature()
            error = abs(beam.q - q_int)
            other_radius_flat = beam.radius_of_curvature()

            if error < best_error:
                best_error = error
                best_waist = waist
                best_Radius = Radius
                best_flat = other_radius_flat

        if best_error < tol:
            best_waist = best_waist
            dummy = Gaussian_Beam(Raylaigh_range_1(best_waist ,  wavelength ) , R1 ,  Position = False)
            waist0 = dummy.Calucalate_w0(wavelength)
            logger.info("The waist has been found %e [microns](error %e) and the waist at the front mirror is %e [mm]",  waist0*1e6  ,  best_error , best_waist*1e6 )
            
            logger.info("The radius of curvature at curved mirror: %e", best_Radius)
        else:
            logger.warning(
                "Closest waist found %e (error %e), radius %e, %e",
                best_waist, best_error, best_Radius, best_flat
            )

        return best_waist , waist0

    def Adjust_for_the_radius(self, waist, tol=0.01):
        My_track = self.Track
        wavelength = self.wavelength

        L = tracks.Find_length_track(My_track)
        R1 = tracks.Track_front_mirror(My_track)
        R2 = tracks.Track_back_mirror(My_track)

        best_waist = waist
        best_error = np.inf
        best_radius = None

        logger.info("The length of the cavity is %e", L)
        logger.info("The curvature of the first mirror is %e", R1)
        logger.info("The curvature of the second mirror is %e", R2)
        logger.info("The wavelength of the laser is %e", wavelength)
        logger.info("The tolerance of the radius of curvature is %e percent", tol * 100)

        sweep = np.linspace(waist - 20e-6, waist + 40e-6, 5000)

        for waist_sweep in sweep:
            zr = Raylaigh_range_1(waist_sweep, wavelength)
            if R1  == np.inf:
                beam = Gaussian_Beam(zr, 0)
            else:
                beam = Gaussian_Beam(zr , R1 ,  Position= False)
            q_int = beam.q
            self._forward_pass( beam , My_track)
            Radius = beam.radius_of_curvature()
            error = abs(R2 - Radius)

            if best_error > error:
                best_error = error
                best_waist = waist_sweep
                best_radius = Radius

        if best_radius is not None and abs(best_radius - R2) <= R2 * tol:
            best_waist = best_waist
            dummy = Gaussian_Beam(Raylaigh_range_1(best_waist ,  wavelength ) , R1 ,  Position = False)
            waist0 = dummy.Calucalate_w0(wavelength)
            
           
            logger.info(
                "Found the best waist that matches the radius of curvature: radius=%e, waist=%e  , waist at the front mirror=%e" ,
                best_radius * 1e3,
                waist0 *1e6,
                best_waist *1e6,
                
            )
            return best_waist , waist0

        logger.warning(
            "The waist has not converged. Best error=%e, best_radius=%e, best_waist=%e",
            best_error,
            best_radius if best_radius is not None else float("nan"),
            best_waist,
        )
        dummy = Gaussian_Beam(Raylaigh_range_1(best_waist ,  wavelength ) , R1 ,  Position = False)
        waist0 = dummy.Calucalate_w0(wavelength)
        
        
        return best_waist , waist0
    
    def _backward_pass(self, beam , My_track):
        """Propagate beam backward through all track elements, applying refractions and propagation."""
        My_track 
        prev = None

        for key in reversed(My_track.keys()):
            elem = My_track[key]
            distance = abs(elem["position_beg"] - elem["position_end"])

            if elem["type"] == "medium":
                if prev is not None:
                    prev_elem = My_track[prev]
                    if elem["surface_right"] == np.inf:
                        beam.flat_refraction(prev_elem['n'], elem['n'])
                    else:
                        beam.curved_refraction(elem["surface_right"], prev_elem['n'], elem['n'])

                beam.propagation(distance, elem["n"])

            elif elem["type"] == "lens":
                beam.grin(elem['g1'], elem['l1'])

            prev = key

    def _forward_pass(self, beam , My_track):
        """Propagate beam forward through all track elements, applying refractions and propagation."""
        My_track 
        prev = None

        for key in My_track.keys():
            elem = My_track[key]
            distance = abs(elem["position_beg"] - elem["position_end"])

            logger.debug(elem["type"])
            if prev is not None:
                prev_elem = My_track[prev]
                if elem["surface_left"] == np.inf:
                    beam.flat_refraction(prev_elem["n"], elem["n"])

                else:
                    beam.curved_refraction(elem["surface_left"], prev_elem["n"], elem["n"])
            if elem["type"] == "medium": 
                beam.propagation(distance, elem["n"])
            elif elem["type"] == "lens":
                beam.grin(elem['g1'], elem['l1'])

            prev = key

    def GRIN_lens_distance_mirror(self, save, R_val=-7e-3, parameters=None, tol=0.01):
        Fiber_waist = 10.4 / 2 * 1e-6

        if parameters is not None:
            grin_index = parameters["grin_index"]
            length = parameters["length"]
            refractive_index = parameters["refractive_index"]
            pitch = np.linspace(0.15, 0.30, 20)

            grin_indexes = grin_index
            pitches = pitch
            tol = 0.02
            distances_to_the_mirror = np.linspace(2e-3, 50e-3, 700)
            distance_to_the_lens = np.linspace(30e-6, 2e-3, 200)
        else:
            grin_indexes = np.linspace(800, 1200, 100)
            pitches = np.linspace(0.15, 0.30, 20)
            distances_to_the_mirror = np.linspace(1e-3, 3e-3, 100)
            distance_to_the_lens = np.linspace(100e-6, 1e-3, 100)

        solutions = {}
        nr1 = 0
        zr = Raylaigh_range_1(Fiber_waist, self.wavelength)

        for d in distances_to_the_mirror:
            for d1 in distance_to_the_lens:
                for i, g in enumerate(grin_indexes):
                    for p in pitches:
                        Beam = Gaussian_Beam(zr, 0)
                        My_track = self.Track(g, Pitch_to_length(g, p), d1, d, R=R_val, refractive_index=refractive_index[i])
                        self._forward_pass(Beam, My_track)

                    # Match the physical condition used later in Match_lens_and_cavity:
                    # the beam must be evaluated AFTER refracting into the mirror
                    # substrate (index n), not on the air side. Otherwise this
                    # coarse filter tests a different quantity than the fine
                    # q-match and can reject valid lens configurations.
                        back_mirror = tracks.Track_back_mirror(My_track, whole=True)
                        Beam.curved_refraction(-back_mirror["surface_right"], back_mirror["n"], 1)
                        mirror_roc = back_mirror["surface_right"]
                        
                        if abs(Beam.radius_of_curvature() - mirror_roc) < tol * abs(mirror_roc):
                            solutions[nr1] = {"g": g, "pitch": p, "d1": d1, "d": d, "refractive_index": refractive_index[i]}
                            nr1 += 1

        if solutions:
            logger.info("We found at least one combination")
            path = save
            with open(path, "w") as f:
                json.dump(solutions, f)
            return solutions
        else:
            logger.warning("We didn't find anything")
            return None
        
    def Match_lens_and_cavity(self, solutions, lens_track, tol=0.01):
        Cavity_track = self.Track
        R2 = tracks.Track_back_mirror(Cavity_track)
        R1 = tracks.Track_front_mirror(Cavity_track)
        frot_mirror_waist, waist = self.Find_the_waist()
        zr_waist_inside = Raylaigh_range_1(frot_mirror_waist, self.wavelength)
        beam_inside_of_the_cavity = Gaussian_Beam(zr_waist_inside, R1, Position=False)
        self._forward_pass(beam_inside_of_the_cavity, Cavity_track)
        beam_inside_of_the_cavity.curved_reflection(R2)
        print(beam_inside_of_the_cavity.radius_of_curvature())
        possible_set_up = []
        Fiber_waist = 10.4 / 2 * 1e-6
        Flag = False
        zr = Raylaigh_range_1(Fiber_waist, self.wavelength)
        logger.info("The size of the beam waist of at the surface of the mirror is %e ", beam_inside_of_the_cavity.beam_waist(self.wavelength))
        waists = []

        # beam_inside_of_the_cavity.q describes the intracavity mode traveling
        # BACK MIRROR -> FRONT MIRROR (backward), because it was built by a
        # forward pass front->back plus a reflection off R2.
        # Beam.q (built below) describes the injected mode traveling
        # FRONT MIRROR -> BACK MIRROR (forward).
        # These are counter-propagating descriptions of the same plane, so a
        # raw q comparison is meaningless. Convert the cavity beam's q into
        # the forward-propagating convention before comparing:
        #   q_reversed = -conj(q_forward)  <=>  q_forward = -conj(q_reversed)

        for key in solutions.keys():
            elem = solutions[key]
            g = float(elem["g"])
            pitch = float(elem["pitch"])
            distance_to_the_lens = float(elem["d1"])
            distance_to_the_mirror = float(elem["d"])
            refractive_index = float(elem["refractive_index"])
            Lens_Track = lens_track(g, Pitch_to_length(g, pitch), distance_to_the_lens, distance_to_the_mirror, R=R1, refractive_index=refractive_index)
            Beam = Gaussian_Beam(zr, 0)
            self._forward_pass(Beam, Lens_Track)

            print("Beam inside of the cavity {}".format(beam_inside_of_the_cavity.radius_of_curvature()))

            back_mirror = tracks.Track_back_mirror(Lens_Track, whole=True)
            print("After the refraction {}".format(Beam.radius_of_curvature()))
            Beam.curved_refraction(-back_mirror["surface_right"], back_mirror["n"], 1)
            print("After the refraction {}".format(Beam.radius_of_curvature()))
            
            error = abs(beam_inside_of_the_cavity.q - Beam.q)
            print(error)
            if error < tol * abs(beam_inside_of_the_cavity.q):
                waists.append(Beam.beam_waist(self.wavelength))
                possible_set_up.append(elem)
                Flag = True
            #apply the matrix of mirror refarction
            ###################
            # Maybe check how whether the radius of curvture matches the
            ################
        if Flag == True:
            logger.info("The lens has been found")
            path = Path(__file__).parent / "save" / " GRIN set_up.json"
            with open(path, "w") as f:
                json.dump(possible_set_up, f)
            return possible_set_up
        else:
            logger.warning("We couldn't find anything with such this relative tollerance %e ", tol)
            return None
        
    def Optimizaiton_of_the_coupling(self):
        waist_at_the_front_mirror ,waist = self.Find_the_waist()
        Propagate()
        points = 5000
        R1=tracks.Track_front_mirror(self.Track())
        positions = np.linspace(0, tracks.Find_length_track(self.Track()) ,  points)
    
        wavelengths = self.wavelength*np.ones(points)
        beam_lists = []
        direction = "right"
        zr =  Raylaigh_range_1(waist_at_the_front_mirror , self.wavelength)
        beam =Gaussian_Beam(zr ,  R1 , Position= False )
        beam_lists.append(copy.deepcopy(beam))
        wavelengths = Propagate(beam_lists ,positions, direction , self.Track(), wavelengths, Refraction_curve= True , )
        #we do foward pass and  
# SO far so good useless





def Check_the_medium(track, position, prev_position, direction):
    logger.debug(f"We are currently at the {position}")
    if direction == "right":
        
        for element in track.keys():
            
            if (track[element]["position_beg"] <= position and 
                track[element]["position_end"] > position and 
                track[element]["position_beg"] > prev_position):

                logger.debug("we found the refractive position in ")

                for other_elements in track.keys():
                    if (track[other_elements]["position_beg"] < prev_position and 
                        track[other_elements]["position_end"] >= prev_position):

                        logger.debug("we managed to get to the loop")
                        return track[element], track[other_elements]
            
            if (track[element]["position_beg"] < position and 
                track[element]["position_end"] >= position):

                logger.debug(f"We are in the {element} ")
                return track[element], track[element]


    elif direction == "left":

        for element in track.keys():
            emergancy = element
            if (track[element]["position_beg"] < position and 
                track[element]["position_end"] > position and 
                track[element]["position_end"] < prev_position):

                logger.debug("we found the refractive position (left)")

                for other_elements in track.keys():
                    if (track[other_elements]["position_beg"] < prev_position and 
                        track[other_elements]["position_end"] > prev_position):

                        logger.debug("we managed to get to the loop (left)")
                        return track[element], track[other_elements]

            if (track[element]["position_beg"] < position and 
                track[element]["position_end"] > position):

                logger.debug("The other medium (left)")
                return track[element], track[element]
        return  track[emergancy], track[emergancy]
    return np.nan , np.nan
def Propagate(Beam_list, positions, direction, My_Track, wavelengths ,Refraction_curve=True):
    """
    Beam_list: the list of beams with one initial beam
    """
    if direction == "right":
        opposite_direction = "left"
    elif direction == "left":
        opposite_direction = "right"
    else:
        logger.error("direction could only be right or left")
        return

    beam = copy.deepcopy(Beam_list[-1])

    for i in range(1, positions.size):
        z_current = positions[i]
        z_prev = positions[i - 1]
        distances = abs(z_current - z_prev)
        logger.debug(f"Currently in the {i} interatio of the loop")
        current, prev = Check_the_medium(My_Track, z_current, z_prev, direction)
        
        if current["type"] == "medium" and prev["type"] == "medium":
            wavelengths[i] =  wavelengths[i]/current["n"]
            if current != prev:

                if direction == "right":
                    d1 = abs(current["position_beg"] - z_prev)
                    d2 = abs(z_current - current["position_beg"])
                else:
                    d1 = abs(current["position_end"] - z_prev)
                    d2 = abs(current["position_end"] - z_current)

                beam.propagation(d1, prev["n"])

                if current[f"surface_{direction}"] != np.inf and Refraction_curve:
                        if current[f"surface_{direction}"] == prev[f"surface_{opposite_direction}"]:
                            logger.error("The surface can't be concave and flat at the same time")

                        beam.curved_refraction(
                        current[f"surface_{opposite_direction}"],
                        prev["n"],
                        current["n"]
                        )
                else:
                    beam.flat_refraction(prev["n"], current["n"])
            
                beam.propagation(d2, current["n"])

            else:
                        beam.propagation(distances, current["n"])
        elif prev["type"] == "lens" and current["type"] == "lens":
                wavelengths[i] =  wavelengths[i]/current["n"]
                beam.grin(current["g1"] , distances)
        elif current["type"] == "medium" and prev["type"] == "lens":
            wavelengths[i] =  wavelengths[i]/current["n"]
            if direction == "right":
                    d1 = abs(current["position_beg"] - z_prev)
                    d2 = abs(z_current - current["position_beg"])
            else:
                    d1 = abs(current["position_end"] - z_prev)
                    d2 = abs(current["position_end"] - z_current)

            
            beam.grin(prev["g1"] , d1)
            beam.flat_refraction(prev["n"], current["n"])
            beam.propagation(d2, current["n"])

        elif current["type"] == "lens" and prev["type"] == "medium":
            wavelengths[i] =  wavelengths[i]/current["n"]
            if direction == "right":
                    d1 = abs(current["position_beg"] - z_prev)
                    d2 = abs(z_current - current["position_beg"])
            else:
                    d1 = abs(current["position_end"] - z_prev)
                    d2 = abs(current["position_end"] - z_current)

            beam.propagation(d1, prev["n"])
            beam.flat_refraction(prev["n"], current["n"])
            beam.grin(current["g1"] , d2)            
        Beam_list.append(copy.deepcopy(beam))
        
    return wavelengths
# For debugging purpouse
def Beam_list_print(beam_list, wavelength=1550e-9, per_line=100):
    for i, beam in enumerate(beam_list, 1):
        print(beam.beam_waist(wavelength), end=" ")
        if i % per_line == 0:
            print()

def plot_waist_colormap(waists_list, radius_of_curvature, distance_to_the_front_mirror, name, 
                        cmap='viridis', log_scale=False, xlabel="Distance to front mirror (m)", ylabel="Radius of curvature (m)"):
    """
    waists_list : list of lists, shape (len(roc), len(d_mirror))
    """
    Z = np.array(waists_list)   # shape: (n_roc, n_d_mirror)
    X, Y = np.meshgrid(distance_to_the_front_mirror, radius_of_curvature)

    # Convert cm to inches (1 inch = 2.54 cm)
    width_cm, height_cm = 8, 6
    figsize = (width_cm / 2.54, height_cm / 2.54)

    fig, ax = plt.subplots(figsize=figsize)

    norm = LogNorm(vmin=Z.min(), vmax=Z.max()) if log_scale else None
    mesh = ax.pcolormesh(X, Y, Z, cmap=cmap, norm=norm, shading='auto')

    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label("Beam waist (microns)")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Beam waist map")

    fig.tight_layout()
    plt.savefig(Path(__file__).parent / "plots" / name, dpi= 300 )
    plt.show()
    return fig, ax



class Cavity:

    def __init__(self, length, R1, R2 ,  c= 3e8):
        self.Length = length
        self.R1 = R1
        self.R2 = R2
        self.c = c
        self.waist = None
        self.z1 = None
        self.z2 = None
        self.Lam = None
        self.zr = None

    def FSR(self):
        return self.c / (2 * self.Length)

    def FSR_adjust(self, FSR):
        self.Length = self.c / (2 * FSR)

    def Position_z2(self, R1, R2):
        L = self.Length
        return L * (R1 - L) / (R2 + R1 - 2 * L)

    def Position_z1(self, z2):
        return z2 + self.Length

    def Rayleigh_range(self, z, R):
        return np.sqrt(R / z - 1) * z

    def Modes(self, Lam):
        return int((2 * self.Length) / Lam)

    def Frequency_mode(self, FSR, n):
        return FSR * n

    def Waist_of_the_beam(self, Lam_eff, zr):
        return np.sqrt(Lam_eff * zr / np.pi)

    def Waist(self, Lam):
        n = 1

        modes = self.Modes(Lam)
        freq = self.Frequency_mode(self.FSR(), modes)

        Lam = self.c / freq
        Lam_eff = Lam / n

        R1 = self.R1
        R2 = self.R2

        if R1 < R2:
            z2 = -self.Position_z2(R1, R2)
            z1 = self.Position_z1(z2)

            zr = self.Rayleigh_range(z1, R1)
            w0 = self.Waist_of_the_beam(Lam_eff, zr)

        else:
            z1 = -self.Position_z2(R2, R1)
            z2 = self.Position_z1(z1)

            zr = self.Rayleigh_range(z2, R2)
            w0 = self.Waist_of_the_beam(Lam_eff, zr)

        self.waist = w0
        self.Lam = Lam
        self.z1 = z1
        self.z2 = z2
        self.zr = zr

        return w0

def sine_func(x, amplitude, frequency, phase, offset):
    return amplitude * np.sin(2 * np.pi * frequency * x + phase) + offset


    

def Find_the_waist_analyticaly(My_track):
    
    L = tracks.Find_length_track(My_track)
    R1 = tracks.Track_front_mirror(My_track)
    R2 = tracks.Track_back_mirror(My_track)
    wavelength = 1550e-9
    tol = 0.5e-6
    Track_analatyical =Analytical(My_track , R1 ,R2, wavelength)
    Track_analatyical.Solve()
    waist1, waist2=Track_analatyical.Get_the_waist()
    Radius1 , Radius2  = Track_analatyical.Get_the_radius()
    logger.info("The analytical soluton of the first beam was %e , and the second %e" , waist1 , waist2)
    logger.debug("The analytical solution of radius of curvture %e,  and the second %e" , Radius1, Radius2)


    if np.isnan(waist1):
        
        front_mirror_waist = Track_analatyical.Calculate_w(wavelength )
        return front_mirror_waist ,waist2
    else:
        front_mirror_waist=Track_analatyical.Calculate_w(wavelength  )
        return  front_mirror_waist ,waist1

def Pitch_to_length(g, pitch ):
    return  np.pi *2 * pitch / g


def Sensitivity_test(target_beam):
    initial_waist = 10.4 / 2 * 1e-6  #distance dosent change_much 
    wavelength = 1550e-9
    zr = Raylaigh_range_1(initial_waist, wavelength)
    g1 = np.linspace(500 , 1500 , 100 ) # varing the first parameter
    g2 = np.linspace(500 , 1500, 100)
    pitches = np.linspace(0.05, 0.25 ,100)
    pitches1 = np.linspace(0.05, 0.25, 100)
    distance_to_the_lens = np.linspace( 0.)
    d_waist=  zr
    ds = np.linspace( d_waist  , d_waist  + 0.5*1e-3, 8)
    My_track = tracks.GRIN_Track(g , 2220 , Pitch_to_length(g ,0.23) , 420e-69 , d)
    






def Visualize_the_beam(save_dir="beam_plots"):
    os.makedirs(save_dir, exist_ok=True)

    initial_waist = 10.4 / 2 * 1e-6
    wavelength = 1550e-9
    zr = Raylaigh_range_1(initial_waist, wavelength)
    print(zr)

    param_sets = [
        {"g": 357.0, "pitch": 0.24473684210526314, "d1": 0.0002279899497487437, "d": 0.007699570815450644, "refractive_index": 1.616},
        {"g": 138.0, "pitch": 0.20526315789473684, "d1": 0.001891105527638191, "d": 0.04100429184549356, "refractive_index": 1.515},
        {"g": 357.0, "pitch": 0.2368421052631579, "d1": 0.0004457788944723618, "d": 0.0027553648068669528, "refractive_index": 1.616},
        {"g": 138.0, "pitch": 0.25263157894736843, "d1": 0.0007724623115577889, "d": 0.0195793991416309, "refractive_index": 1.515}
    ]

    points = 100

    for params in param_sets:
        g = params["g"]
        pitch = params["pitch"]
        d1 = params["d1"]
        d = params["d"]
        ref = params["refractive_index"]
        beam_lists = []
        beam = Gaussian_Beam(zr, 0)
        beam_lists.append(copy.deepcopy(beam))

        My_track = tracks.GRIN_Track_Mirror(g, Pitch_to_length(g, pitch), d1, d, R=-15e-3, refractive_index=ref)
        L = tracks.Find_length_track(My_track)

        positions = np.linspace(0, L, points)
        wavelengths = wavelength * np.ones(points)
        direction = "right"

        wavelengths = Propagate(
            beam_lists, positions, direction, My_track, wavelengths,
            Refraction_curve=True,
        )

        Plot = Plot_the_beam(beam_lists, positions, wavelengths)

        waist_title = f"GRIN lens coupling optimization[waist]"
        roc_title = f"GRIN lens coupling optimization[ROC]"
        Gap_starting_point = Pitch_to_length(g, pitch) + d1
        
        Plot.Plot_the_track(waist_title, mirror_gap=(Gap_starting_point, d) , save = f"g{g}_pitch{round(pitch,3)}_waist.png")
        plt.close()
        Plot.Plot_the_track_radius(roc_title, mirror_gap=(Gap_starting_point, d), save =  f"g{g}_pitch{round(pitch,3)}_ROC.png")
        plt.close()


def Plot_the_sweep(waists, dependencies, name, xlabel, ylabel, title, mark_x=None, mark_y=None):
    fig, ax = plt.subplots(figsize=(7, 5))

    ax.plot(dependencies, np.array(waists) * 1e6, linewidth=1.8)

    if mark_x is not None and mark_y is not None:
        ax.plot(mark_x, mark_y, marker='*', markersize=18, color='red',
                 markeredgecolor='black', markeredgewidth=0.8,
                 linestyle='None', zorder=5)

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.tick_params(labelsize=10)

    plt.tight_layout()
    plt.savefig(Path(__file__).parent / "plots" / f'{name}_.png', dpi=300, bbox_inches='tight')
    plt.show()

def transveral_part(optical_waist , acoustic_wave):
    return np.sqrt((2 / np.pi  )) * (acoustic_wave ) / (optical_waist **2 /2 +   acoustic_wave**2  )

def Coupling_integral(Track,beam_waists_list , Hbar , average, name , wavelength = 1550e-9):
    n=len(beam_waists_list)
    constants =tracks.Material_constants(name)
    angular_acoustic_frequency = 2* np.pi * tracks.Calculate_the_briouilin_frequnecy(constants["refractive_index"] ) 
    Length_of_hbar = Hbar["position_end"] - Hbar["position_beg"]
    Optical_length  = tracks.Find__optical_length_track(Track)
    frequency  = 2*np.pi *speed_of_light/wavelength
    permitivity  = constants["permittivity"]
    c_33 = constants["stiffness_tensor"]["3_3"]*1e9
    p_13 = constants["piezo_tensor"]["1_3"]
    distance_to_Hbar =  tracks.Distance_to_Hbar(Track)
    #acoustic vector of the mode inside of the HBAR
    q = tracks.Acoustic_wave_vector()
    intervalls  = np.linspace(0,  Hbar["position_end"] - Hbar["position_beg"] , n + 1 )
    acoustic_mode  = 1/ np.sqrt(2) * average 

    acoustic_wav_length = (2*np.pi)/q
    k =  int(distance_to_Hbar/acoustic_wav_length) 
    sum = -20
    correction =  np.linspace(0, 3, 30)*1e-6
    for j, shift in enumerate(correction):
        sum_temp = 0
        for i in range(n):
            f = lambda x: (1/2)*(-np.cos(2*q*x *Hbar["n"]+  q*(distance_to_Hbar  ) )/(4*q*Hbar["n"]) - (1/2)*x *np.sin( q*(distance_to_Hbar ) ))
            Integral = f(intervalls[i+1])  - f(intervalls[i])
            sum_temp+= transveral_part(beam_waists_list[i], acoustic_mode) * Integral
        
        if sum < sum_temp:
            sum = sum_temp
    """ 
    model = Model(sine_func)
    params = model.make_params(amplitude=1, frequency=q, phase=0, offset=0)
    result = model.fit(data, params, x=correction)  # ✅
    Amplitude =result.params["amplitude"].value
    frequency_sin= result.params["frequency"].value
    Phase =  result.params["phase"].value
    print(Phase)
    print(Amplitude)
    print(frequency_sin)
    """
    w = np.sqrt( 8 /  (Optical_length**2 *  Length_of_hbar)  )
    average_per  =  1/ ((1-Length_of_hbar/Optical_length) + (1/permitivity) * Length_of_hbar/Optical_length)
    photon=   ((permitivity))* average_per* ( h_bar* frequency /2 ) 
    phonon = p_13*np.sqrt(h_bar  * angular_acoustic_frequency / (2 * c_33))
    return sum *  1/h_bar * phonon * photon * w  / (2*np.pi) #important this is divide by 2 pi my dear I need to make sure that In my report it stays that way

def constants():
    length = np.array([
        1.12, 1.00,
        2.29, 2.07,
        4.40, 3.98,
        4.99, 4.53
    ])

    gradient_max = np.array([
        1.406, 1.427,
        0.685, 0.695,
        0.357, 0.362,
        0.314, 0.319,
        0.138, 0.248
    ]) * 1000

    refractive_indices = np.array([
        1.616, 1.616,
        1.616, 1.616,
        1.616, 1.616,
        1.616, 1.616,
        1.515, 1.515
    ])

    return {
        "length": length,
        "grin_index": gradient_max,
        "refractive_index": refractive_indices
    }

def run_single(args):
    r1, r2, Track, wavelength, points = args

    My_track = Track(front_mirror=r1, back_mirror=r2)
    waist_at_the_front_mirror, waist = Find_the_waist_analyticaly(My_track)

    R1_pos = tracks.Track_front_mirror(My_track)
    z0 = Raylaigh_range_1(waist_at_the_front_mirror, wavelength)
    q_int = complex(1/R1_pos, -1/z0)
    q_inv = 1/q_int
    z = q_inv.real

    positions = np.linspace(0, tracks.Find_length_track(My_track), points)
    #
    # This is the function that alters position crystal please turn it off if you want to fix the position of the crstyal
    My_track = tracks.Adjust_the_set_up(My_track, -z)
    #
    #
    waist_at_the_front_mirror, waist = Find_the_waist_analyticaly(My_track)
    wavelengths = wavelength * np.ones(points)
    zr = Raylaigh_range_1(waist_at_the_front_mirror, wavelength)

    beam = Gaussian_Beam(zr, R1_pos, Position=False)
    beam_lists = [copy.deepcopy(beam)]
    wavelengths = Propagate(
        beam_lists, positions, "right", My_track, wavelengths, Refraction_curve=True
    )

    Hbar, name = tracks.Find_the_HBAR(My_track, name=True)
    beams_inside_hbar = []
    sum_waist = 0
    points_inside = 0

    for k, beam in enumerate(beam_lists):
        if Hbar["position_beg"] < positions[k] < Hbar["position_end"]:
            w = beam.beam_waist(wavelengths[k])
            beams_inside_hbar.append(w)
            sum_waist += w
            points_inside += 1

    if len(beams_inside_hbar) > 0:
        av = sum_waist / points_inside
        coupling_sum = Coupling_integral(
            My_track, beams_inside_hbar, Hbar, av, name
        )
        return (r1, r2, av, coupling_sum)

    return (r1, r2, np.nan, np.nan)  # return nans if no beams inside HBAR



def check_waist_symmetry(Track, R1_vals, R2_vals, wavelength=1550e-9, tol=1e-9):
    """
    Check whether the waist is symmetric under swapping R1 ↔ R2.
    Computes waists on the fly — no matrix needed.
    """
    print("\n── Waist symmetry check (swapping R1 ↔ R2) ──")

    is_symmetric = True
    for i in range(len(R1_vals)):
        for j in range(i + 1, len(R2_vals)):

            _, w_ij = Find_the_waist_analyticaly(Track(front_mirror=R1_vals[i], back_mirror=R2_vals[j]))
            _, w_ji = Find_the_waist_analyticaly(Track(front_mirror=R1_vals[j], back_mirror=R2_vals[i]))

            delta = abs(w_ij - w_ji)
            status = "✓ symmetric" if delta < tol else "✗ ASYMMETRIC"
            print(
                f"  (i={i}, j={j})  "
                f"R1={R1_vals[i]*1e3:.0f}mm R2={R2_vals[j]*1e3:.0f}mm → waist={w_ij*1e6:.4f} µm  |  "
                f"swap → waist={w_ji*1e6:.4f} µm  |  Δ={delta*1e9:.3f} nm  {status}"
            )
            if delta >= tol:
                is_symmetric = False

    print(f"\n  Overall: {'matrix IS symmetric' if is_symmetric else 'matrix is NOT symmetric'}")
    print("─────────────────────────────────────────────")
    return is_symmetric

