import numpy as np 
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def Track(distance_to_the_sample= 4e-3 , R = 1.3  , gap_length = 3e-6 , d4 = 5.5e-3 + 0.07e-3 - 0.009e-3 , mirror_radius_of_curvture = 15e-3):
 
    Vacuum_0 = {
        "type": "medium",
        "n": 1,
        "position_beg": 0,
        "Hbar" : False,
        "position_end": distance_to_the_sample, 
        "surface_left": np.inf,
        "surface_right": np.inf,
    }
    
 
    Quartz = {
        "type": "medium",
        "n": 1.527,
        "Hbar" : True , 
        "position_beg": Vacuum_0["position_end"],
        "position_end": Vacuum_0["position_end"] + 1e-3,  # 1 mm
        "surface_left": np.inf,
        "surface_right": np.inf,
    }
   
 
    Piezo = {
        "type": "medium",
        "n": 2.19,
        "Hbar" : False, 
        "position_beg": Quartz["position_end"],
        "position_end": Quartz["position_end"] + 450e-9,  # 550 nm
        "surface_left": np.inf,
        "surface_right": R,  # positive radius of curvature
    }
    
 
    Vacuum_1 = {
        "type": "medium",
        "n": 1,
        "Hbar" : False,
        "position_beg": Piezo["position_end"],
        "position_end": Piezo["position_end"] + gap_length,  # example spacing
        "surface_left": -R,
        "surface_right": np.inf,
    }
    
 
    chip = {
        "type": "medium",
        "n": 1.74,
        "Hbar" : False,
        "position_beg": Vacuum_1["position_end"],
        "position_end": Vacuum_1["position_end"] + 430e-6,  # 420 µm
        "surface_left": np.inf,
        "surface_right": np.inf,
    }
    
 
    Vacuum_2 = {
        "type": "medium",
        "n": 1,
        "Hbar" : False,
        "position_beg": chip["position_end"],
        "position_end": chip["position_end"] + d4,
        "surface_left": np.inf,
        "surface_right": mirror_radius_of_curvture,
    }
    
 
    return {
        "Vacuum_0": Vacuum_0,
        "Quartz":   Quartz,
        "Piezo":    Piezo,
        "Vacuum_1": Vacuum_1,
        "chip":     chip,
        "Vacuum_2": Vacuum_2,
    }

def Track_concave_concave(distance_to_the_sample= 4e-3 , R = 1.3  , gap_length = 3e-6 , d4 = 5.5e-3 + 0.07e-3 - 0.009e-3 , front_mirror=  -7e-3 , back_mirror = 7e-3  ):
 
    Vacuum_0 = {
        "type": "medium",
        "n": 1,
        "Hbar" : False,
        "position_beg": 0,
        "position_end": distance_to_the_sample, 
        "surface_left":  front_mirror,
        "surface_right": np.inf,
    }
    
 
    Quartz = {
        "type": "medium",
        "n": 1.527,
        "Hbar" : True,
        "position_beg": Vacuum_0["position_end"],
        "position_end": Vacuum_0["position_end"] + 1e-3,  # 1 mm
        "surface_left": np.inf,
        "surface_right": np.inf,
    }
   
 
    Piezo = {
        "type": "medium",
        "n": 2.19,
        "Hbar" : False,
        "position_beg": Quartz["position_end"],
        "position_end": Quartz["position_end"] + 450e-9,  # 550 nm
        "surface_left": np.inf,
        "surface_right": R,  # positive radius of curvature
    }
    
 
    Vacuum_1 = {
        "type": "medium",
        "n": 1,
        "Hbar" : False,
        "position_beg": Piezo["position_end"],
        "position_end": Piezo["position_end"] + gap_length,  # example spacing
        "surface_left": -R,
        "surface_right": np.inf,
    }
    
 
    chip = {
        "type": "medium",
        "n": 1.74,
        "Hbar" : False,
        "position_beg": Vacuum_1["position_end"],
        "position_end": Vacuum_1["position_end"] + 430e-6,  # 420 µm
        "surface_left": np.inf,
        "surface_right": np.inf,
    }
    
 
    Vacuum_2 = {
        "type": "medium",
        "n": 1,
        "Hbar" : False,
        "position_beg": chip["position_end"],
        "position_end": chip["position_end"] + d4,
        "surface_left": np.inf,
        "surface_right": back_mirror,
    }
    
 
    return {
        "Vacuum_0": Vacuum_0,
        "Quartz":   Quartz,
        "Piezo":    Piezo,
        "Vacuum_1": Vacuum_1,
        "chip":     chip,
        "Vacuum_2": Vacuum_2,
    }

def Symmetric_track(distance_to_the_sample = 5.15e-3  , distance_to_the_back_mirror =5.15e-3  , back_mirror = 7e-3 , front_mirror = -7e-3, R =  1.3):
    Vacuum_0 = {
        "type": "medium",
        "n": 1,
        "position_beg": 0,
        "Hbar" : False,
        "position_end": distance_to_the_sample, 
        "surface_left":  front_mirror,
        "surface_right":R,
    }
    Piezo_1 = {
        "type": "medium",
        "n": 2.19,
        "Hbar" : False,
        "position_beg": Vacuum_0["position_end"],
        "position_end": Vacuum_0["position_end"] + 450e-9,  # 550 nm
        "surface_left": -R ,
        "surface_right": np.inf,  # positive radius of curvature
    }
    
 
    Quartz = {
        "type": "medium",
        "n": 1.527,
        "Hbar" : True , 
        "position_beg": Piezo_1["position_end"],
        "position_end": Piezo_1["position_end"] + 1e-3,  # 1 mm
        "surface_left": np.inf,
        "surface_right": np.inf,
    }
    Piezo_2 = {
        "type": "medium",
        "n": 2.19,
        "Hbar" : False,
        "position_beg": Quartz["position_end"],
        "position_end": Quartz["position_end"] + 450e-9,  # 550 nm
        "surface_left": np.inf,
        "surface_right": R,  # positive radius of curvature
            }
    
    
    
    
    Vacuum_1 = {
        "type": "medium",
        "n": 1,
        "Hbar" : False,
        "position_beg": Piezo_2["position_end"],
        "position_end": Piezo_2["position_end"] + distance_to_the_back_mirror,  # example spacing
        "surface_left": -R ,
        "surface_right": back_mirror,
    }

    
    
    
    
    return {
        "Vacuum_0": Vacuum_0,
        "Piezo_1 " : Piezo_1, 
        "Quartz":  Quartz,
        "Piezo_2" : Piezo_2,
        "Vacuum_1": Vacuum_1,
    }

def Tom_2_track():
    Vacuum_0 = {
        "type": "medium",
        "n": 1,
        "position_beg": 0,
        "position_end": 0.0025 ,
        "surface_left": np.inf,
        "surface_right": np.inf,
    }
    
 
    Quartz = {
        "type": "medium",
        "n": 1.54,
        "position_beg": Vacuum_0["position_end"],
        "position_end": Vacuum_0["position_end"] + 1e-3,  # 1 mm
        "surface_left": np.inf,
        "surface_right": np.inf,
    }
   
 
    Piezo = {
        "type": "medium",
        "n": 2.19,
        "position_beg": Quartz["position_end"],
        "position_end": Quartz["position_end"] + 41/ 100000000,  # 550 nm
        "surface_left": np.inf,
        "surface_right": 2/125,  # positive radius of curvature
    }
    
 
    Vacuum_1 = {
        "type": "medium",
        "n": 1,
        "position_beg": Piezo["position_end"],
        "position_end": Piezo["position_end"] + 1/ 500000,  # example spacing
        "surface_left": -2/125,
        "surface_right": np.inf,
    }
    
 
    chip = {
        "type": "medium",
        "n": 1.74,
        "position_beg": Vacuum_1["position_end"],
        "position_end": Vacuum_1["position_end"] + 21/ 50000,  # 420 µm
        "surface_left": np.inf,
        "surface_right": np.inf,
    }
    
 
    Vacuum_2 = {
        "type": "medium",
        "n": 1,
        "position_beg": chip["position_end"],
        "position_end": chip["position_end"] + 0.0093363,
        "surface_left": np.inf,
        "surface_right": 16e-3,
    }
    return {
        "Vacuum_0": Vacuum_0,
        "Quartz":   Quartz,
        "Piezo":    Piezo,
        "Vacuum_1": Vacuum_1,
        "chip":     chip,
        "Vacuum_2": Vacuum_2,
    }

def Tom_Track():
 
    Vacuum_0 = {
        "type": "medium",
        "n": 1,
        "position_beg": 0,
        "position_end": 2.5e-3,
        "surface_left": np.inf,
        "surface_right": np.inf,
    }
    
 
    Quartz = {
        "type": "medium",
        "n": 1.54,
        "position_beg": Vacuum_0["position_end"],
        "position_end": Vacuum_0["position_end"] + 1e-3,  # 1 mm
        "surface_left": np.inf,
        "surface_right": np.inf,
    }
    
 
    Piezo = {
        "type": "medium",
        "n": 2.19,
        "position_beg": Quartz["position_end"],
        "position_end": Quartz["position_end"] + 41/100000000,  # 410 nm
        "surface_left": np.inf,
        "surface_right": 9/20,  # positive radius of curvature
    }
    
 
    Vacuum_1 = {
        "type": "medium",
        "n": 1,
        "position_beg": Piezo["position_end"],
        "position_end": Piezo["position_end"] + 1/500000,  # example spacing
        "surface_left":  -9/20,
        "surface_right": np.inf,
    }
    
 
    chip = {
        "type": "medium",
        "n": 1.74,
        "position_beg": Vacuum_1["position_end"],
        "position_end": Vacuum_1["position_end"] + 21/50000,  # 420 µm
        "surface_left": np.inf,
        "surface_right": np.inf,
    }
   
 
    Vacuum_2 = {
        "type": "medium",
        "n": 1,
        "position_beg": chip["position_end"],
        "position_end": chip["position_end"] + 9.3363e-3,
        "surface_left": np.inf,
        "surface_right": 16e-3,
    }
    
 
    return {
        "Vacuum_0": Vacuum_0,
        "Quartz":   Quartz,
        "Piezo":    Piezo,
        "Vacuum_1": Vacuum_1,
        "chip":     chip,
        "Vacuum_2": Vacuum_2,
    }

def GRIN_Track_Empty(g1, l1,  d1 ,d  , R  = 15e-3):
    Vacuum_0 = {
        'type' : "medium",
        "n": 1,
        "position_beg": 0,
        "position_end": d1, 
        "surface_left": np.inf,
        "surface_right": np.inf
    }
    
    Grin_lens1 = {
        'type' : "lens",
        'g1' : g1,
        'l1'  : l1,
        "position_beg": Vacuum_0["position_end"],
        "position_end": Vacuum_0["position_end"] + l1, 
        "surface_left": np.inf,
        "surface_right": np.inf,
        "n" : 1.616,
    }
    Vacuum_1 = {
        'type' : "medium",
        "n": 1,
        "position_beg": Grin_lens1["position_end"],
        "position_end": Grin_lens1["position_end"] + d + 10e-3,
        "surface_left": np.inf,
        "surface_right": np.inf
        }
   
    return {
        "Vacuum_0" : Vacuum_0,
        "Grin_lens1" : Grin_lens1,
        "Vacuum_1" : Vacuum_1
    }

def GRIN_Track_Mirror(g1, l1,  d1 ,d   , refractive_index = 1.515, R  = 15e-3):

    Vacuum_0 = {
        'type' : "medium",
        "n": 1,
        "position_beg": 0,
        "position_end": d1, 
        "surface_left": np.inf,
        "surface_right": np.inf
    }
    
    Grin_lens1 = {
        'type' : "lens",
        'g1' : g1,
        'l1'  : l1,
        "position_beg": Vacuum_0["position_end"],
        "position_end": Vacuum_0["position_end"] + l1, 
        "surface_left": np.inf,
        "surface_right": np.inf,
        "n" : refractive_index

    }
    Vacuum_1 = {
        'type' : "medium",
        "n": 1,
        "position_beg": Grin_lens1["position_end"],
        "position_end": Grin_lens1["position_end"] + d,
        "surface_left": np.inf,
        "surface_right": np.inf
        }
    Mirror_concave = {
        'type': "medium",
        "n": 1.45 ,
        "position_beg": Vacuum_1["position_end"],  # fixed typo: postion_beg -> position_beg
        "position_end": Vacuum_1["position_end"] + 4e-3,
        "surface_left": np.inf , ##### this is wrong in the case of the mirror 
        "surface_right": R #####
    }
    
    return {
        "Vacuum_0" : Vacuum_0,
        "Grin_lens1" : Grin_lens1,
        "Vacuum_1" : Vacuum_1,
        "Mirror_concave" : Mirror_concave
    }

def Find__optical_length_track(Track):
    distance = 0
    for keys in Track:
        element = Track[keys]
        d = abs(element["position_beg"] -  element["position_end"])
        
        if element["type"] == "medium":
            distance = d*element["n"] + distance
        else:
            distance = d + distance
        logger.debug(f"The distance is {d}")
    return distance

def Find_length_track(Track):
    distance = 0
    for keys in Track:
        element = Track[keys]
        d = abs(element["position_beg"] -  element["position_end"])
        
        if element["type"] == "medium":
            distance = d + distance
        else:
            distance = d + distance
        logger.debug(f"The distance is {d}")
    return distance

def GRIN_Track(g1, g2, l1, l2 , d , include_second_lens = False ):
    if include_second_lens == False:

        Vacuum_0 = {
        'type' : 'medium',
        "n": 1,
        "position_beg": 0,
        "position_end": d, 
        "surface_left": np.inf,
        "surface_right": np.inf
    }
    
        Grin_lens1 = {
        'type' : 'lens',
        'g1' : g1,
        'l1'  : l1,
        "position_beg": Vacuum_0["position_end"],
        "position_end": Vacuum_0["position_end"] + l1, 
        "surface_left": np.inf,
        "surface_right": np.inf

    }
        Vacuum_1 = {
        'type' : 'medium',
        "n": 1,
        "position_beg": Grin_lens1["position_end"],
        "position_end": Grin_lens1["position_end"] + 2e-3,
        "surface_left": np.inf,
        "surface_right": np.inf
        }
        
        return {
        "Vacuum_0" : Vacuum_0,
        "Grin_lens1" : Grin_lens1,
        "Vacuum_1" : Vacuum_1
    }

    else:
        Vacuum_0 = {
        'type' : 'medium',
        "n": 1,
        "position_beg": 0,
        "position_end": d, 
        "surface_left": np.inf,
        "surface_right": np.inf
    }
    
        Grin_lens1 = {
        'type' : 'lens',
        'g1' : g1,
        'l1'  : l1,
        "position_beg": Vacuum_0["position_end"],
        "position_end": Vacuum_0["position_end"] + l1, 
        "surface_left": np.inf,
        "surface_right": np.inf

    }
        Grin_lens2 = {
            'type': 'lens',
            'g1' : g2,
            'l1'  : l2,
            "position_beg": Grin_lens1["position_end"],
            "position_end": Grin_lens1["position_end"] + l2, 
            "surface_left": np.inf,
            "surface_right": np.inf
        }
        return {
            "Vacuum_0" : Vacuum_0,
            "Grin_lens1" : Grin_lens1,
            "Grin_lens2" : Grin_lens2
        }
    
def Concave_conve():
    Mirror_1 = {

    }
    Vacuum_0  =  {

    }
    Quartz = {
        "type": "medium",
        "n": 1.527,
        "position_beg": Vacuum_0["position_end"],
        "position_end": Vacuum_0["position_end"] + 1e-3,  # 1 mm
        "surface_left": np.inf,
        "surface_right": np.inf,
    }
     
    Piezo = {
        "type": "medium",
        "n": 2.19,
        "position_beg": Quartz["position_end"],
        "position_end": Quartz["position_end"] + 450e-9,  # 550 nm
        "surface_left": np.inf,
        "surface_right": 450e-3,  # positive radius of curvature
    }
    
 
    Vacuum_1 = {
        "type": "medium",
        "n": 1,
        "position_beg": Piezo["position_end"],
        "position_end": Piezo["position_end"] + 3e-6,  # example spacing
        "surface_left": -450e-3,
        "surface_right": np.inf,
    }
    
 
    chip = {
        "type": "medium",
        "n": 1.74,
        "position_beg": Vacuum_1["position_end"],
        "position_end": Vacuum_1["position_end"] + 430e-6,  # 420 µm
        "surface_left": np.inf,
        "surface_right": np.inf,
    }
    Vacuum_1 = {

    }
    Mirror_2 = { 

    }

def Track_back_mirror(Track , whole = False):
    keys = list(Track.keys())
    last = keys[-1]
    if whole  == False:
        return Track[last]["surface_right"]
    else:
        return Track[last]
    
def Track_front_mirror(Track):
    keys = list(Track.keys())
    first = keys[0]
    return Track[first]["surface_left"]

def Track_length_adjuster(Track):
    distance = 0
    
    for keys in Track:
        element = Track[keys]
        d = abs(element["position_beg"] -  element["position_end"])
        
        if element["type"] == "medium":
            distance = d*element["n"] + distance
        else:
            distance = d + distance
        logger.debug(f"The distance is {d}")
    to_the_mirror= distance - d 
    return to_the_mirror

def Track_presenter(Track):
    for n, key in enumerate(Track.keys()):
        element = Track[key]
        distance = abs(element["position_beg"] - element["position_end"])
        logger.info("The name of the element nr %d and its name is %s", n, key)
        logger.info("The beginning is %e and the end %e and the distance between them is %e", element["position_beg"], element["position_end"], distance)
        if element["type"] == "medium":
            if element["surface_right"] == np.inf:
                logger.info("The right surface is flat")
            elif element["surface_right"] > 0:
                logger.info("The right surface is concave")
            else:
                logger.info("The right surface is convex")
            if element["surface_left"] == np.inf:
                logger.info("The left surface is flat")
            elif element["surface_left"] > 0:
                logger.info("The left surface is concave")
            else:
                logger.info("The left surface is convex")
        else:
            logger.info("The g factor of the lens is %e and the length %e", element["g1"], element["l1"])

def Track_save(Track, filepath):
    with open(filepath, "w") as f:
        for n, key in enumerate(Track.keys()):
            element = Track[key]
            distance = abs(element["position_beg"] - element["position_end"])
            f.write(f"The name of the element nr {n} and its name is {key}\n")
            f.write(f"The beginning is {element['position_beg']:.6e} and the end {element['position_end']:.6e} and the distance between them is {distance:.6e}\n")
            if element["type"] == "medium":
                if element["surface_right"] == np.inf:
                    f.write("The right surface is flat\n")
                elif element["surface_right"] > 0:
                    f.write("The right surface is concave\n")
                else:
                    f.write("The right surface is convex\n")
                if element["surface_left"] == np.inf:
                    f.write("The left surface is flat\n")
                elif element["surface_left"] > 0:
                    f.write("The left surface is concave\n")
                else:
                    f.write("The left surface is convex\n")
            else:
                f.write(f"The g factor of the lens is {element['g1']:.6e} and the length {element['l1']:.6e}\n")
            f.write("\n")

def Calculate_the_briouilin_frequnecy(n , wavelength = 1550e-9):
    vl = 6400 # longitudal speed of sound
    c = 3*1e8 #speed of light to calculate the frequency of the laser
    cn= 3*1e8

    frequency  =cn / wavelength  
    briouilin_frequnecy =(2* frequency )  * vl  *n  / (c)
    return  briouilin_frequnecy

def Adjust_the_midle_set_up(Track, Z):
    """
    Shift the optical track so that the beam waist coincides with the antenna position( in this case in the middle of the chip).

    The waist position `Z` is computed elsewhere (see [name of the function]).
    This function translates all elements of the track by the offset between the
    waist and the antenna, so that the waist ends up at the antenna position.

    Parameters
    ----------
    Track : [type, e.g. list of dict / np.ndarray]
    Cavity track containing the positions of all elements
    (mirrors, crystal, etc.).
    Z : float
    Position of the beam waist along the optical axis [unit, e.g. mm].

    Returns
    Track
    [type]
    The shifted track, with the beam waist located at the antenna position.
        """
    key_tracker_before = []
    key_tracker_after = []
    chip_elem = None
    chip_found = False
 
    total_length  = sum(elem["position_end"] - elem["position_beg"] for elem in Track.values())
    total_optical = Find__optical_length_track(Track)
 
    for key in Track.keys():
        elem = Track[key]
        if key == "chip":
            chip_elem = elem
            chip_found = True
            continue
        if not chip_found:
            key_tracker_before.append(key)
        else:
            key_tracker_after.append(key)
 
    if chip_elem is None:
        logger.warning("There's no chip inside our cavity")
        return Track
 
    length_chip         = chip_elem["position_end"] - chip_elem["position_beg"]
    optical_length_chip = length_chip * chip_elem["n"]
    half_optical_chip   = optical_length_chip / 2
 
    optical_fixed_before = sum(
        (Track[k]["position_end"] - Track[k]["position_beg"]) * Track[k]["n"]
        for k in key_tracker_before if k != "Vacuum_0"
    )
    optical_fixed_after = sum(
        (Track[k]["position_end"] - Track[k]["position_beg"]) * Track[k]["n"]
        for k in key_tracker_after if k != "Vacuum_2"
    )
 
    # Z is now the optical position of the chip's MIDDLE (not position_beg)
    optical_vacuum_before = Z - optical_fixed_before - half_optical_chip
    optical_vacuum_after  = (total_optical - Z) - half_optical_chip - optical_fixed_after
 
    if optical_vacuum_before < 0:
        logger.warning(
            f"Not enough optical space before chip at Z={Z:.4e}: "
            f"need {optical_fixed_before + half_optical_chip:.4e} but only have {Z:.4e}"
        )
    if optical_vacuum_after < 0:
        logger.warning(
            f"Not enough optical space after chip at Z={Z:.4e}: "
            f"need {optical_fixed_after + half_optical_chip:.4e} "
            f"but only have {total_optical - Z:.4e}"
        )
 
    n_vac = Track["Vacuum_0"]["n"]
    Track["Vacuum_0"]["_new_length"] = optical_vacuum_before / n_vac
    Track["Vacuum_2"]["_new_length"] = optical_vacuum_after  / n_vac
 
    # Z maps to chip midpoint; position_beg is derived from that below
    physical_fixed_before  = sum(
        Track[k]["position_end"] - Track[k]["position_beg"]
        for k in key_tracker_before if k != "Vacuum_0"
    )
    physical_vacuum_before = optical_vacuum_before / n_vac
 
    chip_elem["position_beg"] = physical_vacuum_before + physical_fixed_before
    chip_elem["position_end"] = chip_elem["position_beg"] + length_chip
 
    # Reposition elements BEFORE chip, walking backwards
    prev = chip_elem
    for key in reversed(key_tracker_before):
        elem = Track[key]
        d = elem.pop("_new_length", elem["position_end"] - elem["position_beg"])
        elem["position_end"] = prev["position_beg"]
        elem["position_beg"] = prev["position_beg"] - d
        prev = elem
 
    # Force first element to 0
    if key_tracker_before:
        first_elem = Track[key_tracker_before[0]]
        offset = first_elem["position_beg"]
        if not np.isclose(offset, 0):
            logger.info(f"Shifting entire track by {-offset:.6e} m to align first element to 0")
            for key in Track.keys():
                Track[key]["position_beg"] -= offset
                Track[key]["position_end"] -= offset
 
    # Reposition elements AFTER chip, walking forwards
    prev = chip_elem
    for key in key_tracker_after:
        elem = Track[key]
        d = elem.pop("_new_length", elem["position_end"] - elem["position_beg"])
        elem["position_beg"] = prev["position_end"]
        elem["position_end"] = prev["position_end"] + d
        prev = elem
 
    # Sanity checks
    final_length  = Track[key_tracker_after[-1]]["position_end"] if key_tracker_after else chip_elem["position_end"]
    final_optical = Find__optical_length_track(Track)
 
    if not np.isclose(final_length, total_length):
        logger.warning(f"Physical length changed: was {total_length:.6e}, now {final_length:.6e}")
    if not np.isclose(final_optical, total_optical):
        logger.warning(f"Optical length changed: was {total_optical:.6e}, now {final_optical:.6e}")
 
    first_key = key_tracker_before[0] if key_tracker_before else "chip"
    if not np.isclose(Track[first_key]["position_beg"], 0):
        logger.warning(
            f"First element '{first_key}' does not start at 0: "
            f"{Track[first_key]['position_beg']:.6e}"
        )
 
    return Track

def Adjust_the_set_up(Track, Z):
    """
    Shift the optical track so that the beam waist coincides with the antenna position( in this in the front of the chip).

    The waist position `Z` is computed elsewhere (see [name of the function]).
    This function translates all elements of the track by the offset between the
    waist and the antenna, so that the waist ends up at the antenna position.

    Parameters
    ----------
    Track : [type, e.g. list of dict / np.ndarray]
    Cavity track containing the positions of all elements
    (mirrors, crystal, etc.).
    Z : float
    Position of the beam waist along the optical axis [unit, e.g. mm].

    Returns
    Track
    [type]
    The shifted track, with the beam waist located at the antenna position.
    """


    key_tracker_before = []
    key_tracker_after = []
    chip_elem = None
    chip_found = False

    total_length  = sum(elem["position_end"] - elem["position_beg"] for elem in Track.values())
    total_optical = Find__optical_length_track(Track)

    for key in Track.keys():
        elem = Track[key]
        if key == "chip":
            chip_elem = elem
            chip_found = True
            continue
        if not chip_found:
            key_tracker_before.append(key)
        else:
            key_tracker_after.append(key)

    if chip_elem is None:
        logger.warning("There's no chip inside our cavity")
        return Track

    length_chip         = chip_elem["position_end"] - chip_elem["position_beg"]
    optical_length_chip = length_chip * chip_elem["n"]

    optical_fixed_before = sum(
        (Track[k]["position_end"] - Track[k]["position_beg"]) * Track[k]["n"]
        for k in key_tracker_before if k != "Vacuum_0"
    )
    optical_fixed_after = sum(
        (Track[k]["position_end"] - Track[k]["position_beg"]) * Track[k]["n"]
        for k in key_tracker_after if k != "Vacuum_2"
    )

    # Z is now the optical position of the chip's position_beg (no half-chip offset)
    optical_vacuum_before = Z - optical_fixed_before
    optical_vacuum_after  = (total_optical - Z) - optical_length_chip - optical_fixed_after

    if optical_vacuum_before < 0:
        logger.warning(
            f"Not enough optical space before chip at Z={Z:.4e}: "
            f"need {optical_fixed_before:.4e} but only have {Z:.4e}"
        )
    if optical_vacuum_after < 0:
        logger.warning(
            f"Not enough optical space after chip at Z={Z:.4e}: "
            f"need {optical_fixed_after + optical_length_chip:.4e} "
            f"but only have {total_optical - Z:.4e}"
        )

    n_vac = Track["Vacuum_0"]["n"]
    Track["Vacuum_0"]["_new_length"] = optical_vacuum_before / n_vac
    Track["Vacuum_2"]["_new_length"] = optical_vacuum_after  / n_vac

    # Z maps directly to position_beg of chip — no half-length arithmetic needed
    physical_fixed_before  = sum(
        Track[k]["position_end"] - Track[k]["position_beg"]
        for k in key_tracker_before if k != "Vacuum_0"
    )
    physical_vacuum_before = optical_vacuum_before / n_vac

    chip_elem["position_beg"] = physical_vacuum_before + physical_fixed_before
    chip_elem["position_end"] = chip_elem["position_beg"] + length_chip

    # Reposition elements BEFORE chip, walking backwards
    prev = chip_elem
    for key in reversed(key_tracker_before):
        elem = Track[key]
        d = elem.pop("_new_length", elem["position_end"] - elem["position_beg"])
        elem["position_end"] = prev["position_beg"]
        elem["position_beg"] = prev["position_beg"] - d
        prev = elem

    # Force first element to 0
    if key_tracker_before:
        first_elem = Track[key_tracker_before[0]]
        offset = first_elem["position_beg"]
        if not np.isclose(offset, 0):
            logger.info(f"Shifting entire track by {-offset:.6e} m to align first element to 0")
            for key in Track.keys():
                Track[key]["position_beg"] -= offset
                Track[key]["position_end"] -= offset

    # Reposition elements AFTER chip, walking forwards
    prev = chip_elem
    for key in key_tracker_after:
        elem = Track[key]
        d = elem.pop("_new_length", elem["position_end"] - elem["position_beg"])
        elem["position_beg"] = prev["position_end"]
        elem["position_end"] = prev["position_end"] + d
        prev = elem

    # Sanity checks
    final_length  = Track[key_tracker_after[-1]]["position_end"] if key_tracker_after else chip_elem["position_end"]
    final_optical = Find__optical_length_track(Track)

    if not np.isclose(final_length, total_length):
        logger.warning(f"Physical length changed: was {total_length:.6e}, now {final_length:.6e}")
    if not np.isclose(final_optical, total_optical):
        logger.warning(f"Optical length changed: was {total_optical:.6e}, now {final_optical:.6e}")

    first_key = key_tracker_before[0] if key_tracker_before else "chip"
    if not np.isclose(Track[first_key]["position_beg"], 0):
        logger.warning(
            f"First element '{first_key}' does not start at 0: "
            f"{Track[first_key]['position_beg']:.6e}"
        )

    return Track

def FSR(wavelength = 1550e-9):
    speed_of_sound =  6420 # speed of sound in Aln internet value
    frequency = 12.65 *1e9 # The value that we phonon frequnecy that we target 
    c = 3*1e8
    br = Calculate_the_briouilin_frequnecy( 1.53 , wavelength )
    wavelneth_in_piezo = speed_of_sound/br
    L =(c / (2*br)) 

    return br
    
def Find_the_HBAR(Track, name  = False):
    for key in Track.keys():
        elem = Track[key]
        if elem["Hbar"] == True:
            if name == False:
                return elem
            elif name == True:
                return elem, key
    logger.warning("There's no Hbar inside of our cavity")
    return None

def Acoustic_wave_vector():
    vl  =6400
    br  =  Calculate_the_briouilin_frequnecy(1.53 )
    angular_frequency = 2 * np.pi *br
    q = (angular_frequency/vl)

    return q 

def Distance_to_Hbar(Track):
    distance = 0
    for key in Track.keys():
        elem = Track[key]
        if elem["Hbar"] == True:
               return distance
        distance = elem["position_end"] - elem["position_beg"]
    
                
    logger.warning("There's no Hbar inside of our cavity")
    return None

def Material_constants(name):
    #Crystall Quartz I didn't want to change the entire code so 
    Quartz= {
    "stiffness_tensor": {
        "1_1": 86.76, "1_2": 6.87,  "1_3": 11.85, "1_4": -18.02, "1_5": 0.00, "1_6": 0.00,
        "2_1": 6.87,  "2_2": 86.76, "2_3": 11.85, "2_4": 18.02,  "2_5": 0.00, "2_6": 0.00,
        "3_1": 11.85, "3_2": 11.85, "3_3": 105.46,"3_4": 0.00,   "3_5": 0.00, "3_6": 0.00,
        "4_1": -18.02,"4_2": 18.02, "4_3": 0.00,  "4_4": 58.14,  "4_5": 0.00, "4_6": 0.00,
        "5_1": 0.00,  "5_2": 0.00,  "5_3": 0.00,  "5_4": 0.00,   "5_5": 58.14,"5_6": -18.02,
        "6_1": 0.00,  "6_2": 0.00,  "6_3": 0.00,  "6_4": 0.00,   "6_5": -18.02,"6_6": 39.95
    },

    "permittivity": 2.35,   # for 1550 nm wavelength
    "refractive_index": 1.53,
    "permeability": 1.00,
    "piezo_tensor": {
        "1_1": 0.16,
        "1_2": 0.27,
        "1_3": 0.27,
        "1_4": -0.030,
        "3_1": 0.29,
        "3_3": 0.10,
        "4_1": -0.047,
        "4_4": -0.079
}
    }
    Coupling_mirror=  {
        "Reflexivity" :  0.997,
        
    }
    Coupling_mirror["Coupling"] = 1 - Coupling_mirror["Reflexivity"]

    if name == "Quartz":
        return Quartz
    if name == "Coupling":
        return Coupling_mirror