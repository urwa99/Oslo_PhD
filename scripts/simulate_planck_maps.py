import numpy as np
import healpy as hp
import pysm3
import pysm3.units as u

#simulating foregroud emission and noise at the different Planck frequencies
def simulate_planck_maps(nsides,freq):
    """
    Simulate foreground emission and noise at different Planck frequencies.
    
    Parameters:
    nside (int): Healpix resolution parameter.
    
    Returns:
    tuple: Numpy arrays of simulated maps for noise, dust, synchrotron, free-free, and CMB components.
    """
    # freq = np.array([28.4,  44.1,  70.4,  100.0,  143.0,  217.0,  353.0,  545.0,  857.0])

    
    sens = np.array([150.,  162.,  210.,  77.4,  33.,  46.8, 154,  42.,  5016.])
    
    # Initialize PySM sky models for different foregrounds and CMB
    sky_d = pysm3.Sky(nside=nsides, preset_strings=["d1"])
    sky_s = pysm3.Sky(nside=nsides, preset_strings=["s1"])
    sky_f = pysm3.Sky(nside=nsides, preset_strings=["f1"])
    sky_cmb = pysm3.Sky(nside=nsides, preset_strings=["c1"])
    
       # Initialize storage for maps
    noise_pl, dust_pl, sync_pl, ff_pl, cmb_pl = [], [], [], [], []
    
    for nf, freq in enumerate(freq):
        # Generate noise
        noise = np.random.normal(size=(12 * nsides**2)) * sens[nf] / hp.nside2resol(nsides, True)
        noise_pl.append(noise)
        
        # Convert emissions to uK_CMB
        conversion = u.K_RJ.to(u.K_CMB, equivalencies=u.cmb_equivalencies(freq * u.GHz))
        
        dust_pl.append(sky_d.get_emission(freq * u.GHz)[0] * conversion)
        sync_pl.append(sky_s.get_emission(freq * u.GHz)[0] * conversion)
        ff_pl.append(sky_f.get_emission(freq * u.GHz)[0] * conversion)
        cmb_pl.append(sky_cmb.get_emission(freq * u.GHz)[0] * conversion)
    
    # Convert lists to numpy arrays (N_freq, N_pixels)
    return  np.array(noise_pl), np.array(dust_pl), np.array(sync_pl), np.array(ff_pl), np.array(cmb_pl)
