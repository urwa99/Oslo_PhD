import numpy as np
import healpy as hp

def mask_maps(map:np.ndarray, masks:np.ndarray, almsize:int, l_max:int, n_freq:int, n_pix:int)->tuple[np.ndarray, np.ndarray]:
    """
    Masks input maps and then converts them to harmonic space.

    Parameters:
    - maps: numpy array of shape (n_freqs, n_maps) containing the frequency maps
    - masks: mask to apply to the maps
    - l_max: maximum multipole for the power spectra
    - n_freq: number of frequency channels
    - n_pix: number of pixels

    Returns:
    - masked: masked pixel space maps
    - almmasked: masked maps converted to harmonic space
    
    """
    masked =np.zeros((n_freq, n_pix))
    almmasked = np.zeros((n_freq, almsize), dtype=complex)
    for nf in range(n_freq):
        masked[nf,:] =map[nf,:]*(1-masks)
        almmasked[nf,:]= hp.map2alm(masked[nf,:], lmax=l_max, mmax=None, iter=0, pol=False)
    return masked, almmasked