import numpy as np
import healpy as hp

def compute_cl(almmaps:np.ndarray, l_max: int, n_freq: int)-> tuple[np.ndarray, np.ndarray]:
    r"""
    Compute the power spectrum C_l from alm_maps.

    The power spectrum is computed as:

    \[ \hat{C}_{l}^{ij} = \frac{1}{2l+1}\sum_{m=-l}^{+l} a_{lm}^{i}a^{j \dagger }_{lm} \]
    
    Parameters:
    ----------
    almmaps : np.ndarray
        Array of shape (n_freq, num_alm) containing spherical harmonic coefficients.
    
    l_max : int
        Maximum multipole lmax.
    
    n_freq : int
        Number of frequency channels.
    
    Returns:
    -------
    alm_p : np.ndarray
        Truncated spherical harmonic values for m>0.
    
    C_l : np.ndarray
        Power spectrum array of shape (lmax+1, n_freq, n_freq).
    """
    Cl = np.zeros((l_max + 1, n_freq, n_freq), dtype=float)
    
    for l in range(l_max + 1):
        index = np.zeros((l + 1), dtype=int)
        
        for m in range(l + 1):
            index[m] = hp.Alm.getidx(l_max, l, m)  # Get Alm index
        
        almp = almmaps[:, index]  # Extract relevant alm values
        
        # Compute C_l
        Cl[l, :, :] = np.real(np.outer(almp[:, 0], almp[:, 0]))
        
        for m in range(1, l + 1):  # Sum over m
            Cl[l, :, :] += 2 * np.real(np.outer(almp[:, m], np.conj(almp[:, m])))
    
    return almp, Cl