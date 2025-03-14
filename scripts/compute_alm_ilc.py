import numpy as np
import healpy as hp

def compute_alm_ilc(almmaps:np.ndarray, alm_size:int, wi:np.ndarray, l_max:int):
    r"""
    Computes the weighted spherical harmonic coefficients.

    \(a^{ILC}_{lm} = \sum_i w_i(l)a^i_{lm}\)

    Parameters:
    -------
    almmaps: np.ndarray
        Sky maps converted from pixel space to harmonic space.
    
    alm_size: int
        Number of spherical harmonic coefficients (lmax+1)^2.

    wi : np.ndarray
        Weights.

    l_max: int
        Maximum multipole.
    
    Returns:
    ------
    almilc : np.ndarray
        Shape of (almsize x 1)
    """

    almilc = np.zeros((alm_size), dtype=complex)

    for l in range(l_max +1):
        index= np.zeros((l+1), dtype=int)
        for m in range(l+1):
            index[m]= hp.Alm.getidx(l_max, l, m)
        alm_p = almmaps[:, index]

        alm = np.dot(wi[:,l], alm_p[:,:]) # multiplying (9x384) and (9x1) to get (384x1)=this is for one l (m =0-> l). we iterate over all l
        '''
        #summing over frequencies. and taking ALL m values for a fixed l. Looping over l means 
        so we get an alm weighted by contributions from all frequencies over all scales.
        this should be the same as writing alm_p??
        '''
        almilc[index]= alm #putting alm values back into an array corresponding to their proper positions. flat array
    return almilc