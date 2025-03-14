import numpy as np
import pymaster as nmt

def compute_decoupled_cl_2( maps: np.ndarray, masks:np.ndarray, nsides:int, l_max:int, n_freq:int)-> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the power spectra using NaMaster for a set of frequency maps.

    Parameters:
    - maps: numpy array of shape (n_freqs, n_maps) containing the frequency maps
    - masks: mask to apply to the maps
    - nsides: resolution parameter for the maps
    - l_max: maximum multipole for the power spectra
    - n_freq: number of frequency channels

    Returns:
    - p_cl: raw power spectrum (coupled)
    - matrix: mode coupling matrix
    - cldec: decoupled power spectrum
    """
    p_cl =np.zeros((l_max+1, n_freq, n_freq))
    cldec = np.zeros((l_max+1, n_freq, n_freq))

    f_0 =nmt.NmtField(1-masks, [maps]) #initializing field
    p_cl = nmt.compute_coupled_cell(f_0,f_0) #masked Cl
    # Define a NaMaster binning scheme (no binning)
    b = nmt.NmtBin.from_nside_linear(nsides, 1)
    # Create a NaMaster workspace
    w = nmt.NmtWorkspace.from_fields(f_0, f_0, b)

    # Extract the mode-coupling matrix
    matrix = w.get_coupling_matrix()

    mcm_inv= np.linalg.pinv(matrix) #inverting the mcm
    cldec = np.dot(mcm_inv, p_cl[0]) #decoupling the power spectrum

    return p_cl, matrix, cldec
