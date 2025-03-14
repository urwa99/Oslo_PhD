import numpy as np

def compute_weights(n_freqs:int, l_max=int, Cl= np.ndarray)-> tuple[np.ndarray, np.ndarray]:

    r"""
    Computes the weights for the covariance matrix.

    \(\omega_i(l) = \frac{\sum_j (C_l^{-1})_{ij} b_{j}}{\sum_{ij}b_i (C_l^{-1})_{ij} b_{j}}\)

    Parameters:
    -------
    n_freqs: int
      Number of frequency channels.
    l_max: int
      Maximum number of multipole lmax.
    Cl: np.ndarray
      Covariance matrix.
    
    Returns:
    -------
      Cl_inv : np.ndarray
        Pseudo inverse of the Cl_matrix.
      
      wi : np.ndarray
        Weights dependent on frequency and multipole of shape (n_freqs x m) 
    """
    #CMB spectral energy density
    b= np.ones(n_freqs, dtype=float ) #array filled with values of 1. column vector

    wi= np.zeros((n_freqs, l_max+1), dtype=float) #weights

    Cl_inv =np.zeros((n_freqs, n_freqs), dtype=float)

    for l in range(l_max+1):

        Cl_matrix =   Cl[l,:,:]/(2*l +1) #cov matrix for a specific l (nfreqs x nfreqs)
        

        Cl_inv[:,:] = np.linalg.pinv(Cl_matrix[:,:]) #(nfreqs x nfreqs)
        numerator = np.dot(Cl_inv,  b) #nfreqs x 1
        b_trans = b.T
        denominator = np.dot(b_trans, np.dot( Cl_inv, b)) #scalar?? 1x1


        for nf in range(n_freqs):
            wi[:, l] = numerator[:]/denominator# (nfreqs x m) Since this is inside l loop it is only the weight for a specific l. looping will give us the correct shape?
            #w[:,l]= np.maximum(w[:,l], 0)  # Force weights to be non-negative
            wi[:,l] /= np.sum(wi[:,l], axis=0, keepdims=True)  # Normalize to sum to 1
    return Cl_inv, wi