import numpy as np
import healpy as hp

def cosine_smoothed_mask(nside):
    """Creates a HEALPix mask that smoothly transitions along the Galactic plane as a cosine."""
    
    npix = hp.nside2npix(nside)  # Number of HEALPix pixels
    theta, phi = hp.pix2ang(nside, np.arange(npix))  # Get theta, phi for all pixels
    
    # Convert colatitude to Galactic latitude b
    #HEALPix colatitude 𝜃 to Galactic latitude b

    b = 90 - np.degrees(theta)
    bmin= 2
    bmax=10

    # Initialize mask
    mask = np.zeros(npix)

    # Apply cosine smoothing in the range 8° < |b| < 10°
    transition = (np.abs(b) > bmin) & (np.abs(b) < bmax)
    mask[transition] = 0.5 * (1 - np.cos(np.pi * (np.abs(b[transition]) - bmin) / (bmax-bmin)))

    # Set mask to 1 for |b| ≥ 10°
    mask[np.abs(b) >= bmax] = 1

    return mask