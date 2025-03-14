# ILC in pixel space
Internal Linear Combination of simulated CMB maps and foregrounds in pixel space
Linear combination of component maps:

$$
d_i(p) = a_i s_{CMB}(p) +f_i(p) +n_i(p)
$$

map, cmb, foregrounds, noise respectively

Assumes:

- CMB spectral energy distribution is a perfect blackbody
- CMB is uncorrelated to other sky components

Affected by errors dramatically. The ILC of the CMB is obtained by linearly combining input maps with frequency dependent weights. 

The weights:

- preserve CMB blackbody spectrum
- minimise variance

$$
s_{ILC} (p) = \sum_i w_i . d_i(p)
$$

with $\sum_i w_i .a_i =1$

$$
s_{ILC} (p) = s_{CMB}(p) + \sum_i w_i . (f_i(p) + n_i(p))
$$

$$
\Bigl \langle s^2_{ILC} \Bigr \rangle = \Bigl \langle s^2_{CMB} \Bigr \rangle +  \sum_{i,i^{\prime}} w_i (\langle f_i f_i ^{\prime}\rangle + \langle n_i n_i ^{\prime} \rangle ) w_i^{\prime}
$$

## ILC method

- Compute the covariance $C_{i i^{\prime}} = \langle d_i (p)d_{i^{\prime}}(p)\rangle$
- Compute the inverse of the covariance
- The CMB spectral energy distribution

$$
A_{CMB} = [a_0, a_1,...,a_8]=[1,1,...,1]
$$

- ILC weights can be computed as:

$$
w = \frac{A^T_{CMB} C^{-1}}{A^T_{CMB} C^{-1}A_{CMB}}
$$

- Combine the input maps with the weights:

$$
\hat{S}_{ILC}(p) = \sum_i w_i . d_i(p)
$$

# ILC in harmonic space
Internal Linear Combination of simulated CMB maps and foregrounds in harmonic space
Performing ILC in harmonic space (as opposed to pixel space) allows for better control over spatial scales using spherical harmonic decomposition.

- **Transform Frequency Maps to Harmonic Space**:

Each map $M_i(\hat{n})$ (where $i$ represents a frequency channel and $\hat{n}$ represents a sky direction) is decomposed into spherical harmonics:

$$
M_i(\hat{n}) = \sum_{l,m} a^i_{lm}Y_{lm}(\hat{n})
$$

where $a^i_{lm}$ are the harmonic coefficients

- Compute Covariance Matrix

covariance matrix of the harmonic coefficients:

``` math
C^{ij}_l = \langle a^i_{lm} a^{j*}_{lm} \rangle
```

where $i$, $j$ are the frequency channels

- Compute weights that minimize variance of the combined map while preserving the CMB signal:

```math
\omega_i(l) = \frac{\sum_j (C_l^{-1})_{ij} b_{j}}{\sum_{ij}b_i (C_l^{-1})_{ij} b_{j}}
```

where $b_i$  represents the CMB transfer function. $b_i = [1,...,1]$

- Combine maps

Construct the harmonic coefficients of the cleaned map.

```math 
a^{ILC}_{lm} = \sum_i w_i(l)a^i_{lm}
 ```

- Reconstruct the map
- Compute the power spectrum
