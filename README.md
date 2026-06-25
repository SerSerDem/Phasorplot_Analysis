# Phasorplot_Analysis

This repository contains a Python script for spectral phasor analysis of temperature-dependent fluorescence emission spectra.

The script imports fluorescence emission spectra from CSV files, performs spectral normalization, applies wavelength selection (400–550 nm), and converts the processed spectra into spectral phasor coordinates (G, S). The phasor transformation enables a model-free representation of spectral shifts associated with changes in membrane environments.

For each dataset, the workflow generates normalized emission spectra, calculates temperature-resolved phasor trajectories, and applies interpolation to improve comparability across temperature series. Both spectral and phasor data are exported as CSV files for further quantitative analysis. In addition, the script provides visualization of normalized emission spectra and corresponding phasor plots in G–S space, including temperature-dependent trajectories.

The spectral phasor transformation is defined as:

$$
G = \frac{\int_{\lambda_i}^{\lambda_f} I(\lambda)
\cos\left(\frac{2\pi n (\lambda - \lambda_i)}{\Delta\lambda}\right)\, d\lambda}
{\int_{\lambda_i}^{\lambda_f} I(\lambda)\, d\lambda}
$$

$$
S = \frac{\int_{\lambda_i}^{\lambda_f} I(\lambda)
\sin\left(\frac{2\pi n (\lambda - \lambda_i)}{\Delta\lambda}\right)\, d\lambda}
{\int_{\lambda_i}^{\lambda_f} I(\lambda)\, d\lambda}
$$

where $I(\lambda)$ denotes the fluorescence intensity at wavelength $\lambda$, $n$ the harmonic order (typically $n = 1$), $\lambda_i$ the initial wavelength of the detection range, $\lambda_f$ the final wavelength, and $\Delta\lambda = \lambda_f - \lambda_i$ the spectral width.

Input data must consist of CSV files containing a temperature column and wavelength-resolved fluorescence intensity values encoded as column headers. The script automatically detects valid temperature fields and wavelength columns.

**Features**
- Import fluorescence emission spectra from CSV files
- Automatic detection of temperature and wavelength columns
- Spectral normalization
- Selection of wavelength range (400–550 nm)
- Calculation of spectral phasor coordinates (G, S)
- Interpolation of temperature series for improved comparability
- Export of processed spectra and phasor coordinates as CSV files
- Generation of spectral and phasor trajectory plots

**Input format**
Input data must be provided as CSV files containing:
- one temperature column
- wavelength-resolved fluorescence intensity values as column headers

The script automatically identifies valid temperature fields and wavelength columns.

**Output**
- Normalized spectral CSV files
- Phasor coordinate CSV files
- Emission spectrum plots
- Phasor plots in G–S space

