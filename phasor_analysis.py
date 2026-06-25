import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# Harmonic used for phasor calculation
harmonic = 1

# Fixed wavelength range used for phasor analysis
final_min = 400
final_max = 550

# Folder containing input CSV files
input_folder = "data"

# Output folder
output_folder = os.path.join(
    input_folder,
    "output"
)

# Create output directory if it does not exist
os.makedirs(output_folder, exist_ok=True)

# Find all CSV files inside input folder
csv_files = glob.glob(os.path.join(input_folder, "*.csv"))

print(f"Found {len(csv_files)} CSV files")


def calculate_phasor(intensities, wavelengths):
  
    #Convert input to numpy arrays
    intensities = np.array(intensities, dtype=float)
    wavelengths = np.array(wavelengths, dtype=float)

    # Remove negative intensity values
    intensities[intensities < 0] = 0

    total_intensity = np.sum(intensities)

    # Avoid division by zero
    if total_intensity == 0:
        return np.nan, np.nan

    # Select wavelength range   
    lambda_min = wavelengths.min()
    lambda_max = wavelengths.max()
    
    phi = (
        2
        * np.pi
        * harmonic
        * (wavelengths - lambda_min)
        / (lambda_max - lambda_min)
    )

    # Compute real (G) and imaginary (S) phasor components
    G = np.sum(
        intensities * np.cos(phi)
    ) / total_intensity

    S = np.sum(
        intensities * np.sin(phi)
    ) / total_intensity

    return G, S


for input_file in csv_files:

    #File name
    filename_no_ext = os.path.splitext(
        os.path.basename(input_file)
    )[0]

    data = pd.read_csv(input_file)

    print("\nDetected columns:")
    print(data.columns)

    # Detect temperature column
    possible_temp_cols = [
        "temperature_c",
        "temperature",
        "Temperature"
    ]

    temp_col = None

    for col in possible_temp_cols:

        if col in data.columns:

            temp_col = col
            break

    if temp_col is None:

        raise ValueError(
            "No temperature column found."
        )

    temps = data[temp_col].values

    # Detect wavelength column
    wavelength_cols = []

    for col in data.columns:

        try:

            float(col)
            wavelength_cols.append(col)

        except ValueError:

            pass

    wavelengths = np.array(
        wavelength_cols,
        dtype=float
    )

    # Sort wavelengths in ascending order
    sort_idx = np.argsort(wavelengths)

    wavelengths = wavelengths[sort_idx]

    wavelength_cols = np.array(
        wavelength_cols
    )[sort_idx]

    
    spectra = []

    for _, row in data.iterrows():

        intensities = row[
            wavelength_cols
        ].values.astype(float)

        spectra.append(intensities)

    spectra = np.array(spectra)

    # Normalize spectra to maximum intensity value 
    spectra_norm = []

    for spectrum in spectra:

        max_val = np.max(spectrum)

        # Avoid division by zero
        if max_val == 0:

            spectra_norm.append(spectrum)

        else:

            spectra_norm.append(
                spectrum / max_val
            )

    spectra_norm = np.array(spectra_norm)

    # Plot normalized spectra
    fig, ax = plt.subplots(figsize=(8, 6))

    cmap = plt.cm.plasma

    norm = plt.Normalize(
        temps.min(),
        temps.max()
    )

    for i, spectrum in enumerate(spectra_norm):

        ax.plot(
            wavelengths,
            spectrum,
            color=cmap(norm(temps[i]))
        )

    sm = plt.cm.ScalarMappable(
        cmap=cmap,
        norm=norm
    )

    fig.colorbar(
        sm,
        ax=ax,
        label="Temperature (deg C)"
    )

    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Normalized Intensity")

    ax.set_title(
        f"Normalized Spectra\n{filename_no_ext}"
    )

    ax.grid(True)

    plt.show()

    
    print(
        f"\nUsing wavelength range: "
        f"{final_min:.1f} - {final_max:.1f} nm"
    )

    valid = (
        (wavelengths >= final_min)
        &
        (wavelengths <= final_max)
    )

    wavelengths_final = wavelengths[valid]

    spectra_final = spectra_norm[:, valid]

    # Calculate phasor coordinates
    G_list = []
    S_list = []

    for intensities in spectra_final:

        G, S = calculate_phasor(
            intensities,
            wavelengths_final
        )

        G_list.append(G)
        S_list.append(S)

    G = np.array(G_list)
    S = np.array(S_list)

    # Remove invalid values
    valid_mask = (
        ~np.isnan(G)
        &
        ~np.isnan(S)
    )

    G = G[valid_mask]
    S = S[valid_mask]

    temps = temps[valid_mask]
    spectra_final = spectra_final[valid_mask]

   
    # Sort Data by temperature
    sort_idx = np.argsort(temps)

    temps = temps[sort_idx]

    G = G[sort_idx]
    S = S[sort_idx]

    spectra_final = spectra_final[sort_idx]

    # Remove duplicate temperatures before interpolation
    temps_unique, unique_idx = np.unique(
        temps,
        return_index=True
    )

    temps = temps_unique

    G = G[unique_idx]
    S = S[unique_idx]

    spectra_final = spectra_final[unique_idx]

    # Interpolate G and S
    new_temps = np.arange(
        temps.min(),
        temps.max(),
        0.5
    )

    interp_G = interp1d(temps, G)

    interp_S = interp1d(temps, S)

    G_interp = interp_G(new_temps)

    S_interp = interp_S(new_temps)

    # Interpolate normalized spectra 
    spectra_interp = []

    for i in range(spectra_final.shape[1]):

        interp_spec = interp1d(
            temps,
            spectra_final[:, i],
            bounds_error=False,
            fill_value="extrapolate"
        )

        spectra_interp.append(
            interp_spec(new_temps)
        )

    spectra_interp = np.array(
        spectra_interp
    ).T

    # Export normalized spectra
    spectra_export = pd.DataFrame({

        "Wavelength_nm": wavelengths_final
    })

    for i, T in enumerate(new_temps):

        spectra_export[
            f"T_{T:.1f}C"
        ] = spectra_interp[i]

    spectra_output_file = os.path.join(
        output_folder,
        f"NormalizedSpectra_{filename_no_ext}.csv"
    )

    spectra_export.to_csv(
        spectra_output_file,
        index=False
    )


    # Plot phasor trajectories
    fig, ax = plt.subplots(figsize=(7, 7))

    ax.plot(
        G_interp,
        S_interp,
        marker="o",
        color="red"
    )

    # Plot unit circle
    theta = np.linspace(
        0,
        2 * np.pi,
        300
    )

    ax.plot(
        np.cos(theta),
        np.sin(theta),
        "k--"
    )

    ax.set_xlabel("G")
    ax.set_ylabel("S")

    ax.set_title(
        f"Phasor Plot\n{filename_no_ext}"
    )

    ax.axis("equal")

    ax.grid(True)

    plt.show()

    # Export G and S values
    output_file = os.path.join(
        output_folder,
        f"Phasor_{filename_no_ext}.csv"
    )

    results = pd.DataFrame({

        "Temperature": new_temps,
        "G": G_interp,
        "S": S_interp
    })

    results.to_csv(
        output_file,
        index=False
    )
