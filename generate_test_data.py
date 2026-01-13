import pandas as pd
import numpy as np

# Create a dictionary of test data
data = {
    # --- ID / Metadata ---
    'Date': pd.date_range(start='2024-01-01', periods=10, freq='D'),
    'Location_ID': range(101, 111),

    # --- VEGETATION INDICES VARIABLES (Reflectance 0.0 to 1.0) ---
    'NIR':  [0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.72, 0.75, 0.78, 0.80],
    'Red':  [0.05, 0.06, 0.07, 0.08, 0.04, 0.05, 0.03, 0.04, 0.05, 0.06],
    'Green':[0.08, 0.09, 0.10, 0.11, 0.09, 0.08, 0.07, 0.08, 0.09, 0.10],
    'Blue': [0.02, 0.03, 0.02, 0.04, 0.03, 0.02, 0.01, 0.02, 0.03, 0.02],
    'R531': [0.12, 0.13, 0.14, 0.15, 0.13, 0.12, 0.11, 0.12, 0.13, 0.14],
    'R570': [0.10, 0.11, 0.12, 0.13, 0.11, 0.10, 0.09, 0.10, 0.11, 0.12],
    'R700': [0.15, 0.16, 0.18, 0.20, 0.17, 0.16, 0.15, 0.16, 0.17, 0.18],
    'R670': [0.05, 0.06, 0.07, 0.08, 0.05, 0.04, 0.03, 0.04, 0.05, 0.06], 
    'R550': [0.09, 0.10, 0.11, 0.12, 0.10, 0.09, 0.08, 0.09, 0.10, 0.11],
    'R720': [0.25, 0.28, 0.30, 0.35, 0.27, 0.25, 0.24, 0.26, 0.28, 0.30],
    'R750': [0.40, 0.45, 0.50, 0.55, 0.42, 0.40, 0.38, 0.41, 0.44, 0.46],
    'R790': [0.48, 0.52, 0.58, 0.62, 0.50, 0.48, 0.46, 0.49, 0.53, 0.55],
    'R800': [0.50, 0.55, 0.60, 0.65, 0.52, 0.50, 0.48, 0.51, 0.55, 0.58],
    'R753': [0.42, 0.47, 0.52, 0.57, 0.44, 0.42, 0.40, 0.43, 0.46, 0.48],
    'R708': [0.20, 0.22, 0.25, 0.28, 0.21, 0.20, 0.19, 0.21, 0.23, 0.24],
    'R681': [0.08, 0.09, 0.10, 0.11, 0.09, 0.08, 0.07, 0.08, 0.09, 0.10],

    # --- CONSTANTS FOR VEGETATION INDICES ---
    'L': [0.5] * 10,       # Soil adjustment factor
    'a': [1.2] * 10,       # Slope of soil line
    'b': [0.02] * 10,      # Intercept of soil line
    'NDRE': [0.35, 0.38, 0.40, 0.42, 0.36, 0.35, 0.34, 0.37, 0.39, 0.41], # Pre-calc for CCCI
    'NDRE_min': [0.1] * 10,
    'NDRE_max': [0.9] * 10,

    # --- TEMPERATURE DATA (Celsius) ---
    'T_mean': [25, 26, 24, 28, 30, 32, 29, 27, 26, 25],
    'T_max':  [30, 31, 29, 34, 36, 38, 35, 32, 31, 30],  # Used as TM
    'T_min':  [20, 21, 19, 22, 24, 26, 23, 21, 20, 19],  # Used as Tm
    'T_current': [5, 8, 14, 17, 21, 25, 12, 7, 2, 18],   # Hourly sample for Chill Units
    'Ta':     [25, 26, 24, 28, 30, 32, 29, 27, 26, 25],  # Mean Air Temp (alias)
    'Td':     [15, 16, 14, 18, 20, 22, 19, 17, 16, 15],  # Dew Point

    # --- GDD THRESHOLDS ---
    'TM': [30, 31, 29, 34, 36, 38, 35, 32, 31, 30], # Duplicate of T_max
    'Tm': [20, 21, 19, 22, 24, 26, 23, 21, 20, 19], # Duplicate of T_min
    'Tb': [10] * 10,   # Base Temp (Min threshold)
    'TB': [35] * 10,   # Ceiling Temp (Max threshold)

    # --- EVAPOTRANSPIRATION VARIABLES ---
    'Ra':    [35, 36, 34, 38, 40, 41, 39, 37, 36, 35], # Extraterrestrial Radiation
    'Rn':    [15, 16, 14, 18, 20, 21, 19, 17, 16, 15], # Net Radiation
    'G':     [0.5, 0.6, 0.4, 0.8, 1.0, 1.2, 0.9, 0.7, 0.6, 0.5], # Soil Heat Flux
    'u2':    [2.0, 2.5, 1.8, 3.0, 3.5, 4.0, 3.2, 2.8, 2.5, 2.2], # Wind Speed
    'es_ea': [1.5, 1.6, 1.4, 1.8, 2.0, 2.2, 1.9, 1.7, 1.6, 1.5], # Vapour Pressure Deficit
    'Delta': [0.15, 0.16, 0.14, 0.18, 0.20, 0.22, 0.19, 0.17, 0.16, 0.15], # Slope Vapour Curve
    'Gamma': [0.06] * 10, # Psychrometric Constant
    'Kc':    [0.85] * 10, # Crop Coefficient
    'Rl':    [500, 520, 480, 550, 600, 620, 580, 540, 520, 500], # Solar Rad (Langleys/day)
    'z':     [100] * 10, # Elevation
    'lat':   [30] * 10   # Latitude
}

df = pd.DataFrame(data)

# Save to Excel
df.to_excel("agricultural_test_data.xlsx", index=False)
print("File 'agricultural_test_data.xlsx' created successfully!")
