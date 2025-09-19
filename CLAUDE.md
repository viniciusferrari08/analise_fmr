# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a scientific Python codebase for analyzing Ferromagnetic Resonance (FMR) experimental data from thin films. The software implements the Kittel equation to extract magnetic parameters (Ms, Ha, Hk) from frequency vs. resonance field measurements.

## Core Architecture

### Main Analysis Pipeline
- `analise_fmr_completa.py` - Complete FMR analysis pipeline that orchestrates the entire workflow
- `extrair_campo_ressonancia.py` - FMRSpectrumAnalyzer class for extracting resonance fields from experimental spectra
- `fmr_fitting.py` - FMRFitting class implementing Kittel equation fitting and parameter extraction

### Key Classes and Methods
- **FMRSpectrumAnalyzer**: Processes .dat files from experimental measurements
  - `load_spectrum()` - Loads individual spectrum files
  - `find_resonance_field_derivative()` - Extracts Hr using derivative peak detection methods
  - `_fit_derivative_lorentzian()` - Fits Lorentzian derivative to extract Hr and linewidth ΔH
  - `process_measurement_folder()` - Batch processes all spectra in a folder
  - `plot_spectrum()` - Visualizes individual spectra with resonance field and linewidth info
  - `plot_all_spectra()` - Creates subplot grid showing multiple spectra with ΔH values
- **FMRFitting**: Implements magnetic parameter fitting
  - `kittel_equation()` - Core Kittel equation implementation
  - `fit_perpendicular_data()` and `fit_inplane_data()` - Geometry-specific fitting methods
  - `calculate_errors()` - Statistical error analysis from covariance matrix

## Development Commands

### Installation
```bash
pip install -r requirements.txt
```

### Running Analysis
```bash
# Complete analysis of experimental data
python analise_fmr_completa.py

# Individual spectrum analysis
python extrair_campo_ressonancia.py

# Fitting only (if resonance data exists)
python fmr_fitting.py
```

### Dependencies
- numpy >= 1.20.0 (numerical computations)
- matplotlib >= 3.5.0 (plotting and visualization)
- scipy >= 1.7.0 (curve fitting and optimization)
- pandas >= 1.3.0 (data handling)

## Data Structure

### Experimental Data Location
- `EAFExp2025/` - Root directory for experimental measurements
- `EAFExp2025/NiFe_Cu_6nm/medidas selec/` - Selected NiFe/Cu measurements
- Individual `.dat` files contain: [frequency, current, magnetic_field, fmr_signal]

### Output Files
- `campos_ressonancia_NiFe_Cu.txt` - Extracted resonance fields and linewidths
- `resultado_analise_completa.txt` - Complete analysis report
- `ajuste_fmr_final.png` - Final fitting visualization
- `teste_lorentzian_fit.png` - Example of Lorentzian derivative fitting (when testing)

## Physical Theory Implementation

### Kittel Equation Forms
**Perpendicular geometry (θ = 90°):**
```
f = (γ/2π) × √[(Hr + Ha + Hk)(Hr + Ha + Hk + Ms)]
```

**In-plane geometry (θ = 0°):**
```
f = (γ/2π) × √[(Hr + Ha)(Hr + Ha + Ms - Hk)]
```

Where:
- γ = 2.8×10¹⁰ Hz/T (gyromagnetic ratio)
- Ms = Saturation magnetization (T)
- Ha = In-plane anisotropy field (T)  
- Hk = Perpendicular anisotropy field (T)

## Code Conventions

### File Organization
- Main analysis scripts in root directory
- Reference codes in `códigos FMR/` (legacy implementations)
- Experimental data in `EAFExp2025/` hierarchy

### Data Processing Workflow
1. Load experimental spectra (.dat files) with 4-column format
2. Extract resonance fields (Hr) and linewidths (ΔH) using Lorentzian derivative fitting
3. Apply quality control with R² threshold (>0.7) and fallback methods
4. Store results with both field and linewidth information
5. Fit frequency vs. field data to Kittel equation for magnetic parameters
6. Calculate statistical errors and generate comprehensive visualizations

### Unit Conventions
- Frequencies: GHz (input) → Hz (calculations)
- Magnetic fields: Oe (experimental) → T (calculations) → mT (results)
- Linewidths: Oe (fitted) → mT (displayed)
- Use scipy.optimize.curve_fit for non-linear Lorentzian derivative fitting with corrected mathematical form
- Error propagation via covariance matrix analysis

## Important Notes

- The codebase assumes derivative absorption spectra (dχ"/dH format)
- Experimental files must have 4-column structure: [freq, current, field, signal]
- **New Feature**: Automatic extraction of FMR linewidths (ΔH) from Lorentzian derivative fitting
- Default method now uses `method="fit"` for improved accuracy over zero-crossing
- Quality control: R² > 0.7 threshold with automatic fallback to zero-crossing method
- Physical parameter validation checks for reasonable magnetic values
- Comprehensive visualization includes both Hr and ΔH information in plots and tables

## Usage Examples

### Basic Resonance Field and Linewidth Extraction
```python
from extrair_campo_ressonancia import FMRSpectrumAnalyzer

analyzer = FMRSpectrumAnalyzer()
results = analyzer.process_measurement_folder("EAFExp2025/NiFe_Cu_6nm/medidas selec")

# Results include both Hr and ΔH
for freq, Hr in results.items():
    Delta_H = analyzer.spectra_data[freq]['Delta_H']
    print(f"{freq} GHz: Hr = {Hr:.1f} Oe, ΔH = {Delta_H:.1f} Oe")
```

### Visualization with Linewidth Information
```python
# Plot individual spectrum with fitting overlay
analyzer.plot_spectrum(5.0, show_resonance=True, show_fit=True)

# Plot all spectra with ΔH in subplot titles
analyzer.plot_all_spectra()
```

### Expected Output Format
```
Freq(GHz)    Hr(Oe)    Hr(mT)    ΔH(Oe)    ΔH(mT)
5.0          355.8     35.6      17.8      1.8
6.0          505.8     50.6      21.1      2.1
...
```