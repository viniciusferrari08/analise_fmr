# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a scientific Python codebase for analyzing Ferromagnetic Resonance (FMR) experimental data from thin films. The software implements the Kittel equation to extract magnetic parameters (Ms, H_eff) from frequency vs. resonance field measurements.

## Core Architecture

### Main Analysis Pipeline
- `fmr_gui.py` - **PRIMARY TOOL**: Interactive GUI for complete FMR analysis workflow
- `analise_fmr_completa.py` - Legacy command-line analysis pipeline
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
# GUI Application (RECOMMENDED)
python fmr_gui.py

# Legacy command-line analysis
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

**Simplified in-plane geometry (current implementation):**
```
ω₀ = γ√((Hr - H_eff)(Hr - H_eff + Ms))
```

Where:
- γ = 2.8×10¹⁰ Hz/T (gyromagnetic ratio)
- Ms = Saturation magnetization (T)
- H_eff = Effective field including anisotropies (T)
- Hr = Resonance field (T)

This simplified model provides stable fitting with only 2 parameters, reducing error propagation.

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

## GUI Application Features

### FMR GUI (`fmr_gui.py`)
The graphical interface provides a complete workflow for FMR analysis:

**Key Features:**
- **Automatic Processing**: Uploads and processes spectra automatically
- **Multiple File Selection**: Select and remove multiple files (Ctrl/Shift + click)
- **Real-time Visualization**: View individual spectra with Lorentzian fits
- **Automatic Kittel Fitting**: Fits magnetic parameters when ≥3 spectra loaded
- **Linewidth Analysis**: Automatic extraction and plotting of FMR linewidths
- **Export Options**: Save plots in PNG, PDF, or SVG formats

**Workflow:**
1. Upload `.dat` files → automatic processing extracts Hr and ΔH
2. View individual spectra with resonance field and fit overlay
3. Kittel fitting runs automatically (displays Ms and H_eff)
4. Analyze linewidth vs frequency with linear fit
5. Export all plots with proper formatting

**Interface Tabs:**
- **Espectro Individual**: View selected spectrum with fit
- **Largura de Linha**: ΔH vs frequency analysis
- **Ajuste de Kittel**: Magnetic parameter fitting results

## Usage Examples

### GUI Usage (Recommended)
```bash
python fmr_gui.py
# 1. Click "Upload Arquivos .dat"
# 2. Select experimental files
# 3. View results in tabs automatically
# 4. Export plots as needed
```

### Command-line Usage (Legacy)
```python
from extrair_campo_ressonancia import FMRSpectrumAnalyzer

analyzer = FMRSpectrumAnalyzer()
results = analyzer.process_measurement_folder("EAFExp2025/NiFe_Cu_6nm/medidas selec")

# Results include both Hr and ΔH
for freq, Hr in results.items():
    Delta_H = analyzer.spectra_data[freq]['Delta_H']
    print(f"{freq} GHz: Hr = {Hr:.1f} Oe, ΔH = {Delta_H:.1f} Oe")
```