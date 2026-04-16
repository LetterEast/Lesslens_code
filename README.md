# Lensless Holographic Imaging Reconstruction

This repository provides a complete MATLAB pipeline for multi-distance and multi-angle in-line lensless holographic imaging reconstruction.

## Features
- **Multi-distance/Multi-angle Phase Retrieval**: Robust framework using APRW (MDPRF style) with Nesterov-like weighted feedback momentum.
- **Alias-Free Off-axis Angular Spectrum Propagation**: GPU-accelerated propagation capable of exact frequency demodulation to avoid aliasing artifacts caused by high-angle illuminations.
- **Physical Optimization Integration**: Includes physical bounds constraint mapping and CCTV (Complex Total Variation) capabilities to suppress speckle noise.
- **Data-Driven Image Alignment**: Subpixel shift registration for automatic Z-stack and multi-angle hologram alignment.

## Getting Started

1. Clone or download this repository.
2. Open MATLAB and navigate to the project directory.
3. Open `main.m` or `batch_run.m`.
4. Modify the `RawImgFolder` to point to your holographic data.
5. Set your physical parameters:
   - `WaveLength` (e.g., 514 nm)
   - `PixelSize` (sensor pixel pitch)
   - `DistanceIntervalSet` (propagation distances)
6. Run the script. The script performs reconstruction and will save out amplitude and phase profiles automatically inside a `ResultFolder`.

## Code Structure
- `/core`: Main iterative solver algorithms (`getRec_APRW`), angular spectrum operators (`propGPU`), and illumination modelling.
- `/utils`: Helpful subroutines for image loading, pre-processing, sub-pixel shift bounding (`dftregistration`), and object mapping.
- `/calibration`: Code for calibrating optical distance parameters.
- `main.m`: Standard execution entry point.
- `batch_run.m`: Loops over multiple folders to reconstruct large datasets automatically.

## Requirements
- MATLAB (Tested on recent versions)
- Parallel Computing Toolbox (for `gpuArray` acceleration)
- Image Processing Toolbox

## Acknowledgements
Contains implementations adapted from the CCTV phase retrieval constraints and anti-aliasing Angular Spectrum solutions.
