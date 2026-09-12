# fiberRemodelling

MATLAB code for analysing cell-associated remodelling of fibrous hydrogel networks from 3D fluorescence microscopy images.

This repository contains the matrix-remodelling analysis code associated with:

**Yuan et al. (2023). Synthetic fibrous hydrogels as a platform to decipher cell–matrix mechanical interactions.**  
*Proceedings of the National Academy of Sciences*, 120(15), e2216934120.  
https://doi.org/10.1073/pnas.2216934120

## Overview

The workflow uses separate fluorescence channels for cells and the polymer network to:

- Segment cells in 3D.
- Identify the densified polymer network surrounding cells.
- Generate 3D visualisations of cells and polymer.
- Calculate segmentation statistics and distances between cell and polymer boundaries.
- Save analysis results for subsequent plotting.

The repository also includes scripts for comparing experimental conditions and analysing culture and relaxation experiments.

## Requirements

- MATLAB.
- Image Processing Toolbox.
- Bio-Formats for Leica `.lif` files. A copy is included in `bfmatlab/bfmatlab`.
- A MATLAB-compatible C compiler if recompiling the optional 3D surface-smoothing functions.

Precompiled `.mexw64` files are provided for Windows. Other platforms may require recompilation.

## Setup

Download or clone the repository:

```bash
git clone https://github.com/BorisLouis/fiberRemodelling.git
```

Open MATLAB and navigate to the repository folder. Add the repository root and Bio-Formats folder to the MATLAB path:

```matlab
addpath(pwd);
addpath(fullfile(pwd, 'bfmatlab', 'bfmatlab'));
```

Keep the MATLAB package folders, such as `+Core` and `+Load`, in their original locations.

## Input data

The loading code contains support for Leica `.lif` files and split `.tif` images.

The analysis requires:

- A cell fluorescence channel.
- A polymer fluorescence channel.
- The lateral pixel size and axial slice spacing.

Configure the channel assignments to match the acquisition order:

```matlab
chan.c001 = 'cell';
chan.c002 = 'polymer';
```

For split TIFF images, the loader expects channel identifiers such as `c001` and three-digit slice/time identifiers such as `z001` and `t001` in the filenames.

Use one dataset per input folder. The loader searches for `.lif` files before `.tif` files.

## Single-dataset analysis

Open `main.m` and edit the user-input section:

```matlab
file.path = 'PATH_TO_YOUR_DATASET';
file.ext = '';

info.pxSizeXY = 570;
info.pxSizeZ = 570;

chan.c001 = 'cell';
chan.c002 = 'polymer';
```

Replace the example spatial calibration with the values for your images. Use nanometres consistently; check plot labels before interpreting or exporting physical measurements.

Run the script section by section to inspect the loaded channels and segmentation results before continuing to quantitative analysis.

The main analysis sequence is:

1. Load the microscopy data.
2. Segment the cell volume.
3. Extract the densified polymer network.
4. Visualise the segmented structures.
5. Calculate statistics and cell–polymer distances.
6. Save the results.

Segmentation thresholds and filtering parameters are defined in `+Core/fiberRemodelling.m` and may need adjustment for different imaging conditions.

## Batch analysis

`mainFolder.m` provides a batch-processing workflow. It searches the selected parent directory for subfolders containing a `Split` directory:

```text
ParentDirectory/
├── Sample01/
│   └── Split/
└── Sample02/
    └── Split/
```

Update the parent path, spatial calibration and channel assignments before use.

## Outputs

Depending on the selected workflow, outputs include:

- `channels.mat`: cached image channels.
- `results.mat`: analysis results stored in the variable `res`.
- `Figures/`: generated MATLAB figures and PNG images.
- Additional rendering files produced by the batch script.

Existing `channels.mat` files are reused automatically. If you change the input images or channel assignments, move or rename the cached file before reloading.

## Current limitations

These are research scripts with dataset-specific paths and parameters.

The TIFF-loading path currently references `chan` without receiving it as an argument and needs correction before use. The optional intensity-distribution routine also requires review for compatibility with the current cell-mask storage format.

Inspect segmentation results for each dataset. A tested MATLAB version and example dataset are not specified in this repository.

## Citation

If you use this code, please cite both the repository and the associated paper:

> Yuan, H., et al. (2023). Synthetic fibrous hydrogels as a platform to decipher cell–matrix mechanical interactions. *PNAS*, 120(15), e2216934120. https://doi.org/10.1073/pnas.2216934120

See [CITATION.cff](CITATION.cff) for citation metadata.

<!-- Add the repository's Zenodo DOI here when available. -->

## License

The project is distributed under the [MIT License](LICENSE).

Bundled third-party components retain their respective licenses and notices.
