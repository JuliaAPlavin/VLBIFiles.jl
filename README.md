# 🌌 VLBIFiles.jl

A comprehensive Julia package for reading and writing various data formats commonly used in **Very Long Baseline Interferometry (VLBI)**.

## 🔍 Overview

VLBIFiles.jl provides a unified interface for handling diverse file formats used in VLBI astronomy. The package features automatic format detection and offers flexible data representations tailored to different analysis workflows.

## 📊 Supported Data Types

### 📷 Images
Astronomical images represented as `KeyedArray`s from [AxisKeys.jl](https://github.com/mcabbott/AxisKeys.jl):
- **FITS files**: ✅ Read, ✅ Write

### 📡 Visibilities
Complex interferometric visibility data with multiple representation options:
- **UVFITS files**: ✅ Read, ❌ Write
- **FITS-IDI files**: 🟡 Read (partial support), ❌ Write

### 🔭 Source Models
Parametric astronomical source models powered by [InterferometricModels.jl](https://github.com/JuliaAPlavin/InterferometricModels.jl):
- **Difmap model files**: ✅ Read, ✅ Write
- **CLEAN components** (embedded in FITS images): ✅ Read, ❌ Write

## 🚀 Quick Start

```julia
using VLBIFiles

# 🎯 Automatic format detection and loading
data = VLBI.load("your_data_file")

# 🎛️ Specify output data type explicitly
image = VLBI.load(KeyedArray, "image.fits")           # Load as image
uvdata = VLBI.load(VLBI.UVData, "visibilities.uvf")   # Load as UV data structure
model = VLBI.load(MultiComponentModel, "source.mod")  # Load as source model
clean_model = VLBI.load(MultiComponentModel, "clean_image.fits")  # Extract CLEAN model

# 📋 Convert visibilities to table format for analysis
uv_table = VLBI.load("visibilities.uvfits") |> uvtable
```

## 📚 Documentation & Examples

For more comprehensive documentation and usage examples, explore the [Pluto notebook](https://aplavin.github.io/VLBIFiles.jl/notebooks/examples.html).

> [!NOTE]
> **📜 Package Evolution**: This package was originally named `VLBIData.jl`. In 2025, [VLBIData.jl](https://github.com/JuliaAPlavin/VLBIData.jl) was refocused exclusively on data structure definitions and visibility calculations, while its file I/O functionality was moved to `VLBIFiles.jl`. For continuity, `VLBIFiles.jl` continues the `VLBIData` versioning sequence starting from v0.3.32.
