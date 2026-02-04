# VLBIFiles.jl

Reads and writes file formats used in Very Long Baseline Interferometry (VLBI): images, visibilities, and source models. Automatic format detection.

```julia
using VLBIFiles

# Automatic format detection
data = VLBI.load("your_data_file")

# Specify output type explicitly
image = VLBI.load(KeyedArray, "image.fits")
uvdata = VLBI.load(VLBI.UVData, "visibilities.uvf")
model = VLBI.load(MultiComponentModel, "source.mod")
clean_model = VLBI.load(MultiComponentModel, "clean_image.fits")

# Convert visibilities to table format for analysis
uv_table = VLBI.load("visibilities.uvfits") |> uvtable
```

Part of the [VLBIData](https://github.com/JuliaAPlavin/VLBIData.jl) ecosystem — loaded visibility data is directly compatible with VLBIData.jl operations (averaging, closures, error rescaling, polarization, and more).

## Supported Data Types

### Images
Astronomical images as `KeyedArray`s from [AxisKeys.jl](https://github.com/mcabbott/AxisKeys.jl):
- **FITS files**: ✅ Read, ✅ Write

### Visibilities
Complex interferometric visibility data, read-only:
- **UVFITS files**: ✅ Read
- **FITS-IDI files**: 🟡 Read (partial)
- **HOPS alist files**: ✅ Read

### Source Models
Parametric source models via [InterferometricModels.jl](https://github.com/JuliaAPlavin/InterferometricModels.jl):
- **Difmap model files**: ✅ Read, ✅ Write
- **CLEAN components** (from FITS images): ✅ Read, ❌ Write

## Documentation

More examples in the [Pluto notebook](https://aplavin.github.io/VLBIFiles.jl/notebooks/examples.html). See the [VLBIData.jl README](https://github.com/JuliaAPlavin/VLBIData.jl) for documentation on working with loaded data: visibility tables, averaging, closures, polarization, and more.
