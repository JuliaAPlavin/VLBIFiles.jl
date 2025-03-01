# VLBIFiles.jl

`VLBIFiles.jl` provides functions for reading and writing various data formats typically used in very long baseline interferometry (VLBI).

- Images, represented as `KeyedArray`s (from [AxisKeys.jl](https://github.com/mcabbott/AxisKeys.jl))
  - `fits` files: ✅ read, ✅ write
- Visibilities, represented as (almost raw) `KeyedArray`s or as flat tables
  - `uvfits` files: ✅ read, ❌ write
  - `fits-idi` files: 🟡 read, ❌ write
- Source models, represented by [InterferometricModels.jl](https://github.com/JuliaAPlavin/InterferometricModels.jl)
  - `difmap` model files: ✅ read, ✅ write
  - `clean` components within `fits` images: ✅ read, ❌ write

To read a file, use:
```julia
using VLBIFiles

# determine the file type automatically and read into its default representation
VLBI.load("filename")

# specify the resulting type manually:
VLBI.load(KeyedArray, "filename.fits")  # image as a KeyedArray
VLBI.load(VLBI.UVData, "filename.uvf")  # visibilities as a UVData object
VLBI.load(MultiComponentModel, "filename.mod")  # model from a difmap model file
VLBI.load(MultiComponentModel, "filename.fits")  # CLEAN model from a fits image

# read visibilities and return them as a "uv table"
# see VLBIData.jl for the data structure definition
VLBI.load("filename") |> uvtable
```

See the [notebook](https://aplavin.github.io/VLBIFiles.jl/notebooks/examples.html) for docs and more usage examples.

> [!NOTE]
> This package was originally named `VLBIData.jl`. Then, in 2025, [VLBIData.jl](https://github.com/JuliaAPlavin/VLBIData.jl) was refocused solely on data structure definitions and visibility calculations, and its file I/O functions were transferred to `VLBIFiles.jl` (this package). To maintain continuity, `VLBIFiles.jl` carries on the `VLBIData` versioning starting at v0.3.31.
