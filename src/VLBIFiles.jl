module VLBIFiles

using Reexport
using DataManipulation
using Tables
import Tables: table, columnnames
using FITSIO: FITSHeader, FITS, TableHDU, read_header, colnames, FITSIO
using Dates
using Unitful, UnitfulAstro, UnitfulAngles
using AxisKeys
using AxisKeys.IntervalSets
using StaticArrays
using StructArrays
using AccessorsExtra
using DelimitedFiles: readdlm
using DateFormats: julian_day
using PyFormattedStrings
using Statistics
using Uncertain
@reexport using InterferometricModels
using InterferometricModels: components
@reexport using VLBIData
import VLBIData: frequency, uvtable
using Dictionaries

export VLBI, table, uvtable

include("fitsio/grouphdu.jl")
include("fitsio/fitsutils.jl")
include("fitsio/lazy_table.jl")

include("uvdata.jl")
include("fitsimage.jl")
include("difmap_files.jl")
include("loading.jl")

baremodule VLBI
using Reexport

import ..VLBIData
_names = [n for n in VLBIData.names(VLBIData.VLBI) if Core.:(===)(Core.:(===)(n, :VLBI), false)]
Core.eval(VLBI, Expr(:import, Expr(:(:), Expr(:., :., :., :VLBIData), [Expr(:., n) for n in _names]...)))
Core.eval(VLBI, Expr(:export, _names...))

@reexport import ..VLBIFiles:
    VLBIFiles,
    load, save, guess_type,
    table, read_data_raw, read_data_arrays,
    FitsImage, FrequencyWindow, UVHeader, AntArray, UVData,
    pixel_size, pixel_steps, pixel_area,
    image_stored, image_clean, image_residual, image_clean_with_residual
using ..InterferometricModels
end

# using PrecompileTools
# @compile_workload begin
#     load(joinpath(@__DIR__, "../test/data/map.fits"))
#     load(joinpath(@__DIR__, "../test/data/vis.fits"))
#     load(joinpath(@__DIR__, "../test/data/difmap_model_empty.mod"))
# end

end
