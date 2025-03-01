module VLBIFiles

using Reexport
using DataManipulation
using Tables
import Tables: table
using FITSIO: FITSHeader, FITS, TableHDU, read_header, colnames, FITSIO
using Dates
using Unitful, UnitfulAstro, UnitfulAngles
using AxisKeys
using StaticArrays
using StructArrays
using AccessorsExtra
using DelimitedFiles: readdlm
using DateFormats: julian_day
using PyFormattedStrings
using Statistics
using Uncertain
@reexport using InterferometricModels
@reexport using VLBIData

export VLBI, table, uvtable

include("grouphdu.jl")

include("fitsutils.jl")
include("uvdata.jl")
include("fitsimage.jl")
include("difmap_files.jl")
include("loading.jl")

baremodule VLBI
import ..VLBIData
Core.eval(VLBI, Expr(:import, Expr(:(:), Expr(:., :., :., :VLBIData),
          [Expr(:., n) for n in VLBIData.names(VLBIData.VLBI; imported=true) if Core.:(===)(Core.:(===)(n, :VLBI), false)]...)))

import ..VLBIFiles:
    VLBIFiles,
    load, save, guess_type,
    table, uvtable, read_data_raw, read_data_arrays,
    FitsImage, FrequencyWindow, UVHeader, AntArray, UVData,
    pixel_size, pixel_steps, pixel_area
using ..InterferometricModels
end

using PrecompileTools
@compile_workload begin
    load(joinpath(@__DIR__, "../test/data/map.fits"))
    load(joinpath(@__DIR__, "../test/data/vis.fits"))
    load(joinpath(@__DIR__, "../test/data/difmap_model_empty.mod"))
end

end
