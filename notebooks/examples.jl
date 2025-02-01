### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 644d9d5f-9b1c-4dd2-b688-faac9e7c2473
using Pkg

# ╔═╡ 993027e8-7f5b-4952-8ed7-0345a9522d3a
pkg"status -m --outdated"

# ╔═╡ 4eb97406-7713-11ec-00ca-8be6bef77030
using VLBIFiles

# ╔═╡ ead10789-68f5-4931-b177-622030df739d
using CairoMakie

# ╔═╡ 0c44dd70-5a5c-46b9-b085-1a52c22078ca
using MakieExtra

# ╔═╡ 17a5e5ad-c3a7-41a0-bf56-410675e95475
using VLBIPlots

# ╔═╡ 7a4e96c8-8ab0-4338-82b5-e07f0afdaae5
using DataManipulation

# ╔═╡ 9c6efe8c-385e-4446-acdf-bd19cffe31e2
using DisplayAs: Text as AsText

# ╔═╡ be208a3e-9a03-4194-8be7-05547d9871a7
using RectiGrids

# ╔═╡ 4e7278e0-0d12-4762-8ace-e7cbe8eb371c
using StaticArrays

# ╔═╡ dc179b8b-a5b7-4d2d-bf86-afc0fc5a376b
using Unitful, UnitfulAngles, UnitfulAstro

# ╔═╡ f3898d40-517d-421b-b551-00d8ce1f64dd
using AxisKeys

# ╔═╡ 8f556fb7-c181-4e0f-a42d-1caa85f35f41
using Difmap

# ╔═╡ 398b65dd-67ae-4fa8-92d4-a73cf2992060
using IntervarrrSets

# ╔═╡ 412dc1f1-ea14-40ee-ad3b-28b89990738f
using Statistics

# ╔═╡ 7580f512-0653-479d-854b-f5aa78d85655
using PlutoUI

# ╔═╡ 5491a70e-8aea-42bc-b449-a56fd120d3fd
md"""
!!! info "VLBIFiles.jl"
	`VLBIFiles` reads and writes a range of data formats typically used in very long baseline interferometry (VLBI).

For convenience, it exports itself as `VLBI`: non-exported functions are accessed as `VLBI.load`.

See the table of contents and sections below for the list of supported formats.
The proper file format is guessed automatically, and typically one just needs the `VLBI.load(file)` method.

Some formats, such as "image FITS", may actually contain multiple datasets of different kinds (eg image + model). In such cases, `VLBI.load(file)` reads the main dataset (image) by default. Use `VLBI.load(type, file)` to read the dataset of `type` from the `file`: see examples below.
"""

# ╔═╡ 6cb0d385-5b13-4afe-bf8c-015938efc91d
md"""
`VLBIFiles` is well-tested. It is regularly confirmed to read all files in the [Astrogeo](http://astrogeo.org/vlbi_images/) database successfully, as well as datasets from other telescope networks such as the EHT.
"""

# ╔═╡ b685e8e4-efa5-453c-a27b-4453946a0839
md"""
# Image FITS
"""

# ╔═╡ 9af33a02-67f8-4233-93f5-4797ecf702d9
md"""
## Read image itself
"""

# ╔═╡ 3df0333d-734c-4d4b-bdec-49cac59f4681
md"""
Load the image data from a FITS file:
"""

# ╔═╡ 6bf4ae14-33d7-4faa-a664-80889ece70b7
fimg = VLBI.load("../test/data/map.fits");

# ╔═╡ 50c8dfe9-db54-4239-a14d-c8e047eab443
md"""
The returned object contains both the header information, and the image data as an array:
"""

# ╔═╡ 8b722ca5-3198-4c6c-94b1-4daf8fd2e718
fimg.header["OBJECT"]

# ╔═╡ 9abadd95-1aeb-436f-aab2-64dc2ad78830
KeyedArray(fimg) |> display

# ╔═╡ b104bdf4-cdee-4743-a4e7-31844f5338e7
md"""
As you see, the image data is a `KeyedArray`. It contains information about axes and their units for convenient look up, indexing, plotting, ...:
"""

# ╔═╡ dcd07de0-148d-4641-9cfc-79926e1ba3b1
KeyedArray(fimg)(dec=Near(5u"mas")) |> display

# ╔═╡ 203e7c90-11d3-4fc1-895a-ffcba4247b34
let
	fig, ax, plt = image(KeyedArray(fimg), colormap=:inferno, colorscale=SymLog(5e-4, vmin=0))
	Colorbar(fig[1,2], plt, label="Jy/beam")
	fig
end

# ╔═╡ 3075f803-d583-4764-a7d4-1bc14d8e441b
md"""
Access the so-called "image beam" - the effective point spread function - and add it to the image:
"""

# ╔═╡ b32c9639-c2c3-4a67-8554-57aefea9d196
beam(fimg)

# ╔═╡ 530aa386-84a7-47b2-ba76-a9c01fc5393d
let
	fig, ax, plt = image(KeyedArray(fimg)((-30..20)u"mas", (-25..25)u"mas"), colormap=:inferno, colorscale=SymLog(5e-4, vmin=0))
	beampoly!(ax, ustrip(beam(fimg)), position=(0.1, 0.1), color=:lightgreen)
	# plt.gca().add_artist(ScalebarArtist([(1, "mas")], loc="lower right", color=:w))
	fig
end

# ╔═╡ 6d0bf04e-51ee-43b6-8b7c-912fd051dfb5
md"""
## Read source model
"""

# ╔═╡ c4b58f30-e7db-459e-9d26-2d4bc57e5fbd
md"""
FITS image files commonly contain the source model itself, often as a set of delta functions - "CLEAN components". Components are stored in the `AIPS CC` FITS table.

To load this model instead of the image array, pass the corresponding type to `VLBI.load`:
"""

# ╔═╡ b4e14571-a3fa-41d2-98da-cff594f88202
fimg_mod = VLBI.load(MultiComponentModel, "../test/data/map.fits")

# ╔═╡ c8510200-f820-4b2a-8260-d21e8f7f7880
md"""
The `MultiComponentModel` type, as well as components themselves, are defined in the `InterferometricModels.jl` package. See its help for more details on how to use these models.
"""

# ╔═╡ ebcabc9f-0487-4cad-b7f7-fff448fcf203
md"""
For example, we can convolve the CLEAN model with the image beam:
"""

# ╔═╡ 7994f45f-6e01-4a62-ae3b-3031b46f2d7b
fimg_mod_convolved = convolve(fimg_mod, beam(fimg))

# ╔═╡ 9db795b2-2d79-45b6-a51c-bc09f4d05150
md"""
Convolution turned each delta component into an `EllipticGaussian` component. We can use this model to compute an image: the result should be similar to the image array stored in the FITS file, but without adding residual noise:
"""

# ╔═╡ 2f424454-1264-4887-927c-ad58d8631ecc
fimg_mod_img = intensity(fimg_mod_convolved).(
	grid(SVector,
		ra=range(30, -30, length=200)u"mas",
		dec=range(-30, 30, length=200)u"mas"
	)
);

# ╔═╡ 0a054aab-43d0-4caf-85c9-9303dec7c8fb
fimg_mod_img |> AsText

# ╔═╡ ec73a428-2ecf-4538-afc7-a6e172d4a632
md"""
Or just call `image()` to display the model as an image:
"""

# ╔═╡ c70cf9d4-8231-4591-a8fe-ec91bcbb41ff
image(fimg_mod_convolved, colormap=:inferno, colorscale=SymLog(5e-4, vmin=0))

# ╔═╡ 7014ff5c-8975-4219-bae2-b434bbe8441f
md"""
# Visibility UV FITS
"""

# ╔═╡ b0f8bdd5-131d-43a2-845d-5785ca6274f5
md"""
Another common file format in VLBI is "uvfits": it contains interferometric visibilities.

`VLBI.load` is the only pure-Julia package that can read such files:
"""

# ╔═╡ aad0d6b7-a58e-45f5-83b4-f73a23e5c559
uvfile = VLBI.load("../test/data/vis.fits");

# ╔═╡ 9ae2c18c-84fe-46e2-b19f-0abb1eccac40
md"""
The `uvfile` object only contains metadata extracted from the FITS header and from auxilary tables:
"""

# ╔═╡ d8f74ef4-fc37-4450-87ec-ef7d55075488
uvfile.ant_arrays

# ╔═╡ d971425a-785e-4424-8630-ded6bba69809
uvfile.freq_windows

# ╔═╡ 442febcb-b7de-4fc4-9a1d-53bb67204092
uvfile.header.stokes

# ╔═╡ 141f6f69-04ce-4bce-b053-3c14791347f8
md"""
... and more. Dump the object to see all fields:
"""

# ╔═╡ 269740dc-af7d-4f3b-bc86-cc285549bc82
uvfile

# ╔═╡ 546825e2-5f27-47ad-b03a-546a082cd250
md"""
The actual visibility data can be retrieved in several formats. The most useful is likely the table format:
"""

# ╔═╡ ce74a7f5-06be-44c9-9b45-f64dc83e14a8
uvtbl_full = VLBI.table(uvfile)

# ╔═╡ faf6e8e0-b79d-4e01-94fc-433817c83fd4
md"""
This gives a `Tables.jl`-compatible table. For now it's a `StructArray`, but the concrete type may potentially change.

Fields should be self-explanatory. Let's plot total intensity visibilities:
"""

# ╔═╡ 95c91cc0-5236-42f2-a031-a137313850f4
uvtbl = @p uvtbl_full |> filter(_.stokes ∈ (:RR, :LL))

# ╔═╡ 4db9b6c4-629e-415f-9730-c7f9bef52969
let
	fig, ax, plt = scatter(RadPlot(uvtbl), markersize=1, color=:black, axis=(xaxisposition=:top,))
	ax, plt = scatter(fig[2,1], RadPlot(uvtbl, yfunc=rad2deg∘angle), markersize=1, color=:black)
	rowsize!(fig.layout, 2, Auto(0.5))
	fig
end

# ╔═╡ 77237980-eda9-4e7c-a3ef-25a0425d0149
md"""
# Source model files
"""

# ╔═╡ f946e665-fce7-4b90-a0a0-283ce60eb14b
md"""
Source models, typically consisting of several circular or elliptical Gaussians, can be stored in `.mod` files produced by software such as `difmap`.

They are loaded with `VLBI.load` as well:
"""

# ╔═╡ 5848e406-18cf-4f36-922b-104d6edd2cf5
dmod = VLBI.load("../test/data/difmap_model.mod")

# ╔═╡ b2a45d8a-c0d6-4e27-8f57-0ff9db75ecc8
md"""
The `MultiComponentModel` type, as well as components themselves, are defined in the `InterferometricModels.jl` package. See its help for more details on how to use these models.
"""

# ╔═╡ e7500919-4325-47b8-87a4-c668f4e43576
md"""
Here's how one can display the model components in the image plane:
"""

# ╔═╡ eab32eff-e4df-4465-a19d-d807266fd699
let
	poly(dmod, alpha=0.5)
	scatter!(dmod, color=:black)
	current_figure()
end

# ╔═╡ 7dada6fd-757c-42e9-964b-f62b7b45503a
image(
	convolve(dmod, beam(CircularGaussian, σ=0.05u"mas")),  # model itself: need to use convolve(), otherwise cannot compute brightness for delta functions
	colorscale=SymLog(1e-5), colormap=:inferno,  # regular Makie arguments
)

# ╔═╡ bd6255c9-fe17-4a29-8a5c-d9b734ae7c5c
md"""
We can also plot amplitude and phase envelopes for each component:
"""

# ╔═╡ 9b0f96da-6326-4e04-8956-ffef0cd9a1e7
let
	fig = Figure()
	ax_ampl = Axis(fig[1,1], xautolimitmargin=(0,0))
	ax_phas = Axis(fig[2,1], xautolimitmargin=(0,0))
	rowsize!(fig.layout, 2, Auto(0.5))
	colors = Makie.wong_colors()
	for (i, c) in enumerate(components(dmod))
		bandstroke!(ax_ampl, RadPlot(0..1e9, model=c), color=(colors[i], 0.1), strokecolor=colors[i])
		bandstroke!(ax_phas, RadPlot(0..1e9, model=c, yfunc=rad2deg∘angle), color=(colors[i], 0.1), strokecolor=colors[i])
	end
	fig
end

# ╔═╡ 47d146c8-47f0-45b4-acf7-7e1c080d22f2
md"""
Compare amplitudes and phases from this model with the observed visibilities:
"""

# ╔═╡ fa4d0551-c1c5-451a-9d0b-7d1f82a688b8
let
	fig = Figure()
	ax_ampl = Axis(fig[1,1])
	ax_phas = Axis(fig[2,1])
	rowsize!(fig.layout, 2, Auto(0.5))
	scatter!(ax_ampl, RadPlot(uvtbl), markersize=1, color=:black)
	scatter!(ax_phas, RadPlot(uvtbl, yfunc=rad2deg∘angle), markersize=1, color=:black)
	scatter!(ax_ampl, RadPlot(uvtbl, model=dmod), markersize=3, color=:red)
	scatter!(ax_phas, RadPlot(uvtbl, model=dmod, yfunc=rad2deg∘angle), markersize=3, color=:red)
	fig
end

# ╔═╡ e8e64ce0-d4ba-4745-b971-d170169a92db
md"""
# Compare with `difmap`
"""

# ╔═╡ 50863255-5341-4c79-8e91-3e9d3a5b20c2
md"""
Visually compare the `VLBIFiles` plot with the one produced by `difmap`. They are based on the same UV files and models, so should look the same:
"""

# ╔═╡ 95d32985-b924-43eb-85d2-115c95b3798c
uvtbl_difmap = @p uvtbl_full |> filter(_.stokes == :RR) |> map(
	_.uv.u >= 0 ?  # difmap puts all uv points into the right semiplane
		(;_.uv, _.visibility) :
		(;uv=-_.uv, visibility=conj(_.visibility))
)

# ╔═╡ 616ec291-8ff1-4112-afac-798dd7790ef7
let
	fig = Figure(size=(500, 600))
	ax_ampl = Axis(fig[1,1])
	ax_phas = Axis(fig[2,1])
	scatter!(ax_ampl, RadPlot(uvtbl_difmap), markersize=1, color=:lime)
	scatter!(ax_phas, RadPlot(uvtbl_difmap, yfunc=rad2deg∘angle), markersize=1, color=:lime)
	scatter!(ax_ampl, RadPlot(uvtbl_difmap, model=dmod), markersize=3, color=:red)
	scatter!(ax_phas, RadPlot(uvtbl_difmap, model=dmod, yfunc=rad2deg∘angle), markersize=3, color=:red)
	fig
end

# ╔═╡ ebf1c8c4-d5c4-425b-baff-42252c0d52cf
let
	pltnames = ["radplot", "uvplot"]
	tdir = mktempdir()
	dr = Difmap.execute(
		"""
		observe vis.fits
		select RR
		rmodel mod.mod
		rflags="m3"
		device "radplot.ps/VCPS"
		radplot "", 0, 0, 0, 0, -30, +30
		device "uvplot.ps/VCPS"
		uvplot
		exit
		""";
		in_files=[
			"../test/data/vis.fits" => "vis.fits",
			"../test/data/difmap_model.mod" => "mod.mod",
		],
		out_files=["radplot.ps", "uvplot.ps"] .=> tempname.(),
	)
	@assert success(dr)
	Difmap.plots(dr, `-density 70`)[1]
end

# ╔═╡ 56e2fabb-5671-4405-9933-43e9ea80fada
md"""
Compare some visibility statistics to ensure that `VLBIFiles` reads the same data as `difmap` does:
"""

# ╔═╡ 47da1b40-3624-4dfd-a965-d2fe2286d7c1
@p uvtbl_difmap |>
	map(angle(_.visibility) |> rad2deg) |>
	(length(__), mean(__), std(__), extrema(__))

# ╔═╡ 25dd465b-7e13-4c20-a2d2-a18c8a6ab11c
let
	dr = Difmap.execute(
		"""
		observe vis.fits
		select RR
		vis_stats amplitude
		vis_stats phase
		vis_stats real
		vis_stats imaginary
		exit
		""";
		in_files=[
			"../test/data/vis.fits" => "vis.fits",
		],
	)
	@p Difmap.inout_pairs(dr) |>
		Dict |>
		__["vis_stats phase"]
end

# ╔═╡ 6f4ecab0-2125-4476-bf78-dbe253322233


# ╔═╡ 90b8e37c-1def-438d-981b-e488cc87a3e7


# ╔═╡ c572cf80-5798-4296-97dc-1428b10bb8cb
TableOfContents()

# ╔═╡ Cell order:
# ╟─5491a70e-8aea-42bc-b449-a56fd120d3fd
# ╟─6cb0d385-5b13-4afe-bf8c-015938efc91d
# ╟─b685e8e4-efa5-453c-a27b-4453946a0839
# ╟─9af33a02-67f8-4233-93f5-4797ecf702d9
# ╟─3df0333d-734c-4d4b-bdec-49cac59f4681
# ╠═6bf4ae14-33d7-4faa-a664-80889ece70b7
# ╟─50c8dfe9-db54-4239-a14d-c8e047eab443
# ╠═8b722ca5-3198-4c6c-94b1-4daf8fd2e718
# ╠═9abadd95-1aeb-436f-aab2-64dc2ad78830
# ╟─b104bdf4-cdee-4743-a4e7-31844f5338e7
# ╠═dcd07de0-148d-4641-9cfc-79926e1ba3b1
# ╠═203e7c90-11d3-4fc1-895a-ffcba4247b34
# ╟─3075f803-d583-4764-a7d4-1bc14d8e441b
# ╠═b32c9639-c2c3-4a67-8554-57aefea9d196
# ╠═530aa386-84a7-47b2-ba76-a9c01fc5393d
# ╟─6d0bf04e-51ee-43b6-8b7c-912fd051dfb5
# ╟─c4b58f30-e7db-459e-9d26-2d4bc57e5fbd
# ╠═b4e14571-a3fa-41d2-98da-cff594f88202
# ╟─c8510200-f820-4b2a-8260-d21e8f7f7880
# ╟─ebcabc9f-0487-4cad-b7f7-fff448fcf203
# ╠═7994f45f-6e01-4a62-ae3b-3031b46f2d7b
# ╟─9db795b2-2d79-45b6-a51c-bc09f4d05150
# ╟─0a054aab-43d0-4caf-85c9-9303dec7c8fb
# ╠═2f424454-1264-4887-927c-ad58d8631ecc
# ╟─ec73a428-2ecf-4538-afc7-a6e172d4a632
# ╠═c70cf9d4-8231-4591-a8fe-ec91bcbb41ff
# ╟─7014ff5c-8975-4219-bae2-b434bbe8441f
# ╟─b0f8bdd5-131d-43a2-845d-5785ca6274f5
# ╠═aad0d6b7-a58e-45f5-83b4-f73a23e5c559
# ╟─9ae2c18c-84fe-46e2-b19f-0abb1eccac40
# ╠═d8f74ef4-fc37-4450-87ec-ef7d55075488
# ╠═d971425a-785e-4424-8630-ded6bba69809
# ╠═442febcb-b7de-4fc4-9a1d-53bb67204092
# ╟─141f6f69-04ce-4bce-b053-3c14791347f8
# ╠═269740dc-af7d-4f3b-bc86-cc285549bc82
# ╟─546825e2-5f27-47ad-b03a-546a082cd250
# ╠═ce74a7f5-06be-44c9-9b45-f64dc83e14a8
# ╟─faf6e8e0-b79d-4e01-94fc-433817c83fd4
# ╠═95c91cc0-5236-42f2-a031-a137313850f4
# ╠═4db9b6c4-629e-415f-9730-c7f9bef52969
# ╟─77237980-eda9-4e7c-a3ef-25a0425d0149
# ╟─f946e665-fce7-4b90-a0a0-283ce60eb14b
# ╠═5848e406-18cf-4f36-922b-104d6edd2cf5
# ╟─b2a45d8a-c0d6-4e27-8f57-0ff9db75ecc8
# ╟─e7500919-4325-47b8-87a4-c668f4e43576
# ╠═eab32eff-e4df-4465-a19d-d807266fd699
# ╠═7dada6fd-757c-42e9-964b-f62b7b45503a
# ╟─bd6255c9-fe17-4a29-8a5c-d9b734ae7c5c
# ╠═9b0f96da-6326-4e04-8956-ffef0cd9a1e7
# ╟─47d146c8-47f0-45b4-acf7-7e1c080d22f2
# ╠═fa4d0551-c1c5-451a-9d0b-7d1f82a688b8
# ╟─e8e64ce0-d4ba-4745-b971-d170169a92db
# ╟─50863255-5341-4c79-8e91-3e9d3a5b20c2
# ╠═95d32985-b924-43eb-85d2-115c95b3798c
# ╠═616ec291-8ff1-4112-afac-798dd7790ef7
# ╠═ebf1c8c4-d5c4-425b-baff-42252c0d52cf
# ╟─56e2fabb-5671-4405-9933-43e9ea80fada
# ╠═47da1b40-3624-4dfd-a965-d2fe2286d7c1
# ╠═25dd465b-7e13-4c20-a2d2-a18c8a6ab11c
# ╠═6f4ecab0-2125-4476-bf78-dbe253322233
# ╠═90b8e37c-1def-438d-981b-e488cc87a3e7
# ╠═4eb97406-7713-11ec-00ca-8be6bef77030
# ╠═ead10789-68f5-4931-b177-622030df739d
# ╠═0c44dd70-5a5c-46b9-b085-1a52c22078ca
# ╠═17a5e5ad-c3a7-41a0-bf56-410675e95475
# ╠═7a4e96c8-8ab0-4338-82b5-e07f0afdaae5
# ╠═9c6efe8c-385e-4446-acdf-bd19cffe31e2
# ╠═be208a3e-9a03-4194-8be7-05547d9871a7
# ╠═4e7278e0-0d12-4762-8ace-e7cbe8eb371c
# ╠═dc179b8b-a5b7-4d2d-bf86-afc0fc5a376b
# ╠═f3898d40-517d-421b-b551-00d8ce1f64dd
# ╠═8f556fb7-c181-4e0f-a42d-1caa85f35f41
# ╠═398b65dd-67ae-4fa8-92d4-a73cf2992060
# ╠═412dc1f1-ea14-40ee-ad3b-28b89990738f
# ╠═7580f512-0653-479d-854b-f5aa78d85655
# ╠═c572cf80-5798-4296-97dc-1428b10bb8cb
# ╠═644d9d5f-9b1c-4dd2-b688-faac9e7c2473
# ╠═993027e8-7f5b-4952-8ed7-0345a9522d3a
