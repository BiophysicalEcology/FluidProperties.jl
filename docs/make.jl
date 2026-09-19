using Documenter
using DocumenterVitepress
using FluidProperties
using CairoMakie
using Unitful

# Don't output huge svgs for Makie plots
CairoMakie.activate!(type = "png")

# RasterDataSources downloads go here unless the path is set
get!(ENV, "RASTERDATASOURCES_PATH", joinpath(@__DIR__, "data"))

# Helpers for the figures, loaded in the examples with `using Main.FigureHelpers`
include("figure_helpers.jl")

makedocs(
    modules = [FluidProperties],
    sitename = "FluidProperties.jl",
    authors = "Michael Kearney, Rafael Schouten et al.",
    clean = true,
    doctest = false,
    checkdocs = :exports,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/BiophysicalEcology/FluidProperties.jl", # this must be the full URL!
        devbranch = "main",
        devurl = "dev";
    ),
    source = "src",
    build = "build",
    warnonly = true,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/BiophysicalEcology/FluidProperties.jl",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
