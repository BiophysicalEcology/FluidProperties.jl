```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "FluidProperties.jl"
  text: "Properties of air and water"
  tagline: "physical properties of air, water vapour and liquid water for biophysical ecology, with units."
  actions:
    - theme: brand
      text: Get Started
      link: /get_started
    - theme: alt
      text: View on Github
      link: https://github.com/BiophysicalEcology/FluidProperties.jl
    - theme: alt
      text: API Reference
      link: /api

features:
  - title: 🌬️ Air
    details: Density, viscosity, thermal conductivity, vapour diffusivity, black-body emittance and more for <a class="highlight-link">dry air</a> and <a class="highlight-link">humid air</a> as functions of temperature, pressure and humidity, from the Smithsonian Meteorological Tables.
    link: /manual/dry_air
  - title: 💧 Vapour pressure
    details: Several equations for <a class="highlight-link">saturation vapour pressure</a>, including Goff-Gratch, Tetens, Huang and Bolton, and a fast lookup table.
    link: /manual/vapour_pressure
  - title: 🌡️ Wet bulb temperature
    details: <a class="highlight-link">Wet bulb temperature</a> by the thermodynamic method of Davies-Jones, psychrometric equations, an energy balance and empirical fits, and the natural wet bulb temperature.
    link: /manual/wet_bulb
  - title: 📏 Units
    details: Every input and output uses <a class="highlight-link">Unitful.jl</a>, so any compatible units can be used and unit mistakes are avoided. Units compile away, and the functions are type-stable and non-allocating.
    link: /manual/introduction
  - title: 🗺️ Rasters
    details: The functions broadcast over gridded data, so properties of air can be mapped with <a class="highlight-link">Rasters.jl</a> and <a class="highlight-link">RasterDataSources.jl</a>.
    link: /tutorials/rasters
---
```

## How to install FluidProperties.jl?

FluidProperties.jl can be installed from the Julia REPL:

```julia
julia> using Pkg
julia> Pkg.add("FluidProperties")
# or
julia> ] # ']' should be pressed
pkg> add FluidProperties
```

If you want to use the latest unreleased version, you can run the following command:

```julia
julia> using Pkg
julia> Pkg.add(url = "https://github.com/BiophysicalEcology/FluidProperties.jl")
```

## Manual

This documentation is based on *Properties of Air: A Manual for Use in Biophysical Ecology* by Tracy, Welch, Pinshow, Kearney
and Porter. See the [Introduction](manual/introduction.md) for how it relates to the package.
