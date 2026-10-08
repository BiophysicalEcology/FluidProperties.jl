# Properties of Air: A Manual for Use in Biophysical Ecology

*C. Richard Tracy, William R. Welch, Berry Pinshow, Michael R. Kearney and Warren P. Porter*

## Preface to the Julia port

This manual is the documentation of FluidProperties.jl. It is based on the fifth edition of
*Properties of Air*, which was transcribed into R Markdown as a vignette of the
[NicheMapR](https://github.com/mrke/NicheMapR) R package, with the Fortran subroutines and function converted to R
functions. Here the text has been carried over with the R functions replaced by their Julia equivalents,
which use [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) so that units are attached to
all inputs and outputs. The figures are generated with the Julia functions, and the code that produces each one
is shown below it. Sections on the wet bulb temperature and on liquid water have been added.

## Preface to the Fifth Edition

This fifth edition of Properties of Air has been transcribed into R Markdown, and the Fortran subroutines and function have been converted to R functions, by Michael Kearney. The figures have been generated with the R functions. It is now a vignette as part of the [NicheMapR R package](https://github.com/mrke/NicheMapR).

## Preface to the Fourth Edition

This fourth edition of Properties of Air is largely a reprint of the third edition with errors corrected. This edition was created because Berry Pinshow wanted a digital version for his personal use. He scanned the third edition, and corrected errors in it. The fourth edition exists due to Berry's enthusiasm to develop a digital version with errors
corrected. This preface is largely to say a few things about this small book.

When I was a graduate student, I found myself constantly looking up physical characteristics of air as I developed biophysical ecological models. I got so tired of looking in reference books that I decided to extract what I needed and make tables and graphs of the properties of air. These, and other tables that I accumulated, became my main source of reference material in my research. For years, I carried these reference tools around in a folder of hand-drawn graphs and tables. Various visitors to Warren Porter's lab asked for copies of the graphs and tables in my folder. Bill Welch suggested that we get the artistic staff in the Department of Zoology at the University of Wisconsin to make pen-and-ink drawings from my graphs. Subsequently, I created a set of general subroutines in the program language FORTRAN, which simply translated the graphs. As these tools accumulated, Bill Welch and I continued to add to the reference material, and Warren Porter facilitated getting the book put together in a package that we self-published out of the Porter Lab. Many dozens of people have requested this book, and it seems to remain useful as many continue to ask for copies of it. However, the original printing has now exhausted, so creating this digital version made sense.

One of the original authors of this book, and our good friend, Dr. Bill Welch, died at
much too young an age. To him, we dedicate this fourth edition.

C. Richard Tracy  
Department of Biology  
University of Nevada, Reno  
Reno, NV 89557

First Edition 1973  
Second Edition 1978  
Third Edition 1980  
Fourth Edition 2010

## Preface to the Third Edition

This manual was made possible by the invaluable assistance of Ann Chambers and Cheryle Hughes. Financial support was provided by grants from the Department of Zoology, Wisconsin Alumni Research Foundation, ERDA (Contract EY-76-5-02-2270), and NSF (Grant Nos. 74-19454 and 77-25786 to WPP).

First edition 1973  
Second edition 1978  
Third edition 1980  
Fourth edition 2010

## Introduction

This manual comprises a series of tables and graphs that illustrate selected properties of air as functions of temperature, pressure, and humidity. The particular properties that are displayed were chosen because of their importance in analytical studies of energy (heat) and mass (water) transfer between organisms and their physical environments. Except where otherwise noted, the graphs were drawn from information in the Smithsonian Meteorological Tables (List 1971). These graphs are considered to be a concise references for the values of the indicated properties, and not substitutes for the more accurate Smithsonian tables.

All information in this manual is expressed in terms of the International System of Units (SI), which is based, in part, on the metre, kilogram, second and kelvin. Basic units, derived units, prefixes, and some useful physical constants and conversion factors are listed in [Units and constants](units_constants.md). More extensive lists are provided by List (1971), Mechtly (1973), and the Symbols Committee of the Royal Society (1975). Note that the symbol for a physical quantity is an italicized (slanted) letter of the Roman or Greek alphabet, and the unit for a quantity is indicated by an unitalicized (upright) Roman letter.

An equation is provided as a mathematical description of each graph. In addition, the functions [`dry_air_properties`](@ref) and [`wet_air_properties`](@ref) calculate all of the information in the graphs and are included in FluidProperties.jl, and the Julia code used to obtain the data in the plots is provided below each figure. Two of the equations in [`dry_air_properties`](@ref) (thermal conductivity of air and latent heat of vaporization of water) are linear regressions that were fit to tabular data. All other equations in [`dry_air_properties`](@ref) were taken from the Smithsonian tables.

All humidity calculations in [`wet_air_properties`](@ref) are made with the internationally accepted Goff-Gratch equation ([`GoffGratch`](@ref), used by default by [`vapour_pressure`](@ref)). A mathematically less cumbersome and often more useful alternative was reported by Tetens (Murray 1976), available as [`Teten`](@ref); calculations with this equation differ negligibly from those with the Goff-Gratch expression. Teten's equation is provided in [Figure 3](humid_air.md#Figure-3.-Dew-point-temperature) to calculate the saturation vapor pressure ($e_d^*$) at the dewpoint temperature ($t_d$). Note that this equation can also be used to calculate the saturation vapor pressure ($e^*$) at the dry bulb temperature ($t$) if $t$ replaces $t_d$, and the saturation vapor pressure ($e_w^*$) at the wet bulb temperature ($t_w$) if $t_w$ replaces $t_d$. The relationship of wet bulb temperature to dry bulb temperature, vapor density, and relative humidity is shown in the psychometric diagram [Figure 15](humid_air.md#Figure-15.-Vapour-density-over-water).

Mathematical expressions for some properties of air that are dependent on humidity are not available. The equations for these properties are restricted to dry air in this manual but see Mason and Monchick (1965) for further information on humid air.

## Working with units

Every function takes and returns [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) quantities, so any
compatible units can be used and mistakes such as passing a percentage for a fraction or hectopascals for pascals are avoided.
Temperatures may be in kelvin or degrees Celsius, pressures in any pressure unit, and relative humidity is a fraction
between 0 and 1.

```@example intro
using FluidProperties, Unitful

air = dry_air_properties(20.0u"°C", 101325.0u"Pa")
air.density
```

Converting units is the job of Unitful:

```@example intro
uconvert(u"g/cm^3", air.density)
```

## References

Mason EA, Monchick L. 1965. Survey of the equation of state and transport properties of moist gases. Pages 257-272 in
Wexler A, ed. Humidity and Moisture. Measurement and Control in Science and Industry, vol. 3. New York: Reinhold Publishing Corporation.

Mechtly EA. 1973. The International System of Units: Physical Constants and Conversion Factors. National Aeronautics and Space Administration.

Murray BR. 1976. On the computation of saturation vapor pressure. Journal of Applied Meteorology 6:203-204.

List RJ. 1971. Smithsonian Meteorological Tables. Smithsonian Institution Press.

Symbols Committee of the Royal Society. 1975. Quantities, Units and Symbols. The Royal Society.
