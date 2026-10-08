module FluidPropertiesEnzymeCoreExt

using FluidProperties
using EnzymeCore: EnzymeRules

# Enzyme can't differentiate the Roots.jl search, so skip it.
# Derivatives come from the Newton steps in `_find_temperature`.
EnzymeRules.inactive(::typeof(FluidProperties._bracketed_root), args...) = nothing

end
