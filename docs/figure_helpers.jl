module FigureHelpers

using CairoMakie

export figure_axis

"""
    figure_axis(xlabel, ylabel; size=(700, 500), kw...)

A `Figure` and `Axis` with minor ticks and grid lines, as used for the figures of the manual.
"""
function figure_axis(xlabel, ylabel; size=(700, 500), kw...)
    fig = Figure(; size)
    ax = Axis(fig[1, 1];
        xlabel, ylabel,
        xminorticksvisible=true, yminorticksvisible=true,
        xminorgridvisible=true, yminorgridvisible=true,
        xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
        kw...,
    )
    return fig, ax
end

end
