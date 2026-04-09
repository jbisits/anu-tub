# Setup up plotting style

"""
    function latex_format(values)
Generate latex tick values for a `Colorbar`.
"""
function latex_format(values)
    map(values) do v
        if v == 0
            return L"0"
        else
            # Extract coefficient and exponent
            formatted = @sprintf("%.1e", v)
            coeff, exp = split(formatted, 'e')
            exponent = parse(Int, exp)
            return L"%$(coeff) \times 10^{%$(exponent)}"
        end
    end
end
markersize = 10
publication_theme = Theme(font="CMU Serif", fontsize = 20,
                          Axis=(titlesize = 22,
                                xlabelsize = 20, ylabelsize = 20,
                                xgridstyle = :dash, ygridstyle = :dash,
                                xtickalign = 0, ytickalign = 0,
                                yticksize = 6.5, xticksize = 6.5),
                          Legend=(framecolor = (:black, 0.5),
                                  backgroundcolor = (:white, 0.5),
                                  labelsize = 20),
                          Colorbar=(ticksize=12,
                                    tickalign=1,
                                    spinewidth=0.5))
new_theme = merge(theme_latexfonts(), publication_theme)
set_theme!(new_theme)
