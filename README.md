# EllpacaEscapeAtlas

Here we generate highlighter plots and escape logos for the Ellpca alignments
and then combine them into one svg plot per donor:

First edit the4 Highlighter and EscapeLogo scripts to select the
appropriate epitope, input and output directories and then run as follows:

```
julia HighlighterScript.jl
julia EscapeLogoScript.jl
julia combine_svg_plots.jl CAP256_induced_highlighter_plots CAP256_escape_logo_plots CAP256_induced_combined_plots 20
```
