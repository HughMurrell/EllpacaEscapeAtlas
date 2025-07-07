# EllpacaEscapeAtlas

Here we generate annotated Ellpaca tree, highlighter charts and 
escape logos for each donor Ellpaca functional alignment,
and then combine them into one document using LaTex.

First, install `julia` and the relavent packages as usual then
edit the three scripts below to select the appropriate input and 
output directories and then run the scripts as follows:

```
julia EllpacaTreeScript.jl
julia EscapeAtlasScript.jl
julia EscapeLogoScript.jl  # run this one for each epitope
pdflatex EscapeAtlas.tex
```
