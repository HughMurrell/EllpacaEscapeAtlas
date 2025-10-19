using Pkg
Pkg.activate("./")
Pkg.instantiate()

using CSV, FASTX, DataFrames, Dates, BioSequences
using StatsBase
using Statistics: mean
using Luxor
using Colors
using FreeTypeAbstraction
using FreeType

@show pwd()

# functions from SequenceLogo.jl

# first logo.jl

struct WeightedLetter
    letter::Char
    weight::Float64
end

struct SequenceLogoSite
    weighted_letters::Vector{WeightedLetter}
end

struct SequenceLogo
    sites::Vector{SequenceLogoSite}
end

"""
    sort_letters(site)

Sorts the letters in a site in order of increasing absolute weight.
"""
function sort_letters(site::SequenceLogoSite)
    sorted_letters = sort(site.weighted_letters; by=(l -> abs(l.weight)))
    return SequenceLogoSite(sorted_letters)
end

"""
    sort_letters(logo; remove_duplicates=true)

Sorts the letters in all sites in order of increasing absolute weight.
"""
function sort_letters(logo::SequenceLogo; remove_duplicates=true)
    if remove_duplicates
        logo = remove_duplicate_letters(logo)
    end
    sorted_sites = [sort_letters(site) for site in logo.sites]
    return SequenceLogo(sorted_sites)
end

"""
    remove_duplicate_letters(logo)

Removes duplicate letters at each site (summing their weights).
"""
function remove_duplicate_letters(logo::SequenceLogo)
    dedup_sites = Vector{SequenceLogoSite}()
    for site in logo.sites
        d = Dict{Char,Float64}()
        for weighted_letter in site.weighted_letters
            l, w = weighted_letter.letter, weighted_letter.weight
            d[l] = get(d, l, zero(w)) + w
        end
        site = SequenceLogoSite([WeightedLetter(l,w) for (l,w) in d])
        push!(dedup_sites, site)
    end
    return SequenceLogo(dedup_sites)
end

"""
    logo_from_matrix(weights, alphabet)

Construct a `SequenceLogo` from a weight matrix of size (q, L), where `q` is
the number of possible letters and `L` the length of the sequence.
`alphabet` is a `String` containing the `q` possible letters.
"""
function logo_from_matrix(w::AbstractArray, alphabet::String)
    @assert size(w, 1) ≤ length(alphabet)
    sites = Vector{SequenceLogoSite}(undef, size(w, 2))
    for i = 1:size(w, 2)
        letters = [WeightedLetter(alphabet[a], w[a,i]) for a = 1:size(w, 1)]
        sites[i] = SequenceLogoSite(letters)
    end
    return SequenceLogo(sites)
end


normalize(x::AbstractArray) = x ./ sum(x; dims = 1)

"""
    conservation_scores(P; unit = log(2))

Conservation scores at each site. The default unit = log(2) gives results in bits.
"""
function conservation_scores(P::AbstractMatrix; unit::Real = log(2))
    return log2(size(P, 1)) .+ sum(xlogx.(normalize(P)); dims = 1) ./ unit
end

"""
    conservation_scores(P, M; unit = log(2))

Similar to `conservation_scores(P)`, but instead of a flat measure over
letters at each site, computes a KL divergence to a prior measure `M`.
"""
function conservation_scores(P::AbstractMatrix, M::AbstractMatrix; unit::Real = log(2))
    p = normalize(P)
    m = normalize(M)
    KL = sum(p .* xlogy.(p ./ m); dims=1) ./ unit
    return KL
end

"""
    conservation_matrix(P; unit = log(2))

Given a matrix of frequencies `P` of size (q, L), where `q` is the number
of possible letters and `L` the sequence length, returns the matrix of
conservations (as defined by Schneider and Stephens 1990, 10.1093/nar/18.20.6097).
"""
function conservation_matrix(P::AbstractMatrix; unit = log(2))
    return normalize(P) .* conservation_scores(P; unit = unit)
end

"""
    conservation_matrix(P, M; unit = log(2))

Similar to `conservation_matrix(P)`, but instead of a flat measure over
letters at each site, computes a KL divergence to a prior measure `M`.
"""
function conservation_matrix(P::AbstractMatrix, M::AbstractMatrix; unit = log(2))
    return normalize(P) .* conservation_scores(P, M; unit = unit)
end

function xlogx(x::Real)
    result = x * log(x)
    return ifelse(iszero(x), zero(result), result)
end

function xlogy(x::Real, y::Real)
    result = x * log(y)
    ifelse(iszero(x), zero(result), result)
end

# second color_ftns.jl

function aa_color(aa::Char)
    if aa in ('G', 'S', 'T', 'Y', 'C', 'Q', 'N') # polar
        return "green"
    elseif aa in ('K', 'R', 'H') # basic
        return "blue"
    elseif aa in ('D', 'E') # acidic
        return "red"
    elseif aa in ('A', 'V', 'L', 'I', 'P', 'W', 'F', 'M') # hydrophobic
        return "orange"
    else
        return "black"
    end
end

function nt_color(nt::Char) 
    if nt == 'G'
        return "orange"
    elseif nt in ('T', 'U')
        return "red"
    elseif nt == 'C'
        return "blue"
    elseif nt == 'A'
        return "green"
    else
        return "black"
    end
end

function visit_color(i)
    all_colors = distinguishable_colors(20)
    if ( isnothing(i) || i > 17 )
        return(all_colors[1])
    end
    return(all_colors[i+3])
end

# third bio.jl

aa_letters() = "ACDEFGHIKLMNPQRSTVWY-"

"""
    aa2int(char)

Amino-acid one-letter code to integer (1 to 21, with gap = 21).
"""
function aa2int(c::Union{Char,UInt8})
    # from: https://github.com/carlobaldassi/GaussDCA.jl
    alphabet = Int8.((1,21, 2, 3, 4, 5, 6, 7, 8,21, 9,10,11,12,21,13,14,15,16,17,21,18,19,21,20))
                    # A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y
    i = UInt8(c) - 0x40
    if 1 ≤ i ≤ 25
        return alphabet[i]
    else
        return Int8(21)
    end
end

"""
    int2aa(i)

Get the amino acid (in one letter code) corresponding to the integer `i`.
"""
int2aa(i::Integer) = alphabet_aa()[i]

"""
    aa2int(seq::String)

Convert a string amino-acid sequence to a vector of integers.
"""
aa2int(s::String) = aa2int.(collect(s))

"""
    aa2int(seq::String[])

Given a vector of string amino-acid sequences, converts it
to an integer matrix.
"""
aa2int(ss::AbstractVector{String}) = hcat(aa2int.(ss)...)


const GAP = '∇'
aa_alphabet() = replace(aa_letters(), '-' => GAP)
logo_from_matrix_aa(w::AbstractMatrix) = logo_from_matrix(w, aa_alphabet())



alphabet_dna() = "ACGT-"

aas=aa_alphabet() * "|"

function onehot(s::AbstractString)
    return reshape(collect(s), 1, length(s)) .== collect(aas)
end

function luxor_letter_at(ch::Char, color, yscale::Real, font_size::Real)
    sethue(color)
    w=textextents(string(ch))[3]
    h=textextents(string(ch))[4]
    font_scale=yscale*font_size/h
    # scale(100/w,100/h)
    if ch=='|'
        setdash("solid")
        scale(1,2)
        textoutlines(string(ch), O, :path, valign=:bottom, halign=:left)
        fillpreserve()
        strokepath()
        scale(1,1/2)
        # strokepath()
        # line(Point(font_size/2, 0),Point(font_size/2,-yscale*font_size),action = :stroke)
        setdash("solid")
    else
        scale(1,font_scale)
        textoutlines(string(ch), O, :path, valign=:bottom, halign=:left)
        fillpreserve()
        strokepath()
        scale(1,1/font_scale)
    end
    sethue("black")
    return 
end

function luxor_sequence_logo_aa(ref_logo::Union{Nothing, SequenceLogo},
                                bot_logo::SequenceLogo, 
                                # mid_logo::Union{Nothing, SequenceLogo},
                                freq_logos::Vector{Any},
                                color_fun, logo_length, font_size, path_svg;
                                thresh=0, remove_duplicate_letters=true, annot=[], 
                                scale_factor=1,exploded=true,
                                title="Escape Logo", 
                                ref_annot="",
                                bot_annot="consensus", 
                                freq_annots=[])
    lh=40
    sym_gap=2/5
    freq_scale=2
    offset=14*font_size*(1-sym_gap)  # used to be 13*
    println("generating logo in $(path_svg)")
    p=Drawing((logo_length+14)*font_size*(1-sym_gap), 100 + lh*length(freq_logos) + 75, path_svg, strokescale=true)
    fontface("Menlo")   # Punch, Menlo, Courier
    fontsize(font_size)
    setline(1)
    gsave()
    text(title,Point(font_size,3*font_size))
    fontsize(font_size/2)
    for freq_ind in 1:length(freq_annots)
        color_ind=1+length(freq_annots)-freq_ind
        setcolor(visit_color(color_ind))
        text(freq_annots[freq_ind],Point(font_size,80+lh*freq_ind-font_size))
    end
    @show bot_annot
    text(bot_annot,Point(font_size,80+lh*(length(freq_logos)+1)-lh/2))
    text(ref_annot,Point(font_size,130+font_size))
    fontsize(font_size)
    
    for freq_ind in 1:length(freq_logos)
        # first do the freq logo grid lines
        grestore()
        gsave()
        Luxor.translate(offset,80 + lh*freq_ind + font_size - font_size)
        ygs=[0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
        setline(1)
        setdash("solid")
        setcolor("grey80")
        for yg in ygs
            line(Point(0, -yg*font_size),Point((logo_length)*font_size*(1-sym_gap),-yg*font_size),action = :stroke)
            setdash("dot")
        end
        setdash("solid")
        line(Point(0, -font_size*freq_scale),Point((logo_length)*font_size*(1-sym_gap),-font_size*freq_scale),
                        action = :stroke)
        # now do the freq logo
        sorted_logo = sort_letters(freq_logos[freq_ind]; remove_duplicates=remove_duplicate_letters)
        y_min = y_max = 0.0
        for (x, site) in enumerate(sorted_logo.sites)
            y_pos = y_neg = 0.0
            for weighted_letter in site.weighted_letters
                w = abs(weighted_letter.weight*freq_scale)
                l = weighted_letter.letter
                w ≤ thresh && continue
                if weighted_letter.weight > 0
                    luxor_letter_at(l, color_fun(l), w, font_size)
                    y_pos += w
                    Luxor.translate(0,-w*font_size)
                elseif weighted_letter.weight < 0
                    println("******* negative y shift *************")
                    y_neg -= w
                    luxor_letter_at(l, color_fun(l), w, font_size)
                end
            end
            Luxor.translate(0,y_pos*font_size)
            Luxor.translate((1-sym_gap)*font_size, 0)
        end
    end
    
   
    # now return to origin and do the bottom logo
    grestore()
    gsave()
    Luxor.translate(offset,80 + lh*(length(freq_logos)+1) - lh/2)
    
    sorted_logo = sort_letters(bot_logo; remove_duplicates=remove_duplicate_letters)
    y_min = y_max = 0.0
    annot_stack=reverse(deepcopy(annot))
    for (x, site) in enumerate(sorted_logo.sites)
        ann=pop!(annot_stack)
        if ann > 0
            gsave()
            fontsize(10)
            text(string(ann), Point((font_size/2)-(sym_gap*font_size),5), angle=pi/2)
            grestore()
        end
        y_pos = y_neg = 0.0
        for weighted_letter in site.weighted_letters
            w = abs(weighted_letter.weight) /2
            l = weighted_letter.letter
            w ≤ thresh && continue
            if weighted_letter.weight > 0
                luxor_letter_at(l, color_fun(l), w, font_size)
                y_pos += w
                Luxor.translate(0,-w*font_size)
            elseif weighted_letter.weight < 0
                println("******* negative y shift *************")
                y_neg -= w
                luxor_letter_at(l, color_fun(l), w, font_size)
            end
        end
        Luxor.translate(0,y_pos*font_size)
        Luxor.translate((1-sym_gap)*font_size, 0)
    end

    # now return to origin, translate and do the reference 
    if false
    grestore()
    gsave()
    Luxor.translate(offset,145+font_size)
    # transform([1 0 0 -1 0 0])
    
    sorted_logo = sort_letters(ref_logo; remove_duplicates=remove_duplicate_letters)
    y_min = y_max = 0.0
    annot_stack=reverse(deepcopy(annot))
    for (x, site) in enumerate(sorted_logo.sites)
        ann=pop!(annot_stack)
        if ann > 0
            gsave()
            fontsize(10)
            text(string(ann), Point((font_size/2)-(sym_gap*font_size),5), angle=pi/2)
            grestore()
        end
        y_pos = y_neg = 0.0
        for weighted_letter in site.weighted_letters
            w = abs(weighted_letter.weight) /2
            l = weighted_letter.letter
            w ≤ thresh && continue
            if weighted_letter.weight > 0
                luxor_letter_at(l, color_fun(l), w, font_size)
                y_pos += w
                Luxor.translate(0,-w*font_size)
            elseif weighted_letter.weight < 0
                println("******* negative y shift *************")
                y_neg -= w
                luxor_letter_at(l, color_fun(l), w, font_size)
            end
        end
        Luxor.translate(0,y_pos*font_size)
        Luxor.translate((1-sym_gap)*font_size, 0)
    end
    end

    # grestore()
    # gsave()
    # translate(0,131+font_size)
    # setline(1)
    # setcolor("grey80")
    # line(Point(offset, 0),Point((logo_length+12)*font_size*(1-sym_gap),0),action = :stroke)
    # Point((logo_length)*font_size*(1-sym_gap)
    # grestore()
    # gsave()
    # Luxor.translate(0,131+font_size)
    # setline(1)
    # setcolor("grey80")
    # line(Point(offset, 0),Point((logo_length+14)*font_size*(1-sym_gap),0),action = :stroke)
    grestore()
    gsave()
    Luxor.translate(0,106+font_size)
    setline(1)
    setcolor("grey80")
    # line(Point(offset, 0),Point((logo_length+14)*font_size*(1-sym_gap),0),action = :stroke)
    grestore()
    gsave()
    Luxor.translate(0,80 + lh*(length(freq_logos)+1) -lh/2 + 1)
    setline(1)
    setcolor("grey80")
    line(Point(offset, 0),Point((logo_length+14)*font_size*(1-sym_gap),0),action = :stroke)
    finish()
    return (p)
end

function ref_coords(seq)
    rc = [1]
    for i in 2:length(seq)
        if seq[i] != '-'
            push!(rc,rc[i-1]+1)
        else
            push!(rc,rc[i-1])
        end
    end
    return rc
end

function offset(ref,coords)
    start=coords[1]
    stop=coords[end]
    rcs = ref_coords(ref)
    return (findfirst((x->x==start),rcs):findfirst((x->x==stop),rcs))
end

loop_coords = Dict()
loop_coords["V1"]=131:157



aa_colors = Dict(
     'B'=>RGB{Float64}(1.0,0.630714,0.576563),
     'M'=>RGB{Float64}(0.756869,0.916499,0.965176),
     'I'=>RGB{Float64}(0.187839,0.54561,0.252343),
     'X'=>RGB{Float64}(0.540006,0.493982,0.813567),
     'Y'=>RGB{Float64}(0.0973617,0.285282,0.5329),
     'Z'=>RGB{Float64}(0.0418427,0.156645,0.341597),
     'L'=>RGB{Float64}(0.426131,0.0441442,0.0465628),
     'O'=>RGB{Float64}(0.518954,0.802339,0.930272),
     'F'=>RGB{Float64}(0.587882,0.865532,0.51112),
     'Q'=>RGB{Float64}(0.0,0.225356,0.101282),
     'D'=>RGB{Float64}(0.862653,0.958477,0.981395),
     'V'=>RGB{Float64}(0.188382,0.529206,0.795898),
     'U'=>RGB{Float64}(0.277786,0.635283,0.863472),
     'E'=>RGB{Float64}(0.711814,0.932724,0.629136),
     'T'=>RGB{Float64}(0.394211,0.72627,0.90426),
     'H'=>RGB{Float64}(0.32729,0.673206,0.326717),
     'P'=>RGB{Float64}(0.0232916,0.395886,0.180144),
     'G'=>RGB{Float64}(0.459895,0.779462,0.41097),
     'N'=>RGB{Float64}(0.641543,0.865092,0.94902),
     'K'=>RGB{Float64}(0.825431,0.118066,0.106858),
     'C'=>RGB{Float64}(0.835916,0.980813,0.770886),
     'R'=>RGB{Float64}(0.807625,0.787968,0.949453),
     'W'=>RGB{Float64}(0.137797,0.411028,0.686187),
     'A'=>RGB{Float64}(1.0,0.808314,0.771835),
     'S'=>RGB{Float64}(0.236943,0.0166779,0.407047),
     'J'=>RGB{Float64}(1.0,0.389569,0.336934),
     '*'=>RGB{Float64}(0.3,0.3,0.3),
     '|'=>RGB{Float64}(0.1,0.1,0.1),
     '-'=>RGB{Float64}(1.0,1.0,1.0),
     '.'=>RGB{Float64}(1.0,1.0,1.0),
     '<'=>RGB{Float64}(1.0,1.0,1.0),
     '>'=>RGB{Float64}(1.0,1.0,1.0),
     '∇'=>RGB{Float64}(0.1,0.1,0.1),
     ' '=>RGB{Float64}(1.0,1.0,1.0));

function aa_color_fun(ch)
    if ch=='|' || ch=='∇'
        return("grey80")
    end
    return(aa_color(ch)) #s[ch])
end
(aa_color_fun).(collect(aas))


function generate_logo(fasta_path, ept, ept_name, outdir; file_count=0, prefix="")
    donor=basename(fasta_path)[1:6]
    @show(donor)
    records = collect(FASTX.FASTA.Reader(open(fasta_path)))
    all_seqs = (x->FASTX.sequence(String,x)).(records)
    all_nams = String.(FASTX.identifier.(records))
    gapyness = (x->sum(collect(x).=='-')/length(x)).(all_seqs)
    @show gapyness[1], maximum(gapyness)
    sel_inds = gapyness .<= GAPPY_CUTOFF
    sel_inds[1] = true
    sel_inds[2] = true
    println("keeping $(sum(sel_inds)) out of $(length(all_seqs))")
    all_seqs=all_seqs[sel_inds]
    all_nams=all_nams[sel_inds]
    
    # recompute consensus of first visit and store it
    visits=sort(union((x->split(x,"_")[2][1:4]).(all_nams[3:end])))
    first_inds = (x->split(x,"_")[2][1:4]==visits[1]).(all_nams[3:end])
    consensus_of_first=consensus(all_seqs[3:end][first_inds])
    all_nams[2]="consensus of $(visits[1])"
    all_seqs[2]=consensus_of_first
    
    # close(fasta_path)
    # drop pool 6 sequences in ellpaca data
    # not_pool6=(x->!contains(x,"pool6")).(all_nams)
    # all_nams=all_nams[not_pool6]
    # all_seqs=all_seqs[not_pool6]
    # @show sum(not_pool6), length(all_nams)
    if length(all_seqs) < 3
        println("skipping $(donor), not enough sequences")
        return(logo_frames)
    end
    ept_keys = sort(collect(keys(ept)))
    seqs=["" for seq in all_seqs]
    ept_annot=[]
    for key in ept_keys
        coords=ept[key] #loop_coords["V1"]
        offset_coords = offset(all_seqs[1],coords)
        ept_annot=vcat(ept_annot,coords[1])
        if length(collect(offset_coords)) > 1
            ept_annot=vcat(ept_annot,zeros(Int,length(collect(offset_coords))-2))
            ept_annot=vcat(ept_annot,coords[end])  
        end  
        ept_annot=vcat(ept_annot,[0])  
        seqs= (seqs .* "|") .* ( (s -> s[offset_coords]).(all_seqs) ) 
    end
    ept_annot=vcat([0],ept_annot)
    seqs= (seqs .* "|")
    seqs[2:end] = (x->replace(x,"-"=>GAP)).(seqs[2:end])

    seps = join([ ch=='|' ? '|' : ' ' for ch in seqs[2] ])
    X_sep = reshape(reduce(hcat, onehot.([seps])), length(aas), :, length([seps]));
    
    Xr = reshape(reduce(hcat, onehot.(seqs[1:1])), length(aas), :, length(seqs[1:1]));
    pr = reshape(mean(Xr; dims=3), size(Xr, 1), size(Xr, 2))
    # H = sum(-xlogx.(pc) / log(2); dims=1)
    # cons = pc .* (log2(5) .- H)
    cons = pr  # use this for frequency only
    sites = [ SequenceLogoSite( (x->WeightedLetter(x...)).(zip(collect(aas),cons[:,i])) )  
                for i in 1:size(cons)[2] ]
    ref_logo = SequenceLogo(sites);
    
    Xc = reshape(reduce(hcat, onehot.(seqs[2:2])), length(aas), :, length(seqs[2:2]));
    pc = reshape(mean(Xc; dims=3), size(Xc, 1), size(Xc, 2))
    # H = sum(-xlogx.(pc) / log(2); dims=1)
    # cons = pc .* (log2(5) .- H)
    cons = pc  # use this for frequency only
    sites = [ SequenceLogoSite( (x->WeightedLetter(x...)).(zip(collect(aas),cons[:,i])) )  
                for i in 1:size(cons)[2] ]
    cons_logo = SequenceLogo(sites);
    
    rest_nams=all_nams[3:end]  # use 3:end for ellpaca
    rest_seqs=seqs[3:end]   # use 3:end for ellpaca
    fvsc=0
    visits=union((x->split(x,"_")[2][1:4]).(rest_nams))
    fvsc=sum([ string(split(nm,"_")[2][1:4]) == string(visits[1]) for nm in rest_nams ])
        # visits=union((x->x[12:12]).(rest_nams)) # for ellpaca
        # fvsc=sum([ string(nm[12:12]) == visits[1] for nm in rest_nams ])
    @show visits
    @show fvsc
    
    p=nothing
    freq_logos=[]
    freq_annots=[]
    for visit in visits
        sel_ind = [ ]
        sel_ind = [ split(nm,"_")[2][1:4] == visit for nm in rest_nams ]
        # sel_ind = [ string(nm[12:12]) == visit for nm in rest_nams ]
        lvsc = sum(sel_ind)
        X = (x->max(x,0)).(reshape(reduce(hcat, onehot.(rest_seqs[sel_ind])), 
                length(aas), :, length(rest_seqs[sel_ind])) .+ X_sep);   # .- Xc 
        p = (x->max(x,0)).(reshape(mean(X; dims=3), size(X, 1), size(X, 2)) .- pc )
        # H = sum(-xlogx.(p) / log(2); dims=1)
        # cons = p .* (log2(5) .- H)
        cons = p  # use this for scaled frequency only
        sites = [ SequenceLogoSite( (x->WeightedLetter(x...)).(zip(collect(aas),cons[:,i])) )  
                    for i in 1:size(cons)[2] ]
        freq_logo = SequenceLogo(sites);
        push!(freq_logos,freq_logo)
        freq_annot=" Variants at $(visit) ($(lvsc))";
        push!(freq_annots, freq_annot)
    end

    path_svg="$(outdir)$(donor)_$(ept_name)_escape_logo.svg"
        
    title="$(donor) escape at $(ept_name) epitope"
       
    pf=luxor_sequence_logo_aa(ref_logo,cons_logo,reverse(freq_logos),
                aa_color_fun,size(cons)[2],20,path_svg,
                annot=ept_annot,
                title=title, 
                # ref_annot="$(all_nams[1][1:15])",
                bot_annot=" Consensus at $(visits[1]) ($(fvsc))",
                freq_annots=reverse(freq_annots));
            
    return(pf)
end

################# functions for generating alignments ##########################

"""
    mafft(inpath, outpath; path="", flags::Vector{String}=String[], kwargs...)
Julia wrapper for mafft.
"""
function my_mafft(inpath, outpath)
    cmd = `mafft-fftns --quiet --thread 2 --ep 3 --op 5 --localpair --out $outpath $inpath`
    println(cmd)
    run(cmd)
end

function consensus(seqs)
    cons = join([mode([seqs[i][j]
                    for i in 1:length(seqs)])
                        for j in 1:length(seqs[1])])
    return(cons)
end

function agreement(ref,seq)
    if length(ref) == length(seq)
        return( sum(collect(ref).==collect(seq)) / length(ref) )
    else
        return 0.0
    end
end
        
function my_translate_to_aa(s::String)
    s=s[1:3*div(length(s),3)]
    rna = convert(LongRNA{4}, LongDNA{4}(s))
    return string(BioSequences.translate(rna))
end

        
function pad_description(description)
    desc=description
    while length(desc)<12
        desc="_"*desc
    end
    return(desc)
end
            
function ungap(s::String)
    return replace(s,"-"=>"")
end
            
function my_write_fasta(filename, seqs;
    names=String[], LongSequence = false, append = false, aa = false)
    if !LongSequence
        if aa
            seqs = [BioSequences.LongAA(s) for s in seqs]
        else
            seqs = [BioSequences.LongDNA{4}(s) for s in seqs]
        end
    end
    stream = open(FASTA.Writer, filename, append=append)
    i = 0
    if length(names) != length(seqs)
        names = [string("seq_", i) for i in 1:length(seqs)]
    end
    for (s, n) in zip(seqs, names)
        i += 1
        write(stream, FASTA.Record(n, s))
    end
    close(stream)
end

#### code to massage CAP alignments with tp annotations ####

function tp_annotate(in_file, out_file)
    # in_file="CAP188alignment/CAP188_functional_aa_ali.fasta"
    # out_file="CAP188alignment/aligned_CAP188_2000_4300.fasta"
    records = collect(FASTX.FASTA.Reader(open(in_file)))
    all_seqs = (x->FASTX.sequence(String,x)).(records)
    all_nams = String.(FASTX.identifier.(records))
    visits = sort(union((x->x[8:11]).(all_nams[3:end])))
    tp_count=0
    for visit in visits
        sel_inds = (x->x[8:11]==visit).(all_nams)
        println(visit, " -> ", sum(sel_inds))
        global tp_count += 1
        for i in 1:length(sel_inds)
            sel_inds[i] ? all_nams[i] = all_nams[i] * "_tp$(tp_count)" : nothing
        end
    end
    my_write_fasta(out_file,all_seqs,names=all_nams,aa=true)
    return
end


##################### set the epitope coordinartes ####################

# default epitope
ept_VRC01 = Dict()
ept_VRC01[1]=197:198
ept_VRC01[2]=230:230
ept_VRC01[3]=276:276
ept_VRC01[4]=278:282
ept_VRC01[5]=365:371
ept_VRC01[6]=427:428
ept_VRC01[7]=430:430
ept_VRC01[8]=455:463
ept_VRC01[9]=465:465
ept_VRC01[10]=467:467
ept_VRC01[11]=469:469
ept_VRC01[12]=471:474

# CAP256 epitope
ept_CAP256 = Dict()
ept_CAP256[1] = 156:163
ept_CAP256[2] = 166:167
ept_CAP256[3] = 169:170
ept_CAP256[4] = 178:179
ept_CAP256[5] = 181:184

# 3L6:
ept_3L6 = Dict()
ept_3L6[1]=262:262
ept_3L6[2]=295:297
ept_3L6[3]=301:302
ept_3L6[4]=323:323
ept_3L6[5]=330:330
ept_3L6[6]=332:332
ept_3L6[7]=439:442
ept_3L6[8]=444:446

# 3D14_1E7:
ept_3D14_1E7=Dict()
ept_3D14_1E7[1]=63:65
ept_3D14_1E7[2]=134:134
ept_3D14_1E7[3]=156:156
ept_3D14_1E7[4]=163:163
ept_3D14_1E7[5]=168:173
ept_3D14_1E7[6]=185:185
ept_3D14_1E7[7]=192:194
ept_3D14_1E7[8]=197:198
ept_3D14_1E7[9]=206:207
ept_3D14_1E7[10]=300:308
ept_3D14_1E7[11]=318:318
ept_3D14_1E7[12]=321:325
ept_3D14_1E7[13]=368:368
ept_3D14_1E7[14]=425:425
ept_3D14_1E7[15]=428:430
ept_3D14_1E7[16]=438:441

epitopes = Dict()
epitopes["VRC01"]=ept_VRC01
epitopes["CAP256"]=ept_CAP256
epitopes["3L6"]=ept_3L6
epitopes["3D14_1E7"]=ept_3D14_1E7

epitope_colors = Dict()
epitope_colors["VRC01"]="lightgreen"
epitope_colors["CAP256"]="coral"
epitope_colors["3L6"]="thistle"
epitope_colors["3D14_1E7"]="lightcyan"

v_loops = Dict()
v_loops[1] = 130:156
v_loops[2] = 157:196
v_loops[3] = 296:331
v_loops[4] = 385:418
v_loops[5] = 460:470

v_loop_colors = Dict()
v_loop_colors[1] = "pink"
v_loop_colors[2] = "bisque"
v_loop_colors[3] = "aquamarine"
v_loop_colors[4] = "navajowhite"
v_loop_colors[5] = "tan"

###################### MaxEscatpe starts here ######################

using Luxor, Random

function drawmatrix(AA::Matrix, V::Vector, donor::String; fs=16,
            epitope_names=[], epitope_colors=Dict(), epitope_coords=Dict(),
            v_loop_colors=Dict(), v_loop_coords=Dict(),
            cellsize = (fs, fs) )

    nrow=size(AA)[1]
    ncol=size(AA)[2]
    chunks = vcat(collect(1:180:ncol),[ncol+1])
    ref_ind=0
    col_ind=0
    background("white")
    tab11=Point(0,0)
    
    for i in 2:length(chunks)
        A = AA[:,chunks[i-1]:chunks[i]-1]
        table = Table(size(A)..., cellsize...)
        
        if i == 2
            tab11=table[1,1]
            Luxor.translate(-table[1,1]+Point(6*fs,6*fs))
            gsave()
            s=1
            ll=sum(length.(epitope_names))+length(epitope_names)
            legend = Table(1,ll,cellsize[1],cellsize[2])
            Luxor.translate(table[1,s] - legend[1,s] + Point(2*fs,-3*fs))
            
            for epn in epitope_names
                sethue(epitope_colors[epn])
                for j in 1:length(epn)
                    box(legend, 1, s+j-1, action = :fill)
                end
                sethue("black")
                fontsize(fs)
                text(epn,legend[1,s],halign=:left,valign=:middle)
                s=s+length(epn)+1
            end
            grestore()
            
            gsave()
            s=1
            n_loops=length(keys(v_loop_colors))
            ll=5*n_loops
            legend = Table(1,ll,cellsize[1],cellsize[2])
            Luxor.translate(table[1,s+40] - legend[1,s] + Point(2*fs,-3*fs))
            
            for i in 1:n_loops
                sethue(v_loop_colors[i])
                for j in 1:4
                    box(legend, 1, s+j-1, action = :fill)
                end
                sethue("black")
                fontsize(fs)
                text("V$(i)",legend[1,s+1],halign=:left,valign=:middle)
                s=s+4+1
            end
            grestore()
        else
            Luxor.translate(Point(0.0,(nrow+4)*fs)+(tab11-table[1,1]))
        end

        fontsize(fs)
        sethue("black")
        text(donor, table[1,1] - Point(fs,fs),
            halign=:right,
            valign=:middle)
        fontsize(2*fs/3)
        for j in 1:length(V)
            if j==1
                sethue("black")
            elseif j==2
                sethue("red")
            elseif j==3
                sethue("green")
            else
                sethue("blue")
            end
            text(V[j], table[j,1] - Point(fs,0.0),
                halign=:right,
                valign=:middle)
                # angle=-π/2)
        end
            
        for i in CartesianIndices(A)
            r, c = Tuple(i)
            
            sethue("white")
            for ept_name in epitope_names
                if r>2 && col_ind + c in epitope_coords[ept_name]
                    sethue(epitope_colors[ept_name])
                end
            end
            for loop in keys(v_loop_coords)
                if r<3 && col_ind + c in v_loop_coords[loop]
                    sethue(v_loop_colors[loop])
                end
            end
            box(table, r, c, action = :fill)
            
            if r==1
                if A[r+1,c]!='-'
                    ref_ind += 1
                    sethue("black")
                    fontsize(fs/2)
                    text(string(ref_ind), table[r,c],
                        halign=:centre,
                        valign=:middle)
                    sethue("gray80")
                    box(table, r, c, action = :stroke)
                else
                    sethue("black")
                    text("-", table[r,c],
                        halign=:centre,
                        valign=:middle)
                    sethue("gray80")
                    box(table, r, c, action = :stroke)
                end
            else
                if r==2
                    sethue("red")
                elseif r==3
                    sethue("green")
                else
                    sethue("blue")
                end
                fontsize(fs)
                text(string(A[r, c]), table[r, c],
                    halign=:centre,
                    valign=:middle)
                sethue("gray80")
                box(table, r, c, action = :stroke)
            end
        end
        col_ind += size(A)[2]
    end
end

function escape_plot(in_file, out_dir)

    # read the input file
    records = collect(FASTX.FASTA.Reader(open(in_file)))
    all_seqs = (x->FASTX.sequence(String,x)).(records)
    all_nams = String.(FASTX.identifier.(records))
    
    ref_seq=all_seqs[1]
    
    # first convert the reference into glyfs
    glyfs = vcat(collect(repeat(' ',length(all_seqs[1]))), collect(all_seqs[1]))
    
    # get donor and visit ids
    donor = string(split(basename(in_file),"_")[2])
    @show donor
    visits=sort(union((x->split(x,"_")[3][1:4]).(all_nams[3:end])))
    # first_inds = (x->split(x,"_")[3][1:4]==visits[1]).(all_nams[3:end])
    # consensus_of_first_visit=consensus(all_seqs[3:end][first_inds])
    consensus_of_first_visit=all_seqs[2]
    
    exp_fv = collect(consensus_of_first_visit)
    glyfs = vcat(glyfs, exp_fv)

    for i in 1:length(visits)
        inds = (x->split(x,"_")[3][1:4]==visits[i]).(all_nams[3:end])
        exp_nv = collect( consensus(all_seqs[3:end][inds]) )
        sames = ( exp_fv .== exp_nv )
        exp_nv[sames] .= ' '
        glyfs = vcat(glyfs,exp_nv)
    end
    
    rows=length(visits) + 3
    cols=length(consensus_of_first_visit)
    # @show rows, cols, rows * cols, length(glyfs)
    
    A = reshape(glyfs, cols, rows)
    A = permutedims(A, [2, 1])
    
    mf = 12
 
    V=vcat(["co-ords","HXB2","Cons C"],visits)
    
    ept_names=[]
    ept_coords=Dict()
    
    for ept_name in keys(epitopes)
        ept = epitopes[ept_name]
        ec = vcat( collect.([ offset(ref_seq,ept[key]) for key in keys(ept) ])... )
        ept_coords[ept_name] = ec
        push!(ept_names,ept_name)
    end
    
    v_loop_coords=Dict()
    for loop in keys(v_loops)
        ept = v_loops[loop]
        ec = collect(offset(ref_seq,v_loops[loop]))
        v_loop_coords[loop] = ec
    end
    
    # height = mf * size(A)[1] + 2*mf*8
    # width = mf * size(A)[2] + 2*mf*16
    # Drawing(width, height, out_dir * donor * "_escape.svg")
    Drawing("A1landscape", out_dir * donor * "_escape.pdf")

    # Luxor.translate(width, height)
    fontsize(8)
    setline(0.5)
    sethue("white")
    drawmatrix(A, V, donor, fs=mf,
            epitope_names = ept_names,
            epitope_colors = epitope_colors,
            epitope_coords = ept_coords,
            v_loop_colors = v_loop_colors,
            v_loop_coords = v_loop_coords)

    finish()
    return
end

################################# script to generate escape map of entire env molecule #####

in_dir = "alignments/"
if length(ARGS) > 0
    in_dir=ARGS[1]
    endswith(in_dir,"/") ? nothing : in_dir=in_dir*"/"
end
if ! isdir(in_dir)
    println("ERROR: alignment directory, $(in_dir) does not exist")
    exit()
end
filepaths = readdir(in_dir)
filepaths = in_dir .* filepaths[(x->endswith(x,".fasta")).(filepaths)]
if length(filepaths) == 0
    println("ERROR: no fasta files in $(in_dir)")
    exit()
end

out_dir="escape_charts/"
if length(ARGS) > 1
    out_dir=ARGS[2]
    endswith(out_dir,"/") ? nothing : out_dir=out_dir*"/"
end
mkpath(out_dir)

for in_file in filepaths[1:end]
    escape_plot(in_file, out_dir)
end

exit()


