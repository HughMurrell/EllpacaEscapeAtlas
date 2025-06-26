using Pkg
Pkg.activate("./")
Pkg.instantiate()
using BioSequences, FASTX, YAML, CSV, DataFrames
using StatsBase
using MAFFT_jll, FastTree_jll
using Plots, RecipesBase
using PhyloNetworks

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

function getNameAnnotation(ept,ref_ali,nam_ali)
    ref_nam = map(first,ref_ali)[1]
    ref_seq = map(last, ref_ali)[1]
    all_seqs = map(last,nam_ali)
    all_nams = map(first,nam_ali)
    ept_keys = sort(collect(keys(ept)))
    seqs=["" for seq in all_seqs]
    ept_annot=[]
    for key in ept_keys
        coords=ept[key] #loop_coords["V1"]
        offset_coords = offset(ref_seq,coords)
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
    # seqs[2:end] = (x->replace(x,"-"=>GAP)).(seqs[2:end])

    seps = join([ ch=='|' ? '|' : ' ' for ch in seqs[2] ])
    max_nam_len = maximum( length.(all_nams) )
    anot_dic = Dict()
    for i in 1:length(all_nams)
        if i > 1
            inds = ( collect(seqs[1]) .== collect(seqs[i]) ) .& ( collect(seqs[1]) .!= '|')
            seqs[i] = join( [ inds[j] ? " " : seqs[i][j] for j in 1:length(inds) ])
        end
        anot_dic[all_nams[i]] = all_nams[i] * repeat(" ", max_nam_len - length(all_nams[i]) + 1) * seqs[i]
    end
    return anot_dic, ept_annot
end

function consensus(seqs)
    cons = join([mode([seqs[i][j]
                    for i in 1:length(seqs)])
                        for j in 1:length(seqs[1])])
    return(cons)
end

function dropRefSequenceAndCollapseByVisit(in_file, out_file; nu=false, ept=nothing)
    stream = open(in_file)
    records = collect(FASTX.FASTA.Reader(stream))
    all_seqs = (x->FASTX.sequence(String,x)).(records)
    all_nams = String.(FASTX.description.(records))
    all_nams[2] = "consensus"
    ret_seqs=[]
    ret_nams=[]
    close(stream)
    ref_seq=all_seqs[1]
    ref_nam=all_nams[1]
    
    # recompute consensus of first visit and store it
    visits=sort(union((x->split(x,"_")[2][1:4]).(all_nams[3:end])))
    first_inds = (x->split(x,"_")[2][1:4]==visits[1]).(all_nams[3:end])
    consensus_of_first=consensus(all_seqs[3:end][first_inds])
    all_nams[2]="consensus"
    all_seqs[2]=consensus_of_first
    
    stream = open(FASTA.Writer, out_file, append=false)
    cons_seq=all_seqs[2]
    cons_nam=all_nams[2]
    write(stream, FASTA.Record(cons_nam, cons_seq))
    push!(ret_seqs,cons_seq)
    push!(ret_nams,cons_nam)
    all_nams=all_nams[3:end]
    all_seqs=all_seqs[3:end]
    donor=all_nams[1][1:6]
    if ! isnothing(ept)
        keep_cols=[]
        ept_keys = sort(collect(keys(ept)))
        for key in ept_keys
            coords=ept[key] #loop_coords["V1"]
            keep_cols=vcat(keep_cols, collect(offset(ref_seq,coords)))
        end
        for i in 1:length(all_seqs)
            keep=collect(cons_seq)
            keep[keep_cols]=collect(all_seqs[i])[keep_cols]
            all_seqs[i]=join(keep)
        end
    end
    global visits = sort(unique( (x->x[8:11]).(all_nams) ))
    @show visits
    for visit in visits[1:end]
        visit_seqs=all_seqs[(x->x[8:11]==visit).(all_nams)]
        seq_map=countmap(visit_seqs)
        count=0
        for (key,val) in seq_map
            # write_fasta(out_dir*out_file,[key],names=["$(donor)_$(visit)_$(count)_$(val)"],append=true)
            count+=1
            nam="$(donor)_$(visit)_$(count)_$(val)"
            write(stream, FASTA.Record(nam, key))
            push!(ret_seqs,key)
            push!(ret_nams,nam)
        end
    end
    close(stream)
    return ( zip(ret_nams,ret_seqs), zip([ref_nam],[ref_seq]) )
end

function newickTree(seq_file, tree_file; nu=false)
    cmd=`$(fasttree_double()) -quiet -nosupport -out $(tree_file) $(seq_file)`
    if nu
        cmd=`$(fasttree_double()) -quiet -nt -gtr -nosupport -out $(tree_file) $(seq_file)`
    end
    run(cmd);
    println("tree generated in $(tree_file)")
    return true
end

function my_mafft(in_file, out_file)
    arguments = "$(in_file) > $(out_file) "
    run(`$(mafft_fftns()) $arguments`)
    return
end

# reroot needs PhyloNetworks, try to do this with Phylo
function reroot(tree_file, rooted_file)
    tree = readTopology(read(tree_file,String))
    cons_name="consensus"
    rootatnode!(tree, cons_name)
    directEdges!(tree)
    cladewiseorder!(tree)
    # tree=ladderize!(tree)
    writeTopology(tree, rooted_file) 
    return(tree)
end
        
function getDistanceMatrixFromTree(tree_file)
    tree = readTopology(read(tree_file,String))
    tips = tipLabels(tree)
    dm = pairwiseTaxonDistanceMatrix(tree)
    return(tips, dm)
end
        
function hamming(seq1, seq2)
    if length(seq1) != length(seq2)
        @show "cant compute hamming dist on seqs of different lengths"
    end
    dist=0
    for i in 1:length(seq1)
        if (seq1[i]!='-' || seq2[i]!='-')    
            if (seq1[i]!='N' && seq2[i]!='N')
                if (seq1[i]!=seq2[i])
                    dist+=1
                end
            end
        end
    end
    return dist
    # return sum(collect(seq1) .!= collect(seq2))
end

function get_max_depth(node)
    if isleaf(node)
        return getparentedge(node).length
    else
        return maximum( getparentedge(node).length .+ get_max_depth.(getchildren(node)) )
    end
end

function ladderize!(net,node)
    if isleaf(node)
        return 
    else
        children=getchildren(node)
        if (length(children) == 2)
            depths = get_max_depth.(children)
            if depths[1] > depths[2]
                # println("rotating... at $(node.number)")
                PhyloNetworks.rotate!(net,node.number)
            end
        end
        (x->ladderize!(net,x)).(children)
        return
    end
end

function my_split(s,c,p)
    v=split(s,c)
    if length(v) >= p
        return v[p]
    else
        return v[1]
    end
end

function count_split(s,c)
    v=split(s,c)
    if v[1] == "consensus"
        return "1"
    else
        return v[end]
    end
end

function my_string(i)
    if i > 0
        return string(i)
    end
    return ""
end

function insert_digit(ept_st,ept_anot,d)
    start=findfirst("|",ept_st)[1]
    anots=(x->my_string(x[1])).(ept_anot)
    ret_st=[]
    for i in 1:length(ept_st)
        if i<start || length(anots[i-start+1]) < d
            push!(ret_st, ept_st[i])
        else
            push!(ret_st, anots[i-start+1][d])
        end
    end
    return replace(join(ret_st),"."=>" ")
end

function my_parse(s;default=1)
    ret = default
    try 
        ret = Int(ceil(log2(parse(Int,s))))
    catch e
    end
    return ret
end


function my_color(i)
    all_colors = distinguishable_colors(20)
    if ( isnothing(i) || i > 17 ) 
        return(all_colors[1]) 
    end 
    return(all_colors[i+3])
end

visits = ["consensus"]

visit_color(visit) = my_color(findfirst(visits.==visit))
    
function getnodeheight(node,tree)
    if isrootof(node,tree)
        return 0.0
    else
        return getnodeheight(getparent(node),tree) + getparentedge(node).length
    end
end

function depth_first_heights(node,height,leaf_num)
    if isleaf(node)
        leaf_num += 1
        height[node.number]=leaf_num
        return leaf_num
    else
        children=getchildren(node)
        for child in children
            leaf_num=depth_first_heights(child,height,leaf_num)
        end
        height[node.number]=sum((x->height[x.number]).(children))/length(children)
        return leaf_num
    end
end
    
function findxy(tree::HybridNetwork)
    height, depth, names = Dict(), Dict(), Dict()
    heights=getnodeheights(tree,true)
    i=0
    for n in tree.node
        i=i+1
        height[n.number] = getnodeheight(n,tree)
        names[n.number] = n.name
    end
    leaf_num=0
    depth_first_heights(getroot(tree),depth,leaf_num)
    return height, depth, names
end

@recipe function f(tree::HybridNetwork, 
                    nam_ali::Base.Iterators.Zip{Tuple{Vector{Any}, Vector{Any}}},
                    ref_ali::Base.Iterators.Zip{Tuple{Vector{String}, Vector{String}}},
                    ept::Dict{Any,Any};
                    treetype = :dendrogramhighlighter,
                    marker_group = nothing, line_group = nothing,
                    showtips = true, showtipmarkers = true, tipfont = (8,"Courier Bold"))
    
    linecolor --> :black
    grid --> false
    framestyle --> :none
    legend --> false
    colorbar --> true
    # size --> (1000, 1000)

    d, h, n = findxy(tree)
    adj = 0.1 * maximum(values(d))
    adj <= 0.0 ? adj=0.1 : nothing
    @show adj
    # tipannotations = map(x -> (d[x] + adj, h[x], x), getleafnames(tree))
    tip_numbers = (x->x.number).(tree.leaf)
    tip_names = (x->x.name).(tree.leaf)
    tipmarker_x = map(x -> d[x], tip_numbers)
    tipmarker_y = map(x -> h[x], tip_numbers)
    tipmarker_s = map(x -> count_split(x,'_'), tiplabels(tree))
    tipmarker_g = map(x -> my_split(x,'_',2), tiplabels(tree))
    # global visits = vcat(["consensus"],sort(union(tipmarker_g))[1:end-1])
    max_x = maximum(values(d))
    min_y = minimum(values(h))
    anot_dic, ept_anot = getNameAnnotation(ept, ref_ali, nam_ali)
    tipannotations = map(x -> (max_x + adj, h[x], anot_dic[n[x]]), tip_numbers)
    fw=adj # 0.001
    ept_pieces=split(anot_dic["consensus"],"|")
    ept_st=join( (s->repeat(".",length(s))).(ept_pieces), "|" )
    ept_st_0=replace(ept_st, "."=>"—")
    ept_st_1=insert_digit(ept_st,ept_anot,1)
    ept_st_2=insert_digit(ept_st,ept_anot,2)
    ept_st_3=insert_digit(ept_st,ept_anot,3)
    # ept_st_4=insert_digit(ept_st,ept_anot,4)
    eptannotations = [ (max_x + adj, min_y-1, ept_st_0),
                        (max_x + adj, min_y-2, ept_st_1),
                        (max_x + adj, min_y-3, ept_st_2),
                        (max_x + adj, min_y-4, ept_st_3),
                        (max_x + adj, min_y-5, ept_st_0)]
    # eptannotations = map(i -> (max_x + adj + i*fw, max_y-5, my_string(ept_anot[i])), 1:length(ept_anot))

    # lz = get(plotattributes, :line_z, nothing)
    # mz = get(plotattributes, :marker_z, nothing)
    # isnothing(lz) || (line_z := _handlez(lz, tree, n))
    # isnothing(mz) || (marker_z := _handlez(mz, tree, n))
    # mg = Phylo._handlez(marker_group, tree, n)
    # lg = Phylo._handlez(line_group, tree, n)

    x, y = Float64[], Float64[]
    
    for node in tree.node
        try p = getparent(node) 
            push!(x, d[p.number], d[p.number],
                  d[node.number], NaN)
            push!(y, h[p.number], h[node.number], h[node.number], NaN)
        catch e
        end
    end

    marker_x, marker_y = values(d), values(h) #Phylo._handlemarkers(plotattributes, mg, tree, d, h, n)

    size --> ( length(last(tipannotations[1]))*14 , (6+length(tip_numbers))*16 )
    
    if treetype == :dendrogramhighlighter
        DendrogramHighlighter(x, y, tipannotations, eptannotations, marker_x, marker_y, 
                    tipmarker_x, tipmarker_y, tipmarker_s, tipmarker_g,
                    showtips, showtipmarkers, tipfont)
    elseif treetype == :fanhighlighter
        FanHighlighter(x, y, tipannotations, eptannotations, marker_x, marker_y, 
                    tipmarker_x, tipmarker_y, tipmarker_s, tipmarker_g,
                    showtips, showtipmarkers, tipfont)
    else
        throw(ArgumentError("Unsupported `treetype`; valid values are `:dendrogramhighlighter` or `:fanhighlighter`"))
    end
end

struct DendrogramHighlighter
    x::Any
    y::Any
    tipannotations::Any
    eptannotations::Any
    marker_x::Any
    marker_y::Any
    tipmarker_x::Any
    tipmarker_y::Any
    tipmarker_s::Any
    tipmarker_g::Any
    showtips::Any
    showtipmarkers::Any
    tipfont::Any
    # marker_group::Any
    # line_group::Any
end
struct FanHighlighter
    x::Any
    y::Any
    tipannotations::Any
    eptannotations::Any
    tipmarkers::Any
    marker_x::Any
    marker_y::Any
    tipmarker_x::Any
    tipmarker_y::Any
    tipmarker_s::Any
    tipmarker_g::Any
    showtips::Any
    showtipmarkers::Any
    tipfont::Any
    # marker_group::Any
    # line_group::Any
end

@recipe function f(dend::DendrogramHighlighter)
    ex = extrema(filter(isfinite, dend.x))
    xlims --> (ex[1] - 0.1 * ex[2], ex[2] * 6.0)
        
    ey = extrema(filter(isfinite, dend.y))
    ylims --> (ey[1] - 5.0, ey[2] + 5.0)

    # tip annotations
        
    dend.showtips &&
        (annotations := vcat( map(x -> (x[1], x[2], 
                            text(x[3], visit_color(my_split(x[3],"_",2)), :left, dend.tipfont...)),
                            dend.tipannotations), 
                                map(x -> (x[1], x[2], 
                            text(x[3], :black, :left,  dend.tipfont...)),
                            dend.eptannotations) ) )

    sa = get(plotattributes, :series_annotations, nothing)
        
    @series begin
        seriestype := :path
        markersize := 0
        markershape := :none
        series_annotations := nothing
        label --> ""
        #    primary := false

        # lc = Phylo._extend(get(plotattributes, :linecolor, nothing), dend.x)
        # lc !== nothing && (linecolor := lc)
        # la = Phylo._extend(get(plotattributes, :linealpha, nothing), dend.x)
        # la !== nothing && (linealpha := la)
        # lz = Phylo._extend(get(plotattributes, :line_z, nothing), dend.x)
        # lz !== nothing && (line_z := lz)

        dend.x, dend.y
    end
#    if !isempty(dend.marker_x) || sa !== nothing
#        if isnothing(dend.marker_group)
#            @series begin
#                seriestype := :scatter
#                sa !== nothing && (series_annotations := sa)
#                label --> ""
#                dend.marker_x, dend.marker_y
#            end
#        else
#            groups = sort(unique(dend.marker_group))
#            for group in groups
#                idxs = findall(==(group), dend.marker_group)
#                @series begin
#                    seriestype := :scatter
#                    sa !== nothing && (series_annotations := sa[idxs])
#                    label --> string(group)
#                    dend.marker_x[idxs], dend.marker_y[idxs]
#                end
#            end
#        end
#    end
    if dend.showtipmarkers 
        @series begin
            seriestype := :scatter
            label --> ""
            markersize := (x->3+my_parse(x)).(dend.tipmarker_s)
            markercolor := (x->visit_color(x)).(dend.tipmarker_g)
            markerstrokecolor := (x->:steelblue).(dend.tipmarker_s)
            dend.tipmarker_x, dend.tipmarker_y
        end
    end
            
    primary := false
    label := ""
    return nothing
end

####################### script starts now ###################################

in_dir = "alignments/"
in_files = readdir(in_dir)
in_files=in_files[(x->endswith(x,".fasta")).(in_files)]
in_paths = in_dir .* in_files

work_dir = "working/"
mkpath(work_dir)
out_dir = "highlighter_plots/"

for in_path in in_paths[1:end]
    donor = basename(in_path)[1:6]
    println("$(donor).......")
    col_path = work_dir * "$(donor)_aa_alignment_noref_collapsed.fasta"
    tree_file = work_dir * "$(donor).tree"
    rerooted_tree_file = work_dir * "$(donor)_rerooted.tree"
    plot_file = out_dir * "$(donor)_tree_highlighter.svg"
    ept=ept_CAP256
    ept_name="CAP256"
    title=""
    if isnothing(ept)
        global out_dir = "$(ept_name)_full_highlighter_plots/"
        mkpath(out_dir)
        title = "$(donor) full highlighter plot \n at $(ept_name) epitope"
        plot_file = out_dir * "$(donor)_FULL_tree_highlighter.svg"
    else
        global out_dir = "$(ept_name)_induced_highlighter_plots/"
        mkpath(out_dir)
        title="$(donor) induced highlighter plot \n at $(ept_name) epitope"
        plot_file = out_dir * "$(donor)_$(ept_name)_tree_highlighter.svg"
    end
    nam_ali, ref_ali = dropRefSequenceAndCollapseByVisit(in_path,col_path,ept=ept)
    ept=ept_CAP256
    newickTree(col_path,tree_file )
    tree=reroot(tree_file,rerooted_tree_file)
    PhyloNetworks.resetnodenumbers!(tree; checkpreorder=true, type=:postorder)
    ladderize!(tree,getroot(tree))
    # global visits = sort(union((x->my_split(x,"_",2)).(tiplabels(tree)[1:end])))
    pl=plot(tree, nam_ali, ref_ali, ept_CAP256, linecolor = :orange, linewidth = 3,
        showtips = true, showtipmarkers=true, treetype = :dendrogramhighlighter, aligntips=true,
        title=title)
    savefig(pl,plot_file)
    println("highlighter plot saved in $(plot_file)")
end
println("ALL DONE, exiting...")

exit()

