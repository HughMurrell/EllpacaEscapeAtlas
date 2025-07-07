using Pkg
Pkg.activate("./")
Pkg.instantiate()

using CSV, FASTX, DataFrames, Dates, BioSequences
using StatsBase
using Statistics: mean
using MAFFT_jll, FastTree_jll
using Phylo, Plots
# using PhyloNetworks

################# functions for generating alignments ##########################

function my_mafft(inpath, outpath)
    cmd = `mafft-fftns --thread 2 --ep 3 --op 5 --localpair --out $outpath $inpath`
    run(cmd)
    # mafft-fftns() do exe
    #     run(`$exe --thread 2 --ep 3 --op 5 --localpair --out $outpath $inpath`)
    # end
end

function newickTree(seq_file, tree_file; nu=false)
    cmd=`$(fasttree_double()) -nosupport -out $(tree_file) $(seq_file)`
    if nu
        cmd=`$(fasttree_double()) -quiet -nt -gtr -nosupport -out $(tree_file) $(seq_file)`
    end
    run(cmd);
    println("tree generated in $(tree_file)")
    return true
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

function get_max_depth(node)
    if PhyloNetworks.isleaf(node)
        return getparentedge(node).length
    else
        return maximum( getparentedge(node).length .+ get_max_depth.(PhyloNetworks.getchildren(node)) )
    end
end

function ladderize!(net,node)
    if PhyloNetworks.isleaf(node)
        return
    else
        children=PhyloNetworks.getchildren(node)
        if (length(children) == 2)
            depths = get_max_depth.(children)
            if depths[1] > depths[2]
                println("rotating... at $(node.number)")
                PhyloNetworks.rotate!(net,node.number)
            end
        else
            depths = get_max_depth.(children)
            println("node $(node.number) has $(length(depths)) children at depths $(depths)")
            println(children)
        end
        (x->ladderize!(net,x)).(children)
        return
    end
end

# reroot needs PhyloNetworks, try to do this with Phylo
function reroot(tree_file, rooted_file)
    tree = readTopology(read(tree_file,String))
    cons_name="consensus"
    rootatnode!(tree, cons_name)
    directEdges!(tree)
    # cladewiseorder!(tree)
    # ladderize!(tree,PhyloNetworks.getroot(tree))
    writeTopology(tree, rooted_file)
    return(tree)
end

function my_sort!(tree::T; rev = false) where {T <: AbstractTree}
    function loc!(clade::String)
        if Phylo.isleaf(tree, clade)
            node = getnode(tree, clade)
            return getheight(tree,node)
        end

        sizes = map(loc!, Phylo.getchildren(tree, clade))
        node = getnode(tree, clade)
        if T <: LinkTree
            node.other .= node.other[sortperm(sizes, rev = rev)]
        elseif T <: RecursiveTree
            node.conns .= node.conns[sortperm(sizes, rev = rev)]
        end
        return maximum(sizes) # + getheight(tree,node)
    end

    loc!(first(nodenamefilter(isroot, tree)))
    return tree
end



################################# script to align all consensuses ##########################

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

out_dir="alignment_of_consensuses/"
if length(ARGS) > 1
    out_dir=ARGS[2]
    endswith(out_dir,"/") ? nothing : out_dir=out_dir*"/"
end
mkpath(out_dir)

plots_dir="atlas_location_plots/"
if length(ARGS) > 2
    out_dir=ARGS[3]
    endswith(out_dir,"/") ? nothing : out_dir=out_dir*"/"
end
mkpath(plots_dir)

out_file = out_dir * "ellpaca_consensuses.fasta"
out_file_ali = out_dir * "ellpaca_consensuses_aligned.fasta"
fast_tree_file = out_dir * "ellpaca_consensuses.tree"
# ladder_tree_file = out_dir * "ellpaca_consensus_tree.tree"

my_write_fasta(out_file,[],names=[],aa=true,append=false)
ref_written = true  # set this to false if you want to include the HXB2 reference

if ! isfile(out_file_ali)
    count=0
    for in_file in filepaths[1:end]
        global count+=1
        println(count, " -> ", in_file)
        
        # read the input file
        records = collect(FASTX.FASTA.Reader(open(in_file)))
        all_seqs = (x->FASTX.sequence(String,x)).(records)
        all_nams = String.(FASTX.identifier.(records))
        
        if ! ref_written
            my_write_fasta(out_file,[ungap(all_seqs[1])],names=[all_nams[1]],aa=true,append=true)
            global ref_written = true
        end
    
        # get donor and visit ids
        donor = basename(in_file)[1:6]
        @show donor
        visits=sort(union((x->split(x,"_")[2][1:4]).(all_nams[3:end])))
        for i in 1:length(visits)
            @show i, visits[i]
            inds = (x->split(x,"_")[2][1:4]==visits[i]).(all_nams[3:end])
            cons = consensus(all_seqs[3:end][inds])
            my_write_fasta(out_file,[ungap(cons)],names=[donor*"_"*visits[i]],aa=true,append=true)
        end
    end

    my_mafft(out_file,out_file_ali)
    
    # now generate a consensus of consensuses and append to the alignment
    # cons_records = collect(FASTX.FASTA.Reader(open(out_file_ali)))
    # cons_seqs = (x->FASTX.sequence(String,x)).(cons_records)
    # cons_nams = String.(FASTX.identifier.(cons_records))
    
    # cons_seq=consensus(cons_seqs)
    # my_write_fasta(out_file_ali,[cons_seq],names=["consensus"],aa=true,append=true)
      
    # now create the tree
    newickTree(out_file_ali, fast_tree_file)
elseif ! isfile(fast_tree_file)
    println(" ****** alignment exists, skipping mafft, running fasttree *********")
    newickTree(out_file_ali, fast_tree_file)
else
    println(" ****** alignment and tree exist, skipping mafft and fasttree, drawing trees *********")
end



# use PhyloNetworks to reroot and ladderize tree
# reroot(fast_tree_file, fast_tree_file)

ellpaca = open(parsenewick, fast_tree_file);
my_sort!(ellpaca)

# ladderize!(ellpaca,Phylo.getroot(ellpaca))

donors=sort(union((d->d[1:6]).((x->x.name).(getleaves(ellpaca)))))
@show length(donors)

# leaves=getleaves(ellpaca)
# @show leaves
# for leaf in leaves
#     if ! occursin(donor,leaf.name)
#         @show leaf
#         # leaf.name=""
#     else
#         println(" *********** $(leaf.name) is ok ***********************")
#     end
# end

for donor in donors

@show donor

evolve(tree) = map_depthfirst((val, node) -> getheight(tree,node),
                0.0 , tree, Float64) # getheight(tree,node)
trait = evolve(ellpaca)

sizes=(x->occursin(donor,x) ? 25 : 0).(getnodenames(ellpaca,preorder));

# node_colors= (x->occursin(donor,x) ? :red : :grey).(getnodenames(ellpaca,preorder))
# node_colors=reshape(node_colors,1,length(node_colors))

pl=plot(ellpaca, size = (4200, 2970), showtips=true, treetype = :fan,
                            # line_z = trait,
                            # linecolor = :greys,
                            # linewidth = 3,
                            # colorbar = :none,
                            markersize = sizes,
                            markercolor = :green,
                            markerstrokecolor = :green,
                            # markercolor = node_colors,
                            # markerstrokecolor = node_colors,
                            alpha = 0.5,
                            tipfont = (4, :blue),
                            title="\n\n\n\n"*donor, titlefont=font(64,"Computer Modern"))
                                                       
# plot!(pl,bottom_margin = 1Plots.cm, top_margin=5Plots.cm)
    
savefig(pl,"$(plots_dir)$(donor)_fanplot.pdf")

end

exit()



net = readnewick("(A:3.9,(B:1.3,(C:0.5,D:0.8):0.6):2.0);");
@show net
PhyloNetworks.rotate!(net, -4)
@show net
ladderize!(net,PhyloNetworks.getroot(net))
@show net

plot!(pl, ellpaca, size = (4200, 4200), showtips=true, treetype = :fan,
                            line_z = trait,
                            linecolor = :Greys,
                            linewidth = 3,
                            colorbar=:none,
                            tipfont = (4,:blue),
                            markersize=sizes, markercolor = :grey, markerstrokecolor = :grey,
                            title="\n"*donor, titlefont=font(64,"Computer Modern"))
