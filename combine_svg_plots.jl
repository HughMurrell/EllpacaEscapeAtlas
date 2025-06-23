using Pkg
Pkg.activate("./")

using Plots
using EzXML
using Images
using FileIO

"""
Ensure the xlink namespace is present in an SVG file.

Args:
    svg_path: Path to the SVG file
"""
function ensure_xlink_namespace(svg_path::String)
    content = read(svg_path, String)
    if occursin("xlink:href", content) && !occursin("xmlns:xlink", content)
        content = replace(
            content,
            "xmlns=\"http://www.w3.org/2000/svg\"" =>
            "xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\""
        )
        open(svg_path, "w") do f
            write(f, content)
        end
        println("Fixed xlink namespace in $svg_path")
    end
end

"""
Combine two SVG files horizontally into a single SVG file.

Args:
    svg1_path: Path to the first SVG file
    svg2_path: Path to the second SVG file  
    output_path: Path for the output combined SVG file
    separation: Gap between the two plots (default: 50 pixels)
"""
function combine_svg_plots(svg1_path::String, svg2_path::String, output_path::String; separation::Int=50)
    
    # Ensure xlink namespace is present if needed
    ensure_xlink_namespace(svg1_path)
    ensure_xlink_namespace(svg2_path)
    
    # Read and parse the SVG files
    svg1_content = read(svg1_path, String)
    svg2_content = read(svg2_path, String)
    
    # Parse SVG content to extract dimensions
    svg1_doc = parsexml(svg1_content)
    svg2_doc = parsexml(svg2_content)
    
    # Get the root SVG elements
    svg1_root = root(svg1_doc)
    svg2_root = root(svg2_doc)
    
    # Extract dimensions and viewBox
    width1, height1, viewbox1 = extract_svg_dimensions_and_viewbox(svg1_root)
    width2, height2, viewbox2 = extract_svg_dimensions_and_viewbox(svg2_root)
    
    println("SVG 1 dimensions: $width1 x $height1, viewBox: $viewbox1")
    println("SVG 2 dimensions: $width2 x $height2, viewBox: $viewbox2")
    
    # Extract the inner content of both SVGs
    svg1_inner = extract_svg_content(svg1_content)
    svg2_inner = extract_svg_content(svg2_content)
    
    # Calculate combined dimensions for vertical stacking
    combined_width = max(width1, width2)
    combined_height = height1 + height2 + separation

    println("Combined dimensions: $combined_width x $combined_height")
    
    # Create the combined SVG with proper viewBox (vertical stacking)
    combined_svg = """<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<svg xmlns=\"http://www.w3.org/2000/svg\" 
     xmlns:xlink=\"http://www.w3.org/1999/xlink\"
     width=\"$combined_width\" 
     height=\"$combined_height\" 
     viewBox=\"0 0 $combined_width $combined_height\">
    
    <!-- First SVG positioned at the top, centered horizontally -->
    <g transform=\"translate($((combined_width - width1) / 2), 0)\">
        <svg width=\"$width1\" height=\"$height1\" viewBox=\"$viewbox1\">
$svg1_inner
        </svg>
    </g>
    
    <!-- Second SVG positioned below the first, centered horizontally -->
    <g transform=\"translate($((combined_width - width2) / 2), $(height1 + separation))\">
        <svg width=\"$width2\" height=\"$height2\" viewBox=\"$viewbox2\">
$svg2_inner
        </svg>
    </g>
</svg>"""
    
    # Write the combined SVG to file
    write(output_path, combined_svg)
    println("Combined SVG saved to: $output_path")
end

"""
Extract width, height, and viewBox from an SVG element.
Handles both explicit width/height attributes and viewBox.
"""
function extract_svg_dimensions_and_viewbox(svg_element)
    # Try to get width and height from attributes first
    width_attr = svg_element["width"]
    height_attr = svg_element["height"]
    viewbox_attr = svg_element["viewBox"]
    
    if width_attr !== nothing && height_attr !== nothing
        # Remove "px" if present and convert to Float64
        width = parse(Float64, replace(width_attr, "px" => ""))
        height = parse(Float64, replace(height_attr, "px" => ""))
        
        # Use viewBox if available, otherwise use width/height
        if viewbox_attr !== nothing
            viewbox = viewbox_attr
        else
            viewbox = "0 0 $width $height"
        end
        
        return width, height, viewbox
    end
    
    # If no width/height attributes, try viewBox
    if viewbox_attr !== nothing
        viewbox_parts = split(viewbox_attr)
        if length(viewbox_parts) >= 4
            width = parse(Float64, viewbox_parts[3])
            height = parse(Float64, viewbox_parts[4])
            return width, height, viewbox_attr
        end
    end
    
    # Default dimensions if nothing else works
    return 800.0, 600.0, "0 0 800 600"
end

"""
Extract the inner content of an SVG file (everything between <svg> tags).
"""
function extract_svg_content(svg_content::String)
    # Find the opening <svg> tag
    svg_start = findfirst("<svg", svg_content)
    if svg_start === nothing
        error("Could not find <svg> tag in content")
    end
    
    # Find the closing </svg> tag
    svg_end = findlast("</svg>", svg_content)
    if svg_end === nothing
        error("Could not find </svg> tag in content")
    end
    
    # Find the end of the opening <svg> tag (look for the closing >)
    tag_end = findnext('>', svg_content, svg_start[1])
    if tag_end === nothing
        error("Could not find end of <svg> tag")
    end
    
    # Extract content between <svg> and </svg> (exclude the tags themselves)
    inner_content = svg_content[tag_end+1:svg_end[1]-1]
    
    return inner_content
end


# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 3
        println("Usage: julia combine_svg_plots.jl <svg_dir_1> <svg_dir_2> <output_dir> [separation]")
        println("Example: julia combine_svg_plots.jl highlighter_dir logo_dir combined_dir 50")
        
        # Run example with sample files if they exist
        # sample_svg1 = "VRC01_induced_highlighter_plots/CAP437_VRC01_tree_highlighter.svg"
        # sample_svg2 = "VRC01_induced_highlighter_plots/CAP430_VRC01_tree_highlighter.svg"
        
        # if isfile(sample_svg1) && isfile(sample_svg2)
        #     println("\nRunning example with sample files...")
        #     combine_svg_plots(sample_svg1, sample_svg2, "combined_sample.svg", separation=50)
        # end
    else
        svg1_path = ARGS[1]
        svg2_path = ARGS[2]
        output_path = ARGS[3]
        separation = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 50

        files_1 = readdir(svg1_path)
        files_1 = sort( files_1[(x->endswith(x,".svg")).(files_1)] )
        @show length(files_1)
        files_2 = readdir(svg2_path)
        files_2 = sort( files_2[(x->endswith(x,".svg")).(files_2)] )
        @show length(files_2)

        for (file_1, file_2) in zip(files_1, files_2)
            donor_1 = file_1[1:6]
            donor_2 = file_2[1:6]
            if donor_1 == donor_2
                donor = donor_1
                @show donor 
                combine_svg_plots(svg1_path*"/"*file_1, svg2_path*"/"*file_2, output_path*"/"*"$(donor).svg", 
                    separation=separation)
            else
                @error("files do not match")
                exit()
            end
        end
    end
end

