
using Printf
using DataStructures  # for DefaultDict, or use Dict with get()

function analyze_otus_file(filepath::String)
    otu_count = Dict{String, Int}()
    read_count = Dict{String, Int}()

    open(filepath) do f
        for line in eachline(f)
            if startswith(line, ">")
                m = match(r">(barcode\d+)seq_[^;]+;size=(\d+)", line)
                if m !== nothing
                    bc   = m.captures[1]
                    size = parse(Int, m.captures[2])
                    otu_count[bc]  = get(otu_count, bc, 0) + 1
                    read_count[bc] = get(read_count, bc, 0) + size
                end
            end
        end
    end

    return otu_count, read_count
end

function print_summary(filepath::String, otu_count::Dict, read_count::Dict)
    println("\nFile: $filepath")
    println(@sprintf("%-12s %6s %12s", "Barcode", "OTUs", "Total Reads"))
    println("-"^32)
    for bc in sort(collect(keys(otu_count)))
        println(@sprintf("%-12s %6d %12d", bc, otu_count[bc], read_count[bc]))
    end
    println(@sprintf("Totals: %d barcodes | %d OTUs | %d reads",
        length(otu_count),
        sum(values(otu_count)),
        sum(values(read_count))))
end

# Find all _OTUs.fa files recursively under purc_poolA_ccs
root = "output/purc_poolA_ccs"
root = "output/purc_poolB_ccs"
for (dirpath, _, files) in walkdir(root)
    for file in files
        if endswith(file, "_OTUs.fa")
            filepath = joinpath(dirpath, file)
            otu_count, read_count = analyze_otus_file(filepath)
            print_summary(filepath, otu_count, read_count)
        end
    end
end