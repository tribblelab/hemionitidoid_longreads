#!/usr/bin/env julia

# Script to add reference sequences to OTU files
# For each locus, extracts matching references from reference.fasta
# and appends them to the OTU file, creating *_OTUs_wref.fa

using Printf

# Configuration
REFERENCE_FASTA = joinpath("input", "namedrefsrun", "reference.fasta")

# Check if reference file exists
isfile(REFERENCE_FASTA)

# Locus and pool lists
LOCI = ["APPEFP", "CRY2", "GAPCPSH", "IBR3", "INFA", "MATK", "PGIC", "PSBM", "RPL2", "SQD1"]
POOLS = ["purc_poolA_ccs", "purc_poolB_ccs"]

# Locus name normalization mapping
LOCUS_MAPPING = Dict(
    "GAPCPSH" => "gapCpSh",
    "INFA" => "infA",
    "MATK" => "matK",
    "PGIC" => "pgiC",
    "PSBM" => "psbM",
    "RPL2" => "rpl2"
)

# Read reference FASTA file
function read_reference_fasta(filename)
    """Read FASTA file and return dict of locus => [(header, sequence), ...]"""
    refs = Dict{String, Vector{Tuple{String, String}}}()
    
    current_header = ""
    current_seq = ""
    current_locus = ""
    
    open(filename, "r") do f
        for line in eachline(f)
            line = strip(line)
            
            if startswith(line, ">")
                # Save previous sequence if exists
                if current_header != "" && current_locus != ""
                    if !haskey(refs, current_locus)
                        refs[current_locus] = []
                    end
                    push!(refs[current_locus], (current_header, current_seq))
                end
                
                # Parse new header: >locus=APPEFP/group=H/ref_taxon=Aleuritopteris
                current_header = line
                current_seq = ""
                
                m = match(r"locus=([^/]+)/", line)
                if m !== nothing
                    current_locus = m.captures[1]
                end
            else
                current_seq *= line
            end
        end
        
        # Save last sequence
        if current_header != "" && current_locus != ""
            if !haskey(refs, current_locus)
                refs[current_locus] = []
            end
            push!(refs[current_locus], (current_header, current_seq))
        end
    end
    
    return refs
end

# Read all references
all_refs = read_reference_fasta(REFERENCE_FASTA)

# Process each pool and locus
total_processed = 0
total_skipped = 0

for pool in POOLS
    println("Processing $pool...")
    
    for locus in LOCI
        # Normalize locus name for matching
        locus_key = get(LOCUS_MAPPING, locus, locus)
        
        # Define file paths
        otu_file = joinpath("output", pool, locus, "$(locus)_OTUs.fa")
        output_file = joinpath("output", pool, locus, "$(locus)_OTUs_wref.fa")
        
        # Check if OTU file exists
        if !isfile(otu_file)
            println("  ⊗ $locus: OTU file not found, skipping...")
            total_skipped += 1
            continue
        end
        
        # Check if we have references for this locus
        if !haskey(all_refs, locus_key)
            println("  ⚠ $locus: No references found for locus '$locus_key', copying file without refs...")
            cp(otu_file, output_file, force=true)
            total_processed += 1
            continue
        end
        
        # Count OTU sequences
        otu_count = 0
        open(otu_file, "r") do f
            for line in eachline(f)
                if startswith(line, ">")
                    otu_count += 1
                end
            end
        end
        
        # Write combined file: OTUs + references
        ref_count = length(all_refs[locus_key])
        
        open(output_file, "w") do out
            # First, copy all OTU sequences
            open(otu_file, "r") do f
                for line in eachline(f)
                    println(out, line)
                end
            end
            
            # Then append reference sequences
            for (header, sequence) in all_refs[locus_key]
                println(out, header)
                println(out, sequence)
            end
        end
        
        @printf("  ✓ %-10s: %3d OTUs + %2d refs → %s\n", 
                locus, otu_count, ref_count, basename(output_file))
        total_processed += 1
    end
    
    println()
end
