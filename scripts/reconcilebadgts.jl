using PhyloNetworks

"""
    BarcodeTable
 
Holds the parsed barcode CSV as a vector of named tuples for fast lookup.
Load once with `load_barcode_table(csv_path)` and reuse.
"""
struct BarcodeTable
    rows::Vector{Dict{String,String}}
    gene_cols::Vector{String}
end

"""
    load_barcode_table(csv_path::String) -> BarcodeTable
 
Parse the barcode CSV. The CSV has columns:
  pool, matK, pgic, APPEFP, CRY2, gapCpSh, SQD1, rpl2, psbM, infA, IBR3, purc_id
"""
function load_barcode_table(csv_path::String)::BarcodeTable
    lines = readlines(csv_path)
    # Strip BOM if present
    lines[1] = lstrip(lines[1], '\ufeff')
    header = split(lines[1], ",")
    gene_cols = [h for h in header if h ∉ ("pool", "purc_id")]
    rows = Dict{String,String}[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        vals = split(line, ",")
        row = Dict(header[i] => (i <= length(vals) ? vals[i] : "") for i in eachindex(header))
        push!(rows, row)
    end
    return BarcodeTable(rows, gene_cols)
end

# Taxa names follow the pattern: CollectionCode_Genus_species_ClusterN_size_N
# Extract the genus from a tip label
function get_genus(label::String)::String
    parts = split(label, "_")
    # Guard: expect at least 2 parts (code + genus)
    length(parts) >= 2 ? parts[2] : label
end

is_pentagramma(label::String) = get_genus(label) == "Pentagramma"

# Return all tip names in the subtree rooted at node n
function tips_below(n::PhyloNetworks.Node)::Vector{String}
    tips = String[]
    collect_tips!(n, tips)
    return tips
end

function collect_tips!(n::PhyloNetworks.Node, tips::Vector{String})
    if n.leaf
        push!(tips, n.name)
    else
        for e in n.edge
            child = getchild(e)
            child !== n && collect_tips!(child, tips)
        end
    end
end


# Find the MRCA of all Pentagramma tips by finding the deepest internal node
# whose subtree contains ALL Pentagramma tips and no more than necessary.
function find_mrca(penta_tips::Vector{String}, t::HybridNetwork)::Union{PhyloNetworks.Node,Nothing}
    penta_set = Set(penta_tips)
    best = nothing
    best_size = typemax(Int)
    for n in t.node
        n.leaf && continue
        below = tips_below(n)
        # Node must contain all Pentagramma tips
        all(p -> p in Set(below), penta_tips) || continue
        # Keep the one with fewest total descendants (i.e. deepest/smallest clade)
        if length(below) < best_size
            best = n
            best_size = length(below)
        end
    end
    return best
end

"""
    getparent_node(n::PhyloNetworks.Node) -> Union{PhyloNetworks.Node, Nothing}
 
Return the parent node of `n`, or `nothing` if `n` is the root.
"""
function getparent_node(n::PhyloNetworks.Node)::Union{PhyloNetworks.Node,Nothing}
    for e in n.edge
        getchild(e) === n || continue   # this edge points INTO n, so its other end is the parent
        return getparent(e)
    end
    return nothing  # n is the root
end

# Find the MRCA node of all Pentagramma tips, then return a non-Pentagramma
# tip in the sister clade to use as outgroup.
function find_pentagramma_base_tip(t::HybridNetwork)::Union{String,Nothing}
    penta_tips = [n.name for n in t.node if n.leaf && is_pentagramma(n.name)]
    isempty(penta_tips) && return nothing

    # If every tip is Pentagramma there is no meaningful outgroup
    all_tips = [n.name for n in t.node if n.leaf]
    length(penta_tips) == length(all_tips) && return nothing

    mrca = find_mrca(penta_tips, t)
    mrca === nothing && return nothing

    # Walk up from mrca via its parent edge, then find a non-Pentagramma
    # tip in the sister clade
    for e in mrca.edge
        getchild(e) === mrca || continue   # only the parent edge of mrca
        parent = getparent(e)
        for e2 in parent.edge
            getchild(e2) === parent && continue  # skip parent's own parent edge
            sibling = getchild(e2)
            sibling === mrca && continue
            sister_tips = tips_below(sibling)
            non_penta = filter(l -> !is_pentagramma(l), sister_tips)
            isempty(non_penta) || return non_penta[1]
        end
    end
    return nothing
end

"""
    find_outgroup_clade_edge(t::HybridNetwork, genera::Vector{String})
        -> Union{PhyloNetworks.Edge, Nothing}

Find the edge to root on such that ALL tips belonging to any of `genera`
fall on one side (the outgroup side). Returns the edge whose child node
subtends exactly the outgroup clade — or the largest such subset if the
genera are paraphyletic with Pentagramma interspersed.

Strategy: among all internal edges whose child subtree contains at least
one tip from `genera` and NO Pentagramma tips, pick the one with the
largest child subtree (most inclusive pure-outgroup clade). If no such
edge exists (outgroup tips are mixed with Pentagramma), fall back to the
edge whose child subtree has the most outgroup tips regardless.
"""
function find_outgroup_clade_edge(t::HybridNetwork,
    genera::Vector{String})::Union{PhyloNetworks.Edge,Nothing}
    is_target(label) = get_genus(label) ∈ genera

    best_pure_edge = nothing
    best_pure_size = 0
    best_mixed_edge = nothing
    best_mixed_count = 0

    for e in t.edge
        child = getchild(e)
        below = tips_below(child)
        isempty(below) && continue

        n_target = count(is_target, below)
        n_target == 0 && continue

        if !any(is_pentagramma, below)
            # Pure outgroup subtree — prefer largest
            if n_target > best_pure_size
                best_pure_edge = e
                best_pure_size = n_target
            end
        else
            # Mixed — track most outgroup tips as fallback
            if n_target > best_mixed_count
                best_mixed_edge = e
                best_mixed_count = n_target
            end
        end
    end

    return best_pure_edge !== nothing ? best_pure_edge : best_mixed_edge
end

"""
    root_gene_tree!(t::HybridNetwork) -> Symbol

Root `t` in-place using outgroup priority. For each priority group, finds
the edge subtending the most inclusive pure-outgroup clade and calls
rootonedge!, placing the full clade as outgroup rather than a single tip.

Priority:
  1. Ynesmexia
  2. Myriopteris
  3. Pellaea or Paragymnopteris
  4. Notholaena or Cheiloplecton
  5. Base of the Pentagramma clade (tip in sister group — single-tip fallback)

Returns a symbol indicating which strategy was used.
"""
function root_gene_tree!(t::HybridNetwork)::Symbol
    tip_names = [n.name for n in t.node if n.leaf]
    non_penta = filter(l -> !is_pentagramma(l), tip_names)

    if isempty(non_penta)
        @info "Tree is all Pentagramma — no outgroup rooting applied."
        return :all_pentagramma
    end

    priority_groups = [
        (["Ynesmexia"], :ynesmexia),
        (["Myriopteris"], :myriopteris),
        (["Pellaea", "Paragymnopteris"], :pellaea_paragymnopteris),
        (["Notholaena", "Cheiloplecton"], :notholaena_cheiloplecton),
        (["Gaga", "Aleuritopteris", "Adiantopsis"], :otherhems)
    ]

    for (genera, strategy) in priority_groups
        e = find_outgroup_clade_edge(t, genera)
        if e !== nothing
            PhyloNetworks.rootonedge!(t, e)
            clade_tips = tips_below(getchild(e))
            @info "Rooted on edge subtending $(length(clade_tips))-tip clade (strategy: $strategy): $(join(strip_tip_suffixes.(clade_tips), ", "))"
            return strategy
        end
    end

    # Priority 5: single-tip fallback at base of Pentagramma clade
    tip = find_pentagramma_base_tip(t)
    if tip !== nothing
        PhyloNetworks.rootatnode!(t, tip)
        @info "Rooted at base of Pentagramma clade via tip $tip"
        return :pentagramma_base
    end

    @warn "Could not determine a root — tree unchanged."
    return :unrooted
end

"""
    penta_tips_external_to_main_clade(t::HybridNetwork; min_clade_size::Int=10)
        -> Vector{String}

Find Pentagramma tips that sit in a grade external to the largest purely-
Pentagramma clade.

Strategy:
1. Find the largest purely-Pentagramma internal node (the 'main' clade).
2. Find the MRCA of the main clade node and the nearest non-Pentagramma tip —
   i.e. the deepest node whose subtree contains main_node AND at least one
   outgroup tip. This is the 'boundary' node where the grade begins.
   (Walking up via tips_below fails when rootatnode! placed an outgroup leaf
   as the tree root, since tips_below of that root leaf returns only itself.)
3. Return all Pentagramma tips in the boundary subtree that are NOT in the
   main clade.

`min_clade_size` guards against trivially small pure-Pentagramma nodes
being mistaken for the main clade (default 10).
"""
function penta_tips_external_to_main_clade(t::HybridNetwork;
    min_clade_size::Int=10)::Vector{String}
    # Step 1: find the largest purely-Pentagramma internal node
    main_node = nothing
    main_size = 0
    for n in t.node
        n.leaf && continue
        below = tips_below(n)
        isempty(below) && continue
        all(is_pentagramma, below) || continue
        if length(below) > main_size
            main_node = n
            main_size = length(below)
        end
    end

    if main_node === nothing || main_size < min_clade_size
        @warn "No large purely-Pentagramma clade found (largest = $main_size tips)"
        return String[]
    end
    @info "Main Pentagramma clade: $main_size tips"

    main_tips = Set(tips_below(main_node))

    # Step 2: collect ALL Pentagramma tips in the tree that are not in the main clade.
    # We do not restrict to a boundary subtree — the debug output showed that the
    # grade tips can land on either side of the root depending on which edge was
    # rooted on, so any subtree-based boundary will miss tips on the root's far side.
    # The correct definition is simply: every Pentagramma tip that is not part of
    # the main clade is a candidate for pool C reconciliation.
    all_penta = [n.name for n in t.node if n.leaf && is_pentagramma(n.name)]
    external = filter(tip -> tip ∉ main_tips, all_penta)
    return external
end

# ---------------------------------------------------------------------------
# Convenience wrapper: process a single tree file
# ---------------------------------------------------------------------------
"""
    process_gene_tree(path, pool, bt=nothing; genus_map=nothing, min_clade_size=10)
 
Read a newick gene tree, root it, find any Pentagramma tips sitting in a
grade external to the main (largest pure-Pentagramma) clade, and for each
such tip check whether there is a pool C barcode match.
 
When a barcode match is found **and** `genus_map` is supplied, the function
also:
1. Looks up the pool C individual's purc_id via `find_poolC_match`.
2. Formats a `reference.fasta`-ready entry via `format_reference_entry`,
   using the pool C purc_id as the label.
3. Collects all formatted entries in the fifth return value.
 
- `path`           : path to the .tre file
- `pool`           : `"A"` or `"B"` (which pool this gene tree came from)
- `bt`             : optional `BarcodeTable`; if provided, barcode reconciliation is run
- `genus_map`      : optional `Dict{String,String}` from `load_genus_map`; enables
                     reference formatting on barcode matches
- `min_clade_size` : passed to `penta_tips_external_to_main_clade` (default 10)
 
Returns `(tree, rooting_strategy, external_tips, tip_matches, ref_entries)` where
`tip_matches` is a `Vector{NamedTuple}` with fields
`(penta_tip, poolC_id, ref_entry)` — one entry per external Pentagramma tip —
and `ref_entries` is a flat `Vector{String}` of all formatted reference.fasta
entries (empty when `genus_map` is `nothing` or no pool C match was found).
"""
function process_gene_tree(path::String, pool::String,
    bt::Union{BarcodeTable,Nothing}=nothing;
    genus_map::Union{Dict{String,String},Nothing}=nothing,
    min_clade_size::Int=10)
    t = readTopology(path)
    strategy = root_gene_tree!(t)
    gene = splitext(basename(path))[1]   # e.g. "APPEFP" from "APPEFP.tre"

    external_tips = penta_tips_external_to_main_clade(t; min_clade_size=min_clade_size)
    tip_matches = NamedTuple[]
    ref_entries = String[]

    if isempty(external_tips)
        println("No Pentagramma tips external to main clade found.")
    else
        println("$(length(external_tips)) Pentagramma tip(s) external to main clade:")
        for ptip in external_tips
            poolC_id = bt !== nothing ?
                       find_poolC_match(ptip, pool, gene, bt) : nothing

            ref_entry = nothing
            if genus_map !== nothing && poolC_id !== nothing
                # Fetch the Pentagramma tip's own sequence and label it with
                # the pool C individual's ID
                fetched = fetch_nonpenta_sequence(ptip, pool, gene)
                if fetched !== nothing
                    _, seq = fetched
                    ref_entry = format_reference_entry(poolC_id, seq, gene, genus_map)
                    push!(ref_entries, ref_entry)
                end
            end

            push!(tip_matches, (
                penta_tip=ptip,
                poolC_id=poolC_id,
                ref_entry=ref_entry,
            ))

            if bt !== nothing
                match_str = poolC_id !== nothing ? "✓ pool C: $poolC_id" : "✗ no pool C match"
                println("  $match_str  |  $(strip_tip_suffixes(ptip))")
            else
                println("  $(strip_tip_suffixes(ptip))")
            end
        end
    end

    return (t, strategy, external_tips, tip_matches, ref_entries)
end

# ---------------------------------------------------------------------------
# Batch: process all gene trees in a directory
# ---------------------------------------------------------------------------
"""
    process_gene_tree_dir(dir::String, pool::String; ext=".tre", bt=nothing, genus_map=nothing)
 
Recursively process all gene tree files under `dir`.
 
- `pool`           : `"A"` or `"B"` — which sequencing pool these trees came from
- `bt`             : optional `BarcodeTable` for barcode reconciliation
- `genus_map`      : optional `Dict{String,String}` from `load_genus_map`; when supplied,
                     sequences are fetched and formatted for confirmed barcode matches
- `ext`            : file extension to match (default `.tre`)
- `min_clade_size` : passed to `penta_tips_external_to_main_clade` (default 10)
 
Returns a `Dict` mapping filename => `(tree, rooting_strategy, external_tips, tip_matches, ref_entries)`.
"""
function process_gene_tree_dir(dir::String, pool::String;
    ext::String=".tre",
    bt::Union{BarcodeTable,Nothing}=nothing,
    genus_map::Union{Dict{String,String},Nothing}=nothing,
    min_clade_size::Int=10,
    results::Dict{String,Any}=Dict{String,Any}())
    for entry in readdir(dir, join=true)
        if isdir(entry)
            process_gene_tree_dir(entry, pool; ext=ext, bt=bt, genus_map=genus_map,
                min_clade_size=min_clade_size, results=results)
        elseif endswith(entry, ext)
            println("\n=== $(basename(entry)) ===")
            results[basename(entry)] = process_gene_tree(entry, pool, bt;
                genus_map=genus_map,
                min_clade_size=min_clade_size)
        end
    end
    return results
end

"""
    rename_treefile_to_tre(dir::String)
 
Recursively rename all `.treefile` files in `dir` and its subdirectories
to `.tre`. Skips any file if a `.tre` with the same stem already exists.
Returns a list of renamed file paths.
"""
function rename_treefile_to_tre(dir::String, renamed::Vector{String}=String[])::Vector{String}
    files = readdir(dir)
    for f in files
        fullpath = joinpath(dir, f)
        if isdir(fullpath)
            rename_treefile_to_tre(fullpath, renamed)
        elseif endswith(f, ".treefile")
            newpath = fullpath[1:end-9] * ".tre"
            if isfile(newpath)
                @warn "Skipping $fullpath — $(basename(newpath)) already exists"
                continue
            end
            mv(fullpath, newpath)
            push!(renamed, newpath)
        end
    end
    println("Renamed $(length(renamed)) file(s) under $dir")
    return renamed
end

# ---------------------------------------------------------------------------
# Barcode lookup: find pool C individual sharing barcode with a tip
# ---------------------------------------------------------------------------

# Map gene filename stem (case-insensitive) to the CSV column name
const GENE_TO_COL = Dict(
    "appefp" => "APPEFP",
    "ibr3" => "IBR3",
    "sqd1" => "SQD1",
    "cry2" => "CRY2",
    "matk" => "matK",
    "pgic" => "pgic",
    "gapcpsh" => "gapCpSh",
    "rpl2" => "rpl2",
    "psbm" => "psbM",
    "infa" => "infA",
)

"""
    strip_tip_suffixes(tip::String) -> String
 
Remove `_ClusterN_size_N` suffixes from a tip label, returning just the
collection-code prefix e.g. `KMW013_Pentagramma_triangularis`.
"""
function strip_tip_suffixes(tip::String)::String
    # Match _Cluster<digits>_size_<digits> at the end, case-insensitive
    m = match(r"^(.+?)_[Cc]luster\d+_size_\d+$", tip)
    m !== nothing ? string(m.captures[1]) : tip
end

"""
    collection_code(tip::String) -> String
 
Extract the leading collection code from a tip label, e.g.
`KMW013_Pentagramma_triangularis` → `KMW013`.
"""
function collection_code(tip::String)::String
    split(strip_tip_suffixes(tip), "_")[1]
end

"""
    find_poolC_match(tip::String, pool::String, gene::String,
                     bt::BarcodeTable) -> Union{String, Nothing}
 
Given a tip label from a pool A or pool B gene tree, find the pool C
individual that shares the same barcode for that gene.
 
Arguments:
- `tip`  : full tip label e.g. `KMW013_Pentagramma_triangularis_Cluster132_size_2841`
- `pool` : `"A"` or `"B"` (the pool the gene tree came from)
- `gene` : gene name matching the tree filename stem, e.g. `"APPEFP"` or `"IBR3"`
- `bt`   : `BarcodeTable` loaded with `load_barcode_table`
 
Returns the pool C `purc_id` string, or `nothing` if no match is found.
"""
function find_poolC_match(tip::String, pool::String, gene::String,
    bt::BarcodeTable)::Union{String,Nothing}
    pool = uppercase(pool)
    pool ∈ ("A", "B") || error("pool must be \"A\" or \"B\", got \"$pool\"")

    col = get(GENE_TO_COL, lowercase(gene), nothing)
    if col === nothing
        @warn "Unknown gene name: $gene. Known genes: $(keys(GENE_TO_COL))"
        return nothing
    end
    col ∈ bt.gene_cols || (@warn "Column $col not found in barcode table"; return nothing)

    code = collection_code(tip)

    # Find the row for this tip in the correct pool
    src_row = findfirst(r -> r["pool"] == pool && r["purc_id"] != "" &&
                                 startswith(r["purc_id"], code), bt.rows)
    if src_row === nothing
        @warn "No $pool row found for collection code $code (tip: $tip)"
        return nothing
    end

    barcode = bt.rows[src_row][col]
    if isempty(barcode)
        @warn "No barcode for gene $gene in row for $code (pool $pool)"
        return nothing
    end

    # Find pool C row with same barcode for this gene
    c_row = findfirst(r -> r["pool"] == "C" && r[col] == barcode, bt.rows)
    if c_row === nothing
        # Barcodes >= 44 have no pool C entry by design — suppress the warning
        bnum_m = match(r"_(\d+)", barcode)
        bnum = bnum_m !== nothing ? parse(Int, bnum_m.captures[1]) : 0
        bnum < 44 && @warn "No pool C individual found with barcode '$barcode' for gene $gene"
        return nothing
    end

    return bt.rows[c_row]["purc_id"]
end

# ---------------------------------------------------------------------------
# Sequence retrieval from pool/gene OTU alignment files
# ---------------------------------------------------------------------------

"""
    OTU_DIR
 
Root directory that contains the gene-tree pool directories.
Structure: `<OTU_DIR>/pool<X>/<GENE>/GENE_OTUs.fa`
"""
const OTU_DIR = joinpath("output", "genetrees")

"""
    fetch_nonpenta_sequence(tip::String, pool::String, gene::String;
                            otu_dir=OTU_DIR, prefer_largest=true)
    -> Union{Tuple{String,String}, Nothing}
 
Fetch a single sequence for the given tip (Pentagramma or otherwise) from
the appropriate pool/gene `GENE_OTUs.fa` alignment file.
 
- `tip`            : tip label from a gene tree, e.g.
                     `"KMW096_Pentagramma_maxonii_Cluster32_size_695"`
- `pool`           : `"A"` or `"B"` (which pool the tip came from)
- `gene`           : gene name as it appears in the directory, e.g. `"APPEFP"`
- `prefer_largest` : when multiple clusters exist for the same individual,
                     choose the one with the highest `;size=N` value (default `true`)
 
Only **one** sequence is returned per individual (identified by the leading
collection code). Returns `(header, sequence)` on success, `nothing` if the
file or individual cannot be found.
"""
function fetch_nonpenta_sequence(tip::String, pool::String, gene::String;
    otu_dir::String=OTU_DIR,
    prefer_largest::Bool=true,
    cluster_num::Union{Int,Nothing}=nothing)::Union{Tuple{String,String},Nothing}
    pool = uppercase(pool)
    gene = uppercase(gene)
    pool ∈ ("A", "B") || error("pool must be \"A\" or \"B\", got \"$pool\"")

    otu_file = joinpath(otu_dir, "pool$(pool)", gene, "$(gene)_OTUs.fa")
    if !isfile(otu_file)
        @warn "OTU file not found: $otu_file"
        return nothing
    end

    # Extract the collection code from the tip label (e.g. "KMW096")
    code = collection_code(tip)

    # Parse the FASTA file, collecting all entries whose header starts with code
    candidates = Tuple{String,String,Int}[]   # (header, seq, size)
    current_header = ""
    current_seq = IOBuffer()
    in_target = false

    function flush_candidate!()
        if in_target && current_header != ""
            seq = String(take!(current_seq))
            # Parse ;size=N from the header
            m = match(r";size=(\d+)", current_header)
            sz = m !== nothing ? parse(Int, m.captures[1]) : 0
            push!(candidates, (current_header, seq, sz))
        end
        in_target = false
        current_header = ""
        take!(current_seq)   # discard buffer contents
    end

    open(otu_file, "r") do f
        for line in eachline(f)
            line = strip(line)
            if startswith(line, ">")
                flush_candidate!()
                header = lstrip(line, '>')
                if startswith(header, code * "_")
                    current_header = header
                    in_target = true
                end
            elseif in_target
                write(current_seq, line)
            end
        end
        flush_candidate!()
    end

    isempty(candidates) && (@warn "No sequences found for code $code in $otu_file"; return nothing)

    # Automatically parse cluster details from the tip string if not explicitly passed
    if cluster_num === nothing
        m_c = match(r"_[Cc]luster(\d+)", tip)
        if m_c !== nothing
            cluster_num = parse(Int, m_c.captures[1])
        end
    end

    # Filter candidates by cluster_num if specified/parsed
    filtered = filter(candidates) do c
        c_hdr = c[1]
        if cluster_num !== nothing
            mc = match(r"_[Cc]luster(\d+)", c_hdr)
            c_cnum = mc !== nothing ? parse(Int, mc.captures[1]) : nothing
            c_cnum == cluster_num || return false
        end
        return true
    end

    # If we found matches with the specified filters, pick the best one.
    # Otherwise, fallback to all candidates and warn.
    if !isempty(filtered)
        best = prefer_largest ? argmax(c -> c[3], filtered) : first(filtered)
    else
        if cluster_num !== nothing
            @warn "No exact match for cluster_num=$cluster_num for code $code. Falling back to largest cluster."
        end
        best = prefer_largest ? argmax(c -> c[3], candidates) : first(candidates)
    end

    return (best[1], best[2])
end

# ---------------------------------------------------------------------------
# Format a retrieved sequence as a reference.fasta entry
# ---------------------------------------------------------------------------

"""
    format_reference_entry(header_raw::String, sequence::String,
                           locus::String, genus_map::Dict{String,String};
                           max_header_len::Int=50) -> String
 
Format a sequence for addition to `input/namedrefsrun/reference.fasta`.
 
The FASTA header is constructed as:
    `>locus=<locus>/group=<group>/ref_taxon=<SampleID>`
 
where:
- `<locus>`    — lowercased locus string (e.g. `appefp`)
- `<group>`    — looked up from `genus_map` using the genus parsed from `header_raw`
- `<SampleID>` — the collection code from `header_raw` (e.g. `XZ051`)
 
The header is truncated/fitted so it does not exceed `max_header_len` characters
(including the leading `>`). If the full ref_taxon fits, it is kept; otherwise
only the collection code is used.
 
Arguments:
- `header_raw`  : raw FASTA header (without `>`) from the OTU file,
                  e.g. `"XZ051_Doryopteris_angelica_Cluster3;size=82"`
- `sequence`    : nucleotide string (gaps removed before writing)
- `locus`       : gene name to embed, e.g. `"APPEFP"` (lowercased automatically)
- `genus_map`   : `Dict{String,String}` mapping genus → refgroup,
                  as loaded from `notes/genustorefmap.csv`
- `max_header_len` : maximum total header length including `>` (default 50)
 
Returns a two-line `String` of the form `">header\\nsequence"`.
"""
function format_reference_entry(header_raw::String, sequence::String,
    locus::String, genus_map::Dict{String,String};
    max_header_len::Int=50)::String
    locus_lc = lowercase(locus)

    # Strip cluster/size suffixes then split on underscore
    # e.g. "XZ051_Doryopteris_angelica_Cluster3;size=82"
    #   -> "XZ051_Doryopteris_angelica"
    clean = strip_tip_suffixes(header_raw)
    # strip_tip_suffixes handles "_ClusterN_size_N"; also strip ";size=N" variants
    clean = replace(clean, r";size=\d+$" => "")
    clean = replace(clean, r"_[Cc]luster\d+.*$" => "")

    parts = split(clean, "_")
    sample_code = length(parts) >= 1 ? parts[1] : clean
    genus = length(parts) >= 2 ? parts[2] : ""

    # Lookup group
    group = get(genus_map, genus, "")
    if isempty(group)
        @warn "Genus '$genus' not found in genus map; group will be empty"
    end

    # Build the full ref_taxon label (sample + genus + species if it fits)
    full_taxon = clean   # e.g. "XZ051_Doryopteris_angelica"
    code_only = sample_code  # fallback

    # Construct candidates from longest to shortest ref_taxon
    for ref_taxon in (full_taxon, code_only)
        candidate = ">locus=$(locus_lc)/group=$(group)/ref_taxon=$(ref_taxon)"
        if length(candidate) <= max_header_len
            # Remove gap characters from the sequence
            seq_clean = replace(sequence, r"[-]" => "")
            return "$(candidate)\n$(seq_clean)"
        end
    end

    # Hard truncate as a last resort
    prefix = ">locus=$(locus_lc)/group=$(group)/ref_taxon="
    room = max_header_len - length(prefix)
    ref_taxon_trunc = room > 0 ? code_only[1:min(room, length(code_only))] : ""
    seq_clean = replace(sequence, r"[-]" => "")
    return "$(prefix)$(ref_taxon_trunc)\n$(seq_clean)"
end

"""
    load_genus_map(csv_path::String) -> Dict{String,String}
 
Parse `notes/genustorefmap.csv` into a `Dict` mapping genus → refgroup.
"""
function load_genus_map(csv_path::String)::Dict{String,String}
    d = Dict{String,String}()
    lines = readlines(csv_path)
    # skip header
    for line in lines[2:end]
        isempty(strip(line)) && continue
        cols = split(strip(line), ",")
        length(cols) >= 2 && (d[strip(cols[1])] = strip(cols[2]))
    end
    return d
end

"""
    switch_tip_to_poolC_reference(tip::String, pool::String, gene::String,
                                  bt::BarcodeTable, genus_map::Dict{String,String};
                                  cluster_num::Union{Int, Nothing}=nothing) -> String
 
For a manually specified `tip` (e.g. "KMW013_Pentagramma_triangularis") in a given `pool` ("A" or "B")
and `gene` (e.g. "APPEFP"), find its corresponding pool C individual (sharing the same barcode).
Then retrieve that tip's sequence from the OTU alignment and format it for the reference file,
labeling it with the pool C individual's ID (the "other barcode ID").
 
If the tip name contains cluster details (e.g., "_Cluster132_size_2841"), they will be parsed automatically.
You can also manually specify `cluster_num` as a keyword argument.
 
Throws an error if the barcode number is 44 or greater (as pool C is not present).
"""
function switch_tip_to_poolC_reference(tip::String, pool::String, gene::String,
    bt::BarcodeTable, genus_map::Dict{String,String};
    cluster_num::Union{Int,Nothing}=nothing)::String
    pool = uppercase(pool)
    pool ∈ ("A", "B") || error("pool must be \"A\" or \"B\", got \"$pool\"")

    # Resolve cluster_num from tip label if not explicitly provided
    if cluster_num === nothing
        m_c = match(r"_[Cc]luster(\d+)", tip)
        if m_c !== nothing
            cluster_num = parse(Int, m_c.captures[1])
        else
            error("No cluster_num provided and none found in tip label '$tip'. " *
                  "Pass cluster_num explicitly, e.g. cluster_num=25.")
        end
    end

    col = get(GENE_TO_COL, lowercase(gene), nothing)
    if col === nothing
        error("Unknown gene name: $gene. Known genes: $(keys(GENE_TO_COL))")
    end

    code = collection_code(tip)

    # Find the row for this tip in the correct pool A or B
    src_row = findfirst(r -> r["pool"] == pool && r["purc_id"] != "" &&
                                 startswith(r["purc_id"], code), bt.rows)
    if src_row === nothing
        error("No row found in barcode table for tip code '$code' in pool $pool")
    end

    barcode = bt.rows[src_row][col]
    if isempty(barcode)
        error("No barcode registered for gene $gene in row for '$code' (pool $pool)")
    end

    # Error handling for barcodes 44+
    bnum_m = match(r"_(\d+)", barcode)
    bnum = bnum_m !== nothing ? parse(Int, bnum_m.captures[1]) : 0
    if bnum >= 44
        error("Barcode number is $bnum (>= 44) for '$barcode'. Pool C is not present for barcodes >= 44.")
    end

    # Find the corresponding pool C individual
    c_row = findfirst(r -> r["pool"] == "C" && r[col] == barcode, bt.rows)
    if c_row === nothing
        error("No pool C individual found with barcode '$barcode' for gene $gene")
    end
    poolC_id = bt.rows[c_row]["purc_id"]

    # Fetch the sequence for the tip from the OTU alignment
    fetched = fetch_nonpenta_sequence(tip, pool, gene; cluster_num=cluster_num)
    if fetched === nothing
        error("Could not retrieve OTU sequence for tip '$tip' (code '$code') in pool $pool, gene $gene")
    end
    hdr, seq = fetched

    # Format the reference entry using the pool C ID (the other barcode ID)
    return format_reference_entry(poolC_id, seq, gene, genus_map)
end

# ---------------------------------------------------------------------------
# Script-level execution
# ---------------------------------------------------------------------------

rename_treefile_to_tre("output/genetrees")

bt = load_barcode_table("input/barcode_to_samples.csv")
gmap = load_genus_map("notes/genustorefmap.csv")

results_A = process_gene_tree_dir("output/genetrees/poolA/", "A"; bt=bt, genus_map=gmap)
results_B = process_gene_tree_dir("output/genetrees/poolB/", "B"; bt=bt, genus_map=gmap)

#=
manual_refs = String[
    switch_tip_to_poolC_reference("KMW087_Pentagramma_triangularis", "B", "MATK", bt, gmap; cluster_num=25),
    switch_tip_to_poolC_reference("KMW105_Pentagramma_triangularis", "B", "MATK", bt, gmap; cluster_num=2),

    switch_tip_to_poolC_reference("KMW107_Pentagramma_triangularis", "B", "APPEFP", bt, gmap; cluster_num=1),
    switch_tip_to_poolC_reference("KMW107_Pentagramma_triangularis", "B", "APPEFP", bt, gmap; cluster_num=2),
    switch_tip_to_poolC_reference("KMW096_Pentagramma_maxonii", "B", "APPEFP", bt, gmap; cluster_num=32),
    switch_tip_to_poolC_reference("KMW087_Pentagramma_triangularis", "B", "APPEFP", bt, gmap; cluster_num=5),
    switch_tip_to_poolC_reference("KMW098_Pentagramma_triangularis", "B", "APPEFP", bt, gmap; cluster_num=6),
    switch_tip_to_poolC_reference("KMW098_Pentagramma_triangularis", "B", "APPEFP", bt, gmap; cluster_num=3),
    switch_tip_to_poolC_reference("KMW095_Pentagramma_glanduloviscida", "B", "APPEFP", bt, gmap; cluster_num=10), 

    switch_tip_to_poolC_reference("ES296_Pentagramma_triangularis", "A", "INFA", bt, gmap; cluster_num=216),
    switch_tip_to_poolC_reference("ES299_Pentagramma_semipallida", "A", "INFA", bt, gmap; cluster_num=301), 

    switch_tip_to_poolC_reference("ES296_Pentagramma_triangularis", "A", "PSBM", bt, gmap; cluster_num=161), 

    switch_tip_to_poolC_reference("ES296_Pentagramma_triangularis", "A", "RPL2", bt, gmap; cluster_num=72),
    switch_tip_to_poolC_reference("ES299_Pentagramma_semipallida", "A", "RPL2", bt, gmap; cluster_num=13), 

    switch_tip_to_poolC_reference("ES300_Pentagramma_rebmanii", "A", "SQD1", bt, gmap; cluster_num=11),
    switch_tip_to_poolC_reference("ES300_Pentagramma_rebmanii", "A", "SQD1", bt, gmap; cluster_num=115),
    switch_tip_to_poolC_reference("ES296_Pentagramma_triangularis", "A", "SQD1", bt, gmap; cluster_num=44),
    switch_tip_to_poolC_reference("ES281_Pentagramma_semipallida", "A", "SQD1", bt, gmap; cluster_num=33),
    switch_tip_to_poolC_reference("ES284_Pentagramma_viscosa", "A", "SQD1", bt, gmap; cluster_num=243),
    switch_tip_to_poolC_reference("ES284_Pentagramma_viscosa", "A", "SQD1", bt, gmap; cluster_num=249),
    switch_tip_to_poolC_reference("ES284_Pentagramma_viscosa", "A", "SQD1", bt, gmap; cluster_num=241),
    switch_tip_to_poolC_reference("ES280_Pentagramma_glanduloviscida", "A", "SQD1", bt, gmap; cluster_num=10),
    switch_tip_to_poolC_reference("ES277_Pentagramma_rebmanii", "A", "SQD1", bt, gmap; cluster_num=156), 

    switch_tip_to_poolC_reference("ES307_Pentagramma_triangularis", "A", "CRY2", bt, gmap; cluster_num=186),
    switch_tip_to_poolC_reference("ES307_Pentagramma_triangularis", "A", "CRY2", bt, gmap; cluster_num=142),
    switch_tip_to_poolC_reference("ES307_Pentagramma_triangularis", "A", "CRY2", bt, gmap; cluster_num=410),
    switch_tip_to_poolC_reference("ES307_Pentagramma_triangularis", "A", "CRY2", bt, gmap; cluster_num=348),
    switch_tip_to_poolC_reference("ES300_Pentagramma_rebmanii", "A", "CRY2", bt, gmap; cluster_num=200),
    switch_tip_to_poolC_reference("ES300_Pentagramma_rebmanii", "A", "CRY2", bt, gmap; cluster_num=198),
    switch_tip_to_poolC_reference("ES300_Pentagramma_rebmanii", "A", "CRY2", bt, gmap; cluster_num=205),
    switch_tip_to_poolC_reference("ES296_Pentagramma_triangularis", "A", "CRY2", bt, gmap; cluster_num=6),
    switch_tip_to_poolC_reference("ES298_Pentagramma_triangularis", "A", "CRY2", bt, gmap; cluster_num=9),
]
=#

all_refs = vcat(
    [r[5] for r in values(results_A)]...,
    [r[5] for r in values(results_B)]...,
    #manual_refs,
)

# Append to reference.fasta
REF_FASTA = joinpath("input", "namedrefsrun", "reference.fasta")
open(REF_FASTA, "a") do io
    println(io) # Start with a newline to separate from existing file content
    for entry in all_refs
        println(io, entry)
    end
end

println("\nAppended $(length(all_refs)) reference entr$(length(all_refs) == 1 ? "y" : "ies") to $REF_FASTA")