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

# Return the first tip name whose genus matches any of the target genera (in priority order)
function find_outgroup_tip(t::HybridNetwork, priority_genera::Vector{String})::Union{String,Nothing}
    tip_names = [n.name for n in t.node if n.leaf]
    for genus in priority_genera
        match = findfirst(label -> get_genus(label) == genus, tip_names)
        match !== nothing && return tip_names[match]
    end
    return nothing
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
        parent = _get_parent(e)
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
    root_gene_tree!(t::HybridNetwork) -> Symbol
 
Root `t` in-place using outgroup priority:
  1. Ynesmexia
  2. Myriopteris
  3. Pellaea or Paragymnopteris
  4. Notholaena or Cheiloplecton
  5. Base of the Pentagramma clade (tip in sister group)
 
Returns a symbol indicating which strategy was used:
  :ynesmexia, :myriopteris, :pellaea_paragymnopteris,
  :notholaena_cheiloplecton, :pentagramma_base, or :all_pentagramma
"""
function root_gene_tree!(t::HybridNetwork)::Symbol
    tip_names = [n.name for n in t.node if n.leaf]
    non_penta = filter(l -> !is_pentagramma(l), tip_names)

    # If no non-Pentagramma tips, nothing to do
    if isempty(non_penta)
        @info "Tree is all Pentagramma — no outgroup rooting applied."
        return :all_pentagramma
    end

    # Priority 1–4: specific genera (Pellaea and Paragymnopteris share a slot,
    # as do Notholaena and Cheiloplecton)
    priority_groups = [
        ([:Ynesmexia], :ynesmexia),
        ([:Myriopteris], :myriopteris),
        ([:Pellaea, :Paragymnopteris], :pellaea_paragymnopteris),
        ([:Notholaena, :Cheiloplecton], :notholaena_cheiloplecton),
    ]

    for (genera_syms, strategy) in priority_groups
        genera_strs = string.(genera_syms)
        tip = find_outgroup_tip(t, genera_strs)
        if tip !== nothing
            PhyloNetworks.rootatnode!(t, tip)
            @info "Rooted on $tip (strategy: $strategy)"
            return strategy
        end
    end

    # Priority 5: base of Pentagramma clade
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
    mixed_pentagramma_clades(t::HybridNetwork) -> Vector{Vector{String}}
 
After rooting, find every internal node whose descendant tips include BOTH
Pentagramma and non-Pentagramma taxa. Returns a list of tip-label vectors,
one per such clade. Clades are not nested (only the deepest qualifying node
on each path is returned).
"""
function mixed_pentagramma_clades(t::HybridNetwork)::Vector{Vector{String}}
    results = Vector{Vector{String}}()
    root = t.node[t.rooti]
    find_mixed_clades!(root, nothing, results)
    return results
end

function find_mixed_clades!(n::PhyloNetworks.Node, parent::Union{PhyloNetworks.Node,Nothing},
    results::Vector{Vector{String}})
    n.leaf && return

    below = tips_below(n)
    has_penta = any(is_pentagramma, below)
    has_nonpenta = any(!is_pentagramma, below)

    if has_penta && has_nonpenta
        # Only record this node if none of its children already captured the mix
        # (i.e. this is the shallowest mixed node on this path — we recurse into
        # children first and only push if no child already added a superset)
        child_mixed = false
        for e in n.edge
            child = getchild(e)
            child === n && continue
            if !child.leaf
                cb = tips_below(child)
                if any(is_pentagramma, cb) && any(!is_pentagramma, cb)
                    child_mixed = true
                end
            end
        end

        if !child_mixed
            push!(results, below)
        end
    end

    # Recurse into children
    for e in n.edge
        child = getchild(e)
        child !== n && find_mixed_clades!(child, n, results)
    end
end

# ---------------------------------------------------------------------------
# Convenience wrapper: process a single tree file
# ---------------------------------------------------------------------------
"""
    process_gene_tree(path, pool, bt=nothing; genus_map=nothing, max_print=60)

Read a newick gene tree, root it, report any mixed Pentagramma clades,
and for each mixed clade check whether each Pentagramma tip shares its
barcode with any non-Pentagramma tip in the same clade.

When a barcode match is found **and** `genus_map` is supplied, the function
also:
1. Fetches the non-Pentagramma sequence from the appropriate pool/gene
   `GENE_OTUs.fa` file via `fetch_nonpenta_sequence`.
2. Formats a `reference.fasta`-ready entry via `format_reference_entry`.
3. Prints the formatted entry to stdout.
4. Collects all formatted entries in the fifth return value.

- `path`      : path to the .tre file
- `pool`      : `"A"` or `"B"` (which pool this gene tree came from)
- `bt`        : optional `BarcodeTable`; if provided, barcode reconciliation is run
- `genus_map` : optional `Dict{String,String}` from `load_genus_map`; enables
                sequence retrieval and reference formatting on barcode matches
- `max_print` : suppress per-tip output for clades with more than this many
                Pentagramma individuals (default 60); the clade header is always printed

Returns `(tree, rooting_strategy, mixed_clades, barcode_matches, ref_entries)` where
`barcode_matches` is a `Vector{Vector{NamedTuple}}`, one entry per mixed clade,
each entry being a list of `(penta_tip, nonpenta_tip, poolC_id, match, ref_entry)`,
and `ref_entries` is a flat `Vector{String}` of all formatted reference.fasta entries
(empty when `genus_map` is `nothing`).
"""
function process_gene_tree(path::String, pool::String,
    bt::Union{BarcodeTable,Nothing}=nothing;
    genus_map::Union{Dict{String,String},Nothing}=nothing,
    max_print::Int=60)
    t = readTopology(path)
    strategy = root_gene_tree!(t)
    gene = splitext(basename(path))[1]   # e.g. "APPEFP" from "APPEFP.tre"

    clades = mixed_pentagramma_clades(t)
    barcode_matches = Vector{Vector{NamedTuple}}()
    ref_entries = String[]

    if isempty(clades)
        println("No mixed Pentagramma/outgroup clades found.")
    else
        println("$(length(clades)) mixed clade(s) found:")
        for (i, clade) in enumerate(clades)
            penta_tips = filter(is_pentagramma, clade)
            nonpenta_tips = filter(!is_pentagramma, clade)
            verbose = length(penta_tips) <= max_print
            println("  Clade $i ($(length(clade)) tips: $(length(penta_tips)) Pentagramma, $(length(nonpenta_tips)) other)$(verbose ? ":" : " [too large to print — suppressed]:")")

            clade_matches = NamedTuple[]

            for ptip in penta_tips
                # Find the pool C individual for this Pentagramma tip
                poolC_id = bt !== nothing ?
                           find_poolC_match(ptip, pool, gene, bt) : nothing

                for ntip in nonpenta_tips
                    ntip_code = collection_code(ntip)

                    # Exact match: pool C individual for this Pentagramma IS ntip
                    is_match = poolC_id !== nothing &&
                               startswith(poolC_id, ntip_code)

                    # Soft match: pool C genus == ntip genus, different individual
                    soft_match = !is_match && poolC_id !== nothing &&
                                 get_genus(poolC_id) == get_genus(ntip)

                    # Fetch + format on exact or soft match
                    ref_entry = nothing
                    if genus_map !== nothing && (is_match || soft_match)
                        fetched = fetch_nonpenta_sequence(ntip, pool, gene)
                        if fetched !== nothing
                            hdr, seq = fetched
                            # Exact → label with OTU header; soft → label with pool C ID
                            label = is_match ? hdr : poolC_id
                            ref_entry = format_reference_entry(label, seq, gene, genus_map)
                        end
                    end

                    push!(clade_matches, (
                        penta_tip=ptip,
                        nonpenta_tip=ntip,
                        poolC_id=poolC_id,
                        match=is_match,
                        soft_match=soft_match,
                        ref_entry=ref_entry,
                    ))

                    if verbose
                        if bt !== nothing
                            match_str = is_match ? "✓ BARCODE MATCH" :
                                        soft_match ? "~ soft match    " :
                                        "✗ no match      "
                            println("    $match_str  |  Pentagramma: $(strip_tip_suffixes(ptip))")
                            println("               |  Other:       $(strip_tip_suffixes(ntip))")
                            poolC_id !== nothing && println("               |  Pool C ID:   $poolC_id")
                        else
                            println("    Pentagramma: $(strip_tip_suffixes(ptip))")
                            println("    Other:       $(strip_tip_suffixes(ntip))")
                        end
                    end
                end
            end

            # Always print a match summary (and collect all ref_entries)
            n_exact = count(m -> m.match, clade_matches)
            n_soft = count(m -> m.soft_match, clade_matches)
            n_refs = count(m -> m.ref_entry !== nothing, clade_matches)
            if !verbose
                println("    $(n_exact) exact / $(n_soft) soft match(es), $(n_refs) reference entries generated")
            end
            for m in clade_matches
                m.ref_entry !== nothing && push!(ref_entries, m.ref_entry)
            end

            push!(barcode_matches, clade_matches)
        end
    end

    return (t, strategy, clades, barcode_matches, ref_entries)
end

# ---------------------------------------------------------------------------
# Batch: process all gene trees in a directory
# ---------------------------------------------------------------------------
"""
    process_gene_tree_dir(dir::String, pool::String; ext=".tre", bt=nothing, genus_map=nothing)
 
Recursively process all gene tree files under `dir`.

- `pool`      : `"A"` or `"B"` — which sequencing pool these trees came from
- `bt`        : optional `BarcodeTable` for barcode reconciliation
- `genus_map` : optional `Dict{String,String}` from `load_genus_map`; when supplied,
                sequences are fetched and formatted for confirmed barcode matches
- `ext`       : file extension to match (default `.tre`)

Returns a `Dict` mapping filename => `(tree, rooting_strategy, mixed_clades, barcode_matches, ref_entries)`.
"""
function process_gene_tree_dir(dir::String, pool::String;
    ext::String=".tre",
    bt::Union{BarcodeTable,Nothing}=nothing,
    genus_map::Union{Dict{String,String},Nothing}=nothing,
    results::Dict{String,Any}=Dict{String,Any}())
    for entry in readdir(dir, join=true)
        if isdir(entry)
            process_gene_tree_dir(entry, pool; ext=ext, bt=bt, genus_map=genus_map, results=results)
        elseif endswith(entry, ext)
            println("\n=== $(basename(entry)) ===")
            results[basename(entry)] = process_gene_tree(entry, pool, bt; genus_map=genus_map)
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

Root directory that contains the namedrun gene-tree pool directories.
Structure: `<OTU_DIR>/pool<X>/<GENE>/GENE_OTUs.fa`
"""
const OTU_DIR = joinpath("output", "genetrees", "namedrun")

"""
    fetch_nonpenta_sequence(nonpenta_tip::String, pool::String, gene::String;
                            otu_dir=OTU_DIR, prefer_largest=true)
    -> Union{Tuple{String,String}, Nothing}

Given a barcode-match result, fetch a single sequence for the non-Pentagramma
individual from the appropriate pool/gene `GENE_OTUs.fa` alignment file.

- `nonpenta_tip` : tip label from a gene tree, e.g.
                   `"XZ051_Doryopteris_angelica_Cluster3;size=82"`
- `pool`         : `"A"` or `"B"` (which pool the tip came from)
- `gene`         : gene name as it appears in the directory, e.g. `"APPEFP"`
- `prefer_largest` : when multiple clusters exist for the same individual,
                   choose the one with the highest `;size=N` value (default `true`)

Only **one** sequence is returned per individual (identified by the leading
collection code, e.g. `XZ051`). Returns `(header, sequence)` on success,
`nothing` if the file or individual cannot be found.
"""
function fetch_nonpenta_sequence(nonpenta_tip::String, pool::String, gene::String;
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

    # Extract the collection code from the tip label (e.g. "XZ051")
    code = collection_code(nonpenta_tip)

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
        m_c = match(r"_[Cc]luster(\d+)", nonpenta_tip)
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

rename_treefile_to_tre("output/genetrees/namedrun")

bt = load_barcode_table("input/barcode_to_samples.csv")
gmap = load_genus_map("notes/genustorefmap.csv")

results_A = process_gene_tree_dir("output/genetrees/namedrun/poolA/", "A"; bt=bt, genus_map=gmap)
results_B = process_gene_tree_dir("output/genetrees/namedrun/poolB/", "B"; bt=bt, genus_map=gmap)

manual_refs = String[
    switch_tip_to_poolC_reference("KMW087_Pentagramma_triangularis", "B", "MATK", bt, gmap; cluster_num=25),
    switch_tip_to_poolC_reference("KMW105_Pentagramma_triangularis", "B", "MATK", bt, gmap; cluster_num=2),
    #switch_tip_to_poolC_reference("KMW080_Pentagramma_maxonii", "B", "MATK", bt, gmap; cluster_num=4),

    switch_tip_to_poolC_reference("KMW107_Pentagramma_triangularis", "B", "APPEFP", bt, gmap; cluster_num=1),
    switch_tip_to_poolC_reference("KMW107_Pentagramma_triangularis", "B", "APPEFP", bt, gmap; cluster_num=2),
    switch_tip_to_poolC_reference("KMW096_Pentagramma_maxonii", "B", "APPEFP", bt, gmap; cluster_num=32),
    switch_tip_to_poolC_reference("KMW087_Pentagramma_triangularis", "B", "APPEFP", bt, gmap; cluster_num=5),
    switch_tip_to_poolC_reference("KMW098_Pentagramma_triangularis", "B", "APPEFP", bt, gmap; cluster_num=6),
    switch_tip_to_poolC_reference("KMW098_Pentagramma_triangularis", "B", "APPEFP", bt, gmap; cluster_num=3),
    switch_tip_to_poolC_reference("KMW095_Pentagramma_glanduloviscida", "B", "APPEFP", bt, gmap; cluster_num=10), switch_tip_to_poolC_reference("ES296_Pentagramma_triangularis", "A", "INFA", bt, gmap; cluster_num=216),
    switch_tip_to_poolC_reference("ES299_Pentagramma_semipallida", "A", "INFA", bt, gmap; cluster_num=301), switch_tip_to_poolC_reference("ES296_Pentagramma_triangularis", "A", "PSBM", bt, gmap; cluster_num=161), switch_tip_to_poolC_reference("ES296_Pentagramma_triangularis", "A", "RPL2", bt, gmap; cluster_num=72),
    switch_tip_to_poolC_reference("ES299_Pentagramma_semipallida", "A", "RPL2", bt, gmap; cluster_num=13), switch_tip_to_poolC_reference("ES300_Pentagramma_rebmanii", "A", "SQD1", bt, gmap; cluster_num=11),
    switch_tip_to_poolC_reference("ES300_Pentagramma_rebmanii", "A", "SQD1", bt, gmap; cluster_num=115),
    switch_tip_to_poolC_reference("ES296_Pentagramma_triangularis", "A", "SQD1", bt, gmap; cluster_num=44),
    switch_tip_to_poolC_reference("ES281_Pentagramma_semipallida", "A", "SQD1", bt, gmap; cluster_num=33),
    switch_tip_to_poolC_reference("ES284_Pentagramma_viscosa", "A", "SQD1", bt, gmap; cluster_num=243),
    switch_tip_to_poolC_reference("ES284_Pentagramma_viscosa", "A", "SQD1", bt, gmap; cluster_num=249),
    switch_tip_to_poolC_reference("ES284_Pentagramma_viscosa", "A", "SQD1", bt, gmap; cluster_num=241),
    switch_tip_to_poolC_reference("ES280_Pentagramma_glanduloviscida", "A", "SQD1", bt, gmap; cluster_num=10),
    switch_tip_to_poolC_reference("ES277_Pentagramma_rebmanii", "A", "SQD1", bt, gmap; cluster_num=156), switch_tip_to_poolC_reference("ES307_Pentagramma_triangularis", "A", "CRY2", bt, gmap; cluster_num=186),
    switch_tip_to_poolC_reference("ES307_Pentagramma_triangularis", "A", "CRY2", bt, gmap; cluster_num=142),
    switch_tip_to_poolC_reference("ES307_Pentagramma_triangularis", "A", "CRY2", bt, gmap; cluster_num=410),
    switch_tip_to_poolC_reference("ES307_Pentagramma_triangularis", "A", "CRY2", bt, gmap; cluster_num=348),
    switch_tip_to_poolC_reference("ES300_Pentagramma_rebmanii", "A", "CRY2", bt, gmap; cluster_num=200),
    switch_tip_to_poolC_reference("ES300_Pentagramma_rebmanii", "A", "CRY2", bt, gmap; cluster_num=198),
    switch_tip_to_poolC_reference("ES300_Pentagramma_rebmanii", "A", "CRY2", bt, gmap; cluster_num=205),
    switch_tip_to_poolC_reference("ES296_Pentagramma_triangularis", "A", "CRY2", bt, gmap; cluster_num=6),
    switch_tip_to_poolC_reference("ES298_Pentagramma_triangularis", "A", "CRY2", bt, gmap; cluster_num=9),
]

all_refs = vcat(
    [r[5] for r in values(results_A)]...,
    [r[5] for r in values(results_B)]...,
    manual_refs,
)

seen_headers = Set{String}()
unique_refs = String[]
for entry in all_refs
    hdr = first(split(entry, "\n"))   # ">locus=.../group=.../ref_taxon=..."
    if hdr ∉ seen_headers
        push!(seen_headers, hdr)
        push!(unique_refs, entry)
    end
end

# Append to reference.fasta
REF_FASTA = joinpath("input", "namedrefsrun", "reference.fasta")
if !isempty(unique_refs)
    open(REF_FASTA, "a") do io
        println(io) # Start with a newline to separate from existing file content
        for entry in unique_refs
            println(io, entry)
        end
    end
end

println("\nAppended $(length(unique_refs)) unique reference entr$(length(unique_refs) == 1 ? "y" : "ies") to $REF_FASTA")