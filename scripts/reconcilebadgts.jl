using PhyloNetworks

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
function find_outgroup_tip(t::HybridNetwork, priority_genera::Vector{String})::Union{String, Nothing}
    tip_names = [n.name for n in t.node if n.leaf]
    for genus in priority_genera
        match = findfirst(label -> get_genus(label) == genus, tip_names)
        match !== nothing && return tip_names[match]
    end
    return nothing
end

# Find the MRCA of all Pentagramma tips by finding the deepest internal node
# whose subtree contains ALL Pentagramma tips and no more than necessary.
function find_mrca(penta_tips::Vector{String}, t::HybridNetwork)::Union{PhyloNetworks.Node, Nothing}
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
function find_pentagramma_base_tip(t::HybridNetwork)::Union{String, Nothing}
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
        ([:Ynesmexia],                         :ynesmexia),
        ([:Myriopteris],                        :myriopteris),
        ([:Pellaea, :Paragymnopteris],          :pellaea_paragymnopteris),
        ([:Notholaena, :Cheiloplecton],         :notholaena_cheiloplecton),
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
 
function find_mixed_clades!(n::PhyloNetworks.Node, parent::Union{PhyloNetworks.Node, Nothing},
                              results::Vector{Vector{String}})
    n.leaf && return
 
    below = tips_below(n)
    has_penta    = any(is_pentagramma, below)
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
    process_gene_tree(path, pool, bt=nothing; max_print=60)
 
Read a newick gene tree, root it, report any mixed Pentagramma clades,
and for each mixed clade check whether each Pentagramma tip shares its
barcode with any non-Pentagramma tip in the same clade.
 
- `path`      : path to the .tre file
- `pool`      : `"A"` or `"B"` (which pool this gene tree came from)
- `bt`        : optional `BarcodeTable`; if provided, barcode reconciliation is run
- `max_print` : suppress per-tip output for clades with more than this many
                Pentagramma individuals (default 60); the clade header is always printed
 
Returns `(tree, rooting_strategy, mixed_clades, barcode_matches)` where
`barcode_matches` is a `Vector{Vector{NamedTuple}}`, one entry per mixed clade,
each entry being a list of `(penta_tip, nonpenta_tip, poolC_id, match::Bool)`.
"""
function process_gene_tree(path::String, pool::String,
                           bt::Union{BarcodeTable,Nothing}=nothing;
                           max_print::Int=60)
    t = readTopology(path)
    strategy = root_gene_tree!(t)
    gene = splitext(basename(path))[1]   # e.g. "APPEFP" from "APPEFP.tre"
 
    clades = mixed_pentagramma_clades(t)
    barcode_matches = Vector{Vector{NamedTuple}}()
 
    if isempty(clades)
        println("No mixed Pentagramma/outgroup clades found.")
    else
        println("$(length(clades)) mixed clade(s) found:")
        for (i, clade) in enumerate(clades)
            penta_tips    = filter(is_pentagramma, clade)
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
                    # A barcode match means the pool C id for this Pentagramma
                    # starts with the same collection code as the non-Pentagramma tip
                    is_match = poolC_id !== nothing &&
                               startswith(poolC_id, ntip_code)
 
                    push!(clade_matches, (
                        penta_tip    = ptip,
                        nonpenta_tip = ntip,
                        poolC_id     = poolC_id,
                        match        = is_match,
                    ))
 
                    if verbose
                        if bt !== nothing
                            match_str = is_match ? "✓ BARCODE MATCH" : "✗ no match"
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
 
            # Always print a match summary for large clades
            if !verbose
                n_matches = count(m -> m.match, clade_matches)
                println("    $(n_matches)/$(length(clade_matches)) barcode match(es) found")
            end
 
            push!(barcode_matches, clade_matches)
        end
    end
 
    return (t, strategy, clades, barcode_matches)
end
 
# ---------------------------------------------------------------------------
# Batch: process all gene trees in a directory
# ---------------------------------------------------------------------------
"""
    process_gene_tree_dir(dir::String, pool::String; ext=".tre", bt=nothing)
 
Recursively process all gene tree files under `dir`.

- `pool` : `"A"` or `"B"` — which sequencing pool these trees came from
- `bt`   : optional `BarcodeTable` for barcode reconciliation
- `ext`  : file extension to match (default `.tre`)

Returns a `Dict` mapping filename => `(tree, rooting_strategy, mixed_clades, barcode_matches)`.
"""
function process_gene_tree_dir(dir::String, pool::String;
                               ext::String=".tre",
                               bt::Union{BarcodeTable,Nothing}=nothing,
                               results::Dict{String,Any}=Dict{String,Any}())
    for entry in readdir(dir, join=true)
        if isdir(entry)
            process_gene_tree_dir(entry, pool; ext=ext, bt=bt, results=results)
        elseif endswith(entry, ext)
            println("\n=== $(basename(entry)) ===")
            results[basename(entry)] = process_gene_tree(entry, pool, bt)
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
    "appefp"  => "APPEFP",
    "ibr3"    => "IBR3",
    "sqd1"    => "SQD1",
    "cry2"    => "CRY2",
    "matk"    => "matK",
    "pgic"    => "pgic",
    "gapcpsh" => "gapCpSh",
    "rpl2"    => "rpl2",
    "psbm"    => "psbM",
    "infa"    => "infA",
)
 
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
                          bt::BarcodeTable)::Union{String, Nothing}
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
        @warn "No pool C individual found with barcode '$barcode' for gene $gene"
        return nothing
    end
 
    return bt.rows[c_row]["purc_id"]
end

rename_treefile_to_tre("output/genetrees/namedrun")

bt = load_barcode_table("input/barcode_to_samples.csv")

results = process_gene_tree_dir("output/genetrees/namedrun/poolA/", "A"; bt=bt)
results = process_gene_tree_dir("output/genetrees/namedrun/poolB/", "B"; bt=bt)