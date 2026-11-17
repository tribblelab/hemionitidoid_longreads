"""
rename_tree_tips.jl

Renames barcode tip labels in .tre files to sample purc_id names.

Directory structure expected:
    <root>/
        poolA/
            <gene_folder>/
                *.tre
        poolB/
            <gene_folder>/
                *.tre

Gene folder name must contain one of the gene keys (case-insensitive) to
determine which spreadsheet column to use for the barcode->number lookup.

Usage:
    julia rename_tree_tips.jl <root_dir> <spreadsheet.xlsx>

Dependencies (install once):
    using Pkg; Pkg.add(["XLSX", "DataFrames"])

Output:
    Each .tre file is renamed in-place (original backed up as *.tre.bak).
    A summary is printed to stdout.
"""

using XLSX
using DataFrames

# ── Spreadsheet column name → prefix used inside that column ────────────────
const GENE_COL_PREFIXES = Dict(
    "matK"     => "Ma",
    "pgic"     => "Pg",
    "APPEFP-C" => "Ap",
    "CRY-2"    => "Cr",
    "gapCpSh"  => "Ga",
    "SQD1"     => "Sq",
    "rpl2"     => "Rp",
    "psbM"     => "Ps",
    "infA"     => "In",
    "IBR3"     => "Ib",
)

# ── Lowercase folder-name fragment → spreadsheet column ─────────────────────
const FOLDER_TO_GENE_COL = Dict(
    "matk"    => "matK",
    "pgic"    => "pgic",
    "appefp"  => "APPEFP-C",
    "cry"     => "CRY-2",
    "gapcpsh" => "gapCpSh",
    "sqd1"    => "SQD1",
    "rpl2"    => "rpl2",
    "psbm"    => "psbM",
    "infa"    => "infA",
    "ibr3"    => "IBR3",
)

# ── Detect gene column from folder name ─────────────────────────────────────
function detect_gene_col(folder_name::String)::Union{String, Nothing}
    lower = lowercase(folder_name)
    for (fragment, col) in FOLDER_TO_GENE_COL
        if occursin(fragment, lower)
            return col
        end
    end
    return nothing
end

# ── Load spreadsheet, skipping the first header row ─────────────────────────
function load_spreadsheet(path::String)::DataFrame
    xf = XLSX.readxlsx(path)
    sheet = xf[XLSX.sheetnames(xf)[1]]      # first sheet
    data  = XLSX.eachtablerow(sheet) |> collect

    # XLSX.eachtablerow picks up the first non-empty row as the header.
    # Our sheet has a dummy top row ("pool", "matK", …) that IS the real header
    # when read with header=1 in pandas — XLSX.jl handles this the same way by
    # default, so we just read it straight into a DataFrame.
    df = DataFrame(XLSX.readtable(path, XLSX.sheetnames(xf)[1]; first_row=2)...)
    return df
end

# ── Build barcode_number → new_label for one pool ───────────────────────────
"""
    build_barcode_map(df, gene_col, pool) -> Dict{Int, String}

Rules
─────
• Barcodes 2-42: pool C row may also carry a sample.
  - Pool C non-missing in this gene column  →  "<poolAB_purc_id>-<poolC_purc_id>"
  - Pool C missing in this gene column      →  "<poolAB_purc_id>"
• Barcodes 43+: no pool C exists            →  "<poolAB_purc_id>"
"""
function build_barcode_map(df::DataFrame, gene_col::String, pool::String)::Dict{Int, String}
    prefix = GENE_COL_PREFIXES[gene_col]
    pattern = Regex("^$(prefix)_(\\d+)")

    # per_barcode[bnum][pool_letter] = purc_id
    per_barcode = Dict{Int, Dict{String, String}}()

    for row in eachrow(df)
        val = row[Symbol(gene_col)]
        ismissing(val) && continue
        val_str = string(val)
        m = match(pattern, val_str)
        isnothing(m) && continue

        bnum  = parse(Int, m.captures[1])
        p     = string(row[:pool])
        purc  = string(row[:purc_id])

        if !haskey(per_barcode, bnum)
            per_barcode[bnum] = Dict{String, String}()
        end
        per_barcode[bnum][p] = purc
    end

    mapping = Dict{Int, String}()
    for (bnum, pools) in per_barcode
        pool_label = get(pools, pool, nothing)
        isnothing(pool_label) && continue          # barcode not in this pool

        c_label = get(pools, "C", nothing)
        mapping[bnum] = isnothing(c_label) ? pool_label : "$(pool_label)-$(c_label)"
    end
    return mapping
end

# ── Rename tips in a tree string ─────────────────────────────────────────────
"""
    rename_tips(tree_text, barcode_map) -> (new_text, n_replaced, warnings)

Replaces every `barcodeNNseq_...` tip label with the mapped name.
Labels with no mapping are left unchanged and added to warnings.
"""
function rename_tips(tree_text::String, barcode_map::Dict{Int, String})
    replaced  = 0
    warnings  = String[]

    # Pattern: barcode<digits>seq_<anything up to a Newick delimiter>
    pattern = r"barcode(\d+)(seq_[^(),;:]+)"

    new_text = replace(tree_text, pattern => function(m)
        bnum   = parse(Int, m.captures[1])
        suffix = m.captures[2]
        if haskey(barcode_map, bnum)
            replaced += 1
            return barcode_map[bnum] * suffix
        else
            push!(warnings, "  barcode$(lpad(bnum, 2, '0')) — no mapping found, left unchanged")
            return m.match
        end
    end)

    return new_text, replaced, warnings
end

# ── Walk the directory tree and process every .tre file ──────────────────────
function process_directory(root::String, spreadsheet::String)
    df = load_spreadsheet(spreadsheet)

    pool_folders = [("poolA", "A"), ("poolB", "B")]
    total_files    = 0
    total_replaced = 0

    for (folder_name, pool_letter) in pool_folders
        pool_dir = joinpath(root, folder_name)
        if !isdir(pool_dir)
            println("[SKIP] $pool_dir not found")
            continue
        end

        for entry in sort(readdir(pool_dir; join=true))
            isdir(entry) || continue
            gene_dir  = entry
            folder    = basename(gene_dir)

            gene_col = detect_gene_col(folder)
            if isnothing(gene_col)
                println("[SKIP] Cannot determine gene for folder: $gene_dir")
                continue
            end

            tre_files = filter(f -> endswith(f, ".tre"), readdir(gene_dir; join=true))
            if isempty(tre_files)
                println("[SKIP] No .tre files in $gene_dir")
                continue
            end

            barcode_map = build_barcode_map(df, gene_col, pool_letter)

            for tre_path in tre_files
                tree_text = read(tre_path, String)

                # Back up original
                bak_path = tre_path * ".bak"
                cp(tre_path, bak_path; force=true)

                new_text, n, warns = rename_tips(tree_text, barcode_map)
                write(tre_path, new_text)

                total_files    += 1
                total_replaced += n
                rel = relpath(tre_path, root)
                println("[OK] $rel: $n tip(s) renamed")
                for w in warns
                    println(w)
                end
            end
        end
    end

    println("\nDone. $total_files file(s) processed, $total_replaced tip(s) renamed total.")
    println("Originals backed up as *.tre.bak")
end

# ── Entry point ──────────────────────────────────────────────────────────────
function main()
    if length(ARGS) != 2
        println("""
Usage: julia rename_tree_tips.jl <root_dir> <spreadsheet.xlsx>

Dependencies (install once):
    using Pkg; Pkg.add(["XLSX", "DataFrames"])
        """)
        exit(1)
    end

    root_dir  = ARGS[1]
    xlsx_path = ARGS[2]

    isdir(root_dir)   || (println("Error: $root_dir is not a directory"); exit(1))
    isfile(xlsx_path) || (println("Error: $xlsx_path not found");          exit(1))

    process_directory(root_dir, xlsx_path)
end

main()
