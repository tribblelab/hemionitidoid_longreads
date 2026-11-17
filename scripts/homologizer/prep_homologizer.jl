#!/usr/bin/env julia
"""
prep_homologizer.jl

Turns PURC OTU alignments (pulled from Hyak with sync_hyak.sh) into
homologizer-ready aligned FASTA files + the Rev snippets that set the initial phase
and the phasing moves.

Steps per locus:
  1. read the pool A and pool B OTU alignments, strip gaps
  2. for each included sample keep the top-k OTUs by read size (size >= MIN_SIZE)
       nuclear: k = n_copies from taxa.csv; plastid: k = 1
  3. rename: diploids -> <purc_id>, polyploids -> <code>_copy1..k
  4. realign pools together with MAFFT, write aligned FASTA

Usage (from project root):
    julia --project=. scripts/homologizer/prep_homologizer.jl

First run writes input/homologizer/taxa.csv from the sampling sheet and exits;
edit it (include / n_copies) and rerun.

Output (output/homologizer/):
    aln/<LOCUS>.fa         alignments read by homologizer.Rev
    rev/loci.Rev           alignment paths + nuclear/plastid flags
    rev/init_phase.Rev     setHomeologPhase / addMissingTaxa calls
    rev/phase_moves.Rev    mvHomeologPhase moves (nuclear loci only)
    copy_summary.csv       which OTU each copy name refers to, per locus
"""

using CSV
using DataFrames
using Printf

# ── Configuration ───────────────────────────────────────────────────────────
const RAW_DIR      = joinpath("output", "homologizer", "raw")
const OUT_DIR      = joinpath("output", "homologizer")
const TAXA_CSV     = joinpath("input", "homologizer", "taxa.csv")
const SAMPLING_CSV = joinpath("input", "pentagramma_sampling_lauren_2026-04-14.csv")
const POOLS        = ["purc_poolA", "purc_poolB"]
const MIN_SIZE     = 10      # drop OTUs supported by fewer reads than this

const NUCLEAR_LOCI = ["APPEFP", "CRY2", "GAPCPSH", "IBR3", "PGIC", "SQD1"]
const PLASTID_LOCI = ["MATK", "PSBM", "RPL2", "INFA"]

const SUBGENOME_LABELS = 'A':'Z'

# ── Taxon-name helpers ──────────────────────────────────────────────────────
# keep '-' out of every label (tip names are reused in the Rev snippets)
clean(s::AbstractString) = replace(s, "-" => "_")
sample_code(purc_id::AbstractString) = clean(first(split(purc_id, "_")))

# OTU headers look like  <purc_id>_Cluster85;size=1237  (or _size_1237 after IQ-TREE)
function parse_otu_name(name::AbstractString)
    m = match(r"^(.+?)_Cluster\d+", name)
    sample = isnothing(m) ? first(split(name, ';')) : m.captures[1]
    s = match(r"size[=_](\d+)", name)
    size = isnothing(s) ? 0 : parse(Int, s.captures[1])
    return String(sample), size
end

# ── FASTA I/O ──────────────────────────────────────────────────────────────
function read_fasta(path::String)
    recs = Pair{String,String}[]
    name, seq = "", IOBuffer()
    for line in eachline(path)
        line = strip(line)
        isempty(line) && continue
        if startswith(line, ">")
            isempty(name) || push!(recs, name => String(take!(seq)))
            name = line[2:end]
        else
            print(seq, line)
        end
    end
    isempty(name) || push!(recs, name => String(take!(seq)))
    return recs
end

function write_fasta(path::String, recs)
    open(path, "w") do io
        for (n, s) in recs
            println(io, ">", n)
            println(io, s)
        end
    end
end

# ── taxa.csv template ───────────────────────────────────────────────────────
default_copies(ploidy) = ploidy == "4" ? 2 : ploidy == "3" ? 3 : 1

function write_taxa_template()
    s = CSV.read(SAMPLING_CSV, DataFrame; types=String)
    rename!(s, strip.(lstrip.(names(s), '﻿')))
    s.ploidy = coalesce.(s.ploidy, "")
    t = DataFrame(
        purc_id  = s.purc_id,
        pool     = s.pool,
        ploidy   = s.ploidy,
        n_copies = default_copies.(s.ploidy),
        include  = [p in ("2", "3", "4") && pool != "C" ? "yes" : "no"
                    for (p, pool) in zip(s.ploidy, s.pool)],
    )
    sort!(t, [:include, :purc_id]; rev=[true, false])
    mkpath(dirname(TAXA_CSV))
    CSV.write(TAXA_CSV, t)
    println("""
    Wrote template $TAXA_CSV ($(nrow(t)) samples, $(count(==("yes"), t.include)) included).
    Edit `include` (yes/no) and `n_copies` (subgenomes per sample; 1 = diploid),
    then rerun this script.""")
end

# ── Locus file lookup (prefer reclustered alignment) ────────────────────────
function find_alignment(pool::String, locus::String)
    for p in (joinpath(RAW_DIR, pool, locus, "$(locus)_recluster", "$(locus)_OTUs.aligned.fa"),
              joinpath(RAW_DIR, pool, locus, "$(locus)_OTUs.aligned.fa"))
        isfile(p) && return p
    end
    return nothing
end

# ── Main ────────────────────────────────────────────────────────────────────
function main()
    isfile(TAXA_CSV) || (write_taxa_template(); return)
    Sys.which("mafft") === nothing && error("mafft not found on PATH")

    taxa = CSV.read(TAXA_CSV, DataFrame; types=Dict(:purc_id => String, :include => String))
    taxa = taxa[lowercase.(strip.(taxa.include)) .== "yes", :]
    ncopies = Dict(r.purc_id => Int(r.n_copies) for r in eachrow(taxa))
    polyploids = sort([p for (p, k) in ncopies if k > 1])
    println("$(nrow(taxa)) samples included ($(length(polyploids)) polyploid)")

    aln_dir   = joinpath(OUT_DIR, "aln");   mkpath(aln_dir)
    rev_dir   = joinpath(OUT_DIR, "rev");   mkpath(rev_dir)
    tmp_dir   = mktempdir()

    loci     = String[]
    nuclear  = Bool[]
    present  = Dict{String,Set{String}}()       # locus => copy names with data
    seen     = Set{String}()                    # samples with data at ≥1 locus
    summary  = DataFrame(locus=String[], purc_id=String[], name=String[],
                         otu=String[], size=Int[], n_otus=Int[])

    for locus in vcat(NUCLEAR_LOCI, PLASTID_LOCI)
        is_nuc = locus in NUCLEAR_LOCI
        files = filter(!isnothing, [find_alignment(p, locus) for p in POOLS])
        if isempty(files)
            println("[SKIP] $locus: no alignments in $RAW_DIR")
            continue
        end

        # sample => [(size, otu_name, ungapped_seq), ...]
        by_sample = Dict{String,Vector{Tuple{Int,String,String}}}()
        for f in files, (name, seq) in read_fasta(f)
            sample, size = parse_otu_name(name)
            haskey(ncopies, sample) || continue
            push!(get!(by_sample, sample, []), (size, name, uppercase(replace(seq, r"[-.]" => ""))))
        end

        recs = Pair{String,String}[]
        present[locus] = Set{String}()
        for (sample, otus) in sort(collect(by_sample); by=first)
            k = is_nuc ? ncopies[sample] : 1
            keep = first(filter(o -> o[1] >= MIN_SIZE, sort(otus; rev=true)), k)
            for (i, (size, otu, seq)) in enumerate(keep)
                name = ncopies[sample] == 1 ? clean(sample) : "$(sample_code(sample))_copy$i"
                push!(recs, name => seq)
                push!(present[locus], name)
                push!(summary, (locus, sample, name, otu, size, length(otus)))
                push!(seen, sample)
            end
        end
        if length(recs) < 4
            println("[SKIP] $locus: only $(length(recs)) sequences after filtering")
            continue
        end

        raw = joinpath(tmp_dir, "$locus.fa")
        aln = joinpath(tmp_dir, "$locus.aln.fa")
        write_fasta(raw, recs)
        run(pipeline(`mafft --auto --quiet $raw`; stdout=aln))
        aligned = [n => uppercase(s) for (n, s) in read_fasta(aln)]
        write_fasta(joinpath(aln_dir, "$locus.fa"), aligned)

        push!(loci, locus); push!(nuclear, is_nuc)
        @printf("[OK] %-8s %-7s %3d seqs, %4d sites  (%s)\n", locus,
                is_nuc ? "nuclear" : "plastid", length(aligned),
                length(last(first(aligned))), join(basename.(dirname.(files)), ", "))
    end
    isempty(loci) && error("no loci written")

    missing_all = setdiff(keys(ncopies), seen)
    isempty(missing_all) || println("[WARN] no data at any locus (drop from taxa.csv): ",
                                    join(sort(collect(missing_all)), ", "))
    polyploids = filter(in(seen), polyploids)

    # ── Rev snippets ────────────────────────────────────────────────────────
    open(joinpath(rev_dir, "loci.Rev"), "w") do io
        println(io, "# generated by prep_homologizer.jl — do not edit")
        println(io, "locus_names = [", join(repr.(loci), ", "), "]")
        println(io, "alignments = [", join(repr.(joinpath.(aln_dir, loci .* ".fa")), ",\n              "), "]")
        println(io, "is_nuclear = [", join(nuclear, ", "), "]")
    end

    open(joinpath(rev_dir, "init_phase.Rev"), "w") do io
        println(io, "# generated by prep_homologizer.jl — do not edit")
        println(io, "# Polyploid copies -> mul-tree tips. Plastid copy1 is pinned to tip A (no moves),")
        println(io, "# so tip A is the maternal subgenome. Missing copies become BLANK sequences.")
        for (i, locus) in enumerate(loci), sample in polyploids
            code, tipbase = sample_code(sample), clean(sample)
            println(io, "\n# $locus  $sample")
            for c in 1:ncopies[sample]
                tip = "$(tipbase)_$(SUBGENOME_LABELS[c])"
                name = "$(code)_copy$c"
                if !(name in present[locus])
                    name = "$(code)_BLANK$c"
                    println(io, "data[$i].addMissingTaxa(\"$name\")")
                end
                println(io, "data[$i].setHomeologPhase(\"$name\", \"$tip\")")
            end
        end
    end

    open(joinpath(rev_dir, "phase_moves.Rev"), "w") do io
        println(io, "# generated by prep_homologizer.jl — do not edit")
        for (i, locus) in enumerate(loci)
            nuclear[i] || continue
            println(io, "\n# $locus")
            for sample in polyploids
                tips = ["$(clean(sample))_$(SUBGENOME_LABELS[c])" for c in 1:ncopies[sample]]
                for a in 1:length(tips)-1, b in a+1:length(tips)
                    println(io, "moves[++mvi] = mvHomeologPhase(ctmc[$i], \"$(tips[a])\", \"$(tips[b])\", weight=phase_weight)")
                end
            end
        end
    end

    CSV.write(joinpath(OUT_DIR, "copy_summary.csv"), summary)
    rm(tmp_dir; recursive=true)
    println("\nWrote $(length(loci)) alignments to $aln_dir and Rev snippets to $rev_dir")
end

main()
