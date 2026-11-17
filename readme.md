# hemionitidoid_longreads

Analyzing PacBio longread data for Hemionitidoid ferns (Pentagramma, Doryopteris) using the PURC pipeline

## Layout

Run everything from the project root.

```text
input/            # everything scripts read: PURC configs + mapfiles, reference seqs,
                  #   sampling sheet, genus -> ref group map, homologizer taxa
scripts/
  purc.qmd        # notes/workflow for running PURC on Hyak
  homologizer.qmd # notes/workflow for phasing copies with homologizer
  purc/           # slurm jobs (CCS conversion, PURC pools A/B, recluster + IQ-TREE)
                  #   and my patched purc.py
  homologizer/    # homologizer prep / run / Hyak sync scripts
  utils/          # small julia/shell helpers (add refs, align, rename tips, ...)
  debugging/      # one-off debugging scripts
output/           # results (gitignored)
notes/            # human notes and reference spreadsheets; no script depends on these
archive/          # old/superseded material (dummy run, old genetrees, etc.; gitignored)
bin/              # local clone of upstream PURC (gitignored)
log/              # job logs (gitignored)
```
