#!/bin/zsh

# Pull PURC OTU alignments from Hyak for homologizer (run from project root).
#   ./scripts/homologizer/sync_hyak.sh pull   # OTU alignments  Hyak -> output/homologizer/raw/
# One rsync, so only one Duo prompt.

RUN_NAME="june3"   # ← run folder under output/ on Hyak
HYAK="laufran@klone.hyak.uw.edu"
REMOTE="/mmfs1/gscratch/tribblelab/hemionitidoid_longreads"
LOCAL="output/homologizer"

case "$1" in
  pull)
    mkdir -p "${LOCAL}/raw"
    # grabs <pool>/<LOCUS>/<LOCUS>_OTUs.aligned.fa and the _recluster/ version if present
    rsync -av --prune-empty-dirs \
      --include='purc_pool*/' --include='purc_pool*/*/' --include='purc_pool*/*/*_recluster/' \
      --include='*_OTUs.aligned.fa' --exclude='*' \
      "${HYAK}:${REMOTE}/output/${RUN_NAME}/" "${LOCAL}/raw/"
    ;;
  *)
    echo "usage: $0 pull"; exit 1
    ;;
esac
