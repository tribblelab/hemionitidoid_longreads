#!/bin/zsh

# Script to align all *_OTUs_wref.fa files using MAFFT
# Creates *_OTUs_wref.aligned.fa for each file

# Locus list
LOCI=(
  APPEFP
  CRY2
  GAPCPSH
  IBR3
  INFA
  MATK
  PGIC
  PSBM
  RPL2
  SQD1
)

# Pool list
POOLS=(
  purc_poolA_ccs
  purc_poolB_ccs
)

total_aligned=0
total_skipped=0

# Process each pool and locus
for pool in "${POOLS[@]}"; do
    echo "Processing $pool..."
    
    for locus in "${LOCI[@]}"; do
        INPUT_FILE="output/${pool}/${locus}/${locus}_OTUs_wref.fa"
        OUTPUT_FILE="output/${pool}/${locus}/${locus}_OTUs_wref.aligned.fa"
        
        # Check if input file exists
        if [ ! -f "$INPUT_FILE" ]; then
            echo "  ⊗ $locus: Input file not found, skipping..."
            ((total_skipped++))
            continue
        fi
        
        # Count sequences in input
        SEQ_COUNT=$(grep -c "^>" "$INPUT_FILE")
        
        echo "  → Aligning $locus ($SEQ_COUNT sequences)..."
        
        # Run MAFFT
        mafft --auto "$INPUT_FILE" > "$OUTPUT_FILE" 2> /dev/null
        
        # Check if alignment succeeded
        if [ $? -eq 0 ] && [ -s "$OUTPUT_FILE" ]; then
            ALIGNED_COUNT=$(grep -c "^>" "$OUTPUT_FILE")
            echo "  ✓ $locus: $ALIGNED_COUNT sequences aligned → $(basename $OUTPUT_FILE)"
            ((total_aligned++))
        else
            echo "  ✗ $locus: MAFFT failed!"
            rm -f "$OUTPUT_FILE"
            ((total_skipped++))
        fi
    done
    
    echo ""
done

echo "=================================================="
echo "Summary"
echo "=================================================="
echo "Files aligned: $total_aligned"
echo "Files skipped: $total_skipped"
echo ""
echo "Output files: *_OTUs_wref.aligned.fa"
echo "=================================================="
