#!/bin/bash
outputdir="$2"
mismatchespath="$3"
mkdir -p $outputdir

# Check if an argument is provided
if [ "$1" == "nofilter" ]; then
    for file in $(ls $mismatchespath | grep "mismatches"); do
        cp $mismatchespath/"$file" $outputdir/no_filter_"$file"
    done
elif [ "$1" == "bed" ]; then
    for file in $(ls $mismatchespath | grep "mismatches"); do
        python args_bed_filtering.py $mismatchespath/"$file" $outputdir/filter_bed_"$file" twist1.bed
    done
elif [ "$1" == "rsID" ]; then
    for file in $(ls $mismatchespath | grep "mismatches"); do
        python args_rsID_filtering.py $mismatchespath/"$file" $outputdir/filter_rsID_"$file" listrs.txt
    done
elif [ "$1" == "genes-transcripts" ]; then
    for file in $(ls $mismatchespath | grep "mismatches"); do
        python args_gene-transcript_filtering.py $mismatchespath/"$file" $outputdir/filter_genes-transcripts_"$file" genes-transcripts.txt
    done
fi
