## ---------------------------
## Purpose of script: Script to run MMseqs2 on input fasta file against balanced assembly database
## Author: Henri Chung
## Date Created: 2021-01-12
## Date Modified: 2024-01-04
## ---------------------------
INPUT=$1
START=$(date +%s.%N)
FASTA_FOLDER="/data/protein_fa"
OUTPUT_FOLDER="/data/mmseqs_results/$INPUT"

# Create a folder with the same name as the input, if it does not already exist
if [[ ! -d "$OUTPUT_FOLDER" ]]; then
    mkdir "$OUTPUT_FOLDER"
    echo "Created folder: $OUTPUT_FOLDER"
    echo "Query Accession: $INPUT "
else
    echo "Folder already exists: $OUTPUT_FOLDER"
fi
# Find the first file in the specified folder that partially matches the input
FILE=$(find "$FASTA_FOLDER" -type f -name "*$INPUT*" | head -n 1)
OUTFILE="$OUTPUT_FOLDER/${INPUT}_S7_5.mmseq"
if [ ! -f "$OUTFILE" ]; then
    echo "Running mmseq at sensitive"
    zgrep -f compare/data/balanced_fusion_dtm_protein_names.txt $FILE > $FILE"temp"
    sed -i 's/>//g' $FILE"temp"
    seqkit grep -n -j 16 --delete-matched -f $FILE"temp" $FILE -o $FILE"temp".match
    mmseqs easy-search -s 7.5 -v 0 --k-score 85 -e 10000.0 --max-seqs 4000 --format-output "query,target,bits" $FILE"temp".match compare/mmseqsDB/balanced_assemblyDB "$OUTPUT_FOLDER/${INPUT}_S7_5.mmseq" compare/mmseqsDB/balanced_assemblyDB_tmp
    rm $FILE"temp".match
    rm $FILE"temp"
else
    echo "Match found: $OUTFILE"
fi

END=$(date +%s.%N)
RUNTIME=$(echo "$END - $START" | bc)
FORMATTED_RUNTIME=$(printf "%.1f\n" $RUNTIME)
echo "Runtime: $FORMATTED_RUNTIME seconds"
exit

