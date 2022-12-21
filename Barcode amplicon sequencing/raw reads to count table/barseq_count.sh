#this code will extract the read count of cell barcodes from fastq files.
#it uses the bartender pipeline. For instructions on how to install it see: https://github.com/LaoZZZZZ/bartender-1.1
#after monterey update the path /usr/bin/python was removed in MacOS and thus, the bartender scripts stopped workinig. 
#To fix it, edit the scripts first by replacing the first line with '/usr/bin/env python' instead of '/usr/bin/python'


#variables
WORK_DIR="/Volumes/LaCie/sequencing_data/MPU_data/lineage_tracing/per_experiment/MIL_FS_BPK282_barcoded_dorien/fastp_outputs/processed_reads"
INPUTS_PATH="/Volumes/LaCie/sequencing_data/MPU_data/lineage_tracing/per_experiment/MIL_FS_BPK282_barcoded_dorien/fastp_outputs/processed_reads"
SPLIT_RIGHT="_merged"
READ_ID="_merged"

#parameters
#BC_PATTERN is used to find the barcodes. It has a specific structure: nucleotides preceeding the barcode, the number of random nucleotides between brackets,
#the fixed nucleotides in between, and the succeeding nucleotides. Thus, for a barcode like "GCCTNNNNNAANNNNNAANNNNNTTNNNNNAGGT" the pattern should be: "GCCT[5]AA[5]AA[5]TT[5]AGGT"
#It must be always preceeded and succeeded by fixed nucleotides! Maximum of 5 fixed nucleotides in both positions.
#You can allow for mismatches in the random sequences by giving a numerical interval instead: "GCCT[4-5]AA[4-5]AA[4-5]TT[4-5]AGGT" will accepts barcodes with 4 or 5 random nucleotides
#in each position. 

#mCerBC-CapSeq-V1 patterns 
BC_PATTERN="AGCCT[5]AA[5]AA[5]TT[5]AGGTC"
BC_PATTERN="AGCCT[26]AGGTC"
#dorien library pattern
#BC_PATTERN="GACAG[29-31]AGCAG"
#more relaxed pattern
BC_PATTERN="GACAG[28-39]AGCAG"

#number of mismatches allowed in the preceeding and succeeding fixed nucleotides. Default is 2. 
MISMATCHES=2
#minimum number of different nucleotides between two barcodes to not cluster them. 
MIN_DISTANCE=2
#Barcodes with read count lower than MIN_COUNT will be removed.
MIN_COUNT=2
#DIRECTION tells the code to read the fastq reads in the 'forward', 'reverse_complement' or 'both' directions. 
DIRECTION="forward"
#NUMBER of CPU cores
THREADS=4

#moving to working directory
cd "$WORK_DIR"

#making output directories
mkdir bartender_outputs
mkdir bartender_outputs/read_count_tables
mkdir bartender_outputs/other_files


#getting samples
SAMPLES=($(find "$INPUTS_PATH" -maxdepth 1 -iname '*.fastq*' -exec basename {} \;| \
sed "s/$SPLIT_RIGHT/+/" | \
cut -d "+" -f 1 | \
sort | \
uniq))

echo "found '${#SAMPLES[@]}' samples"

#finding the fastq file of each sample to extract the barcodes
for SAMPLE in ${SAMPLES[@]}
do
FILE=$(find "$INPUTS_PATH" -maxdepth 1 -name "*$SAMPLE*$READ_ID*fastq*")
echo "Extracting barcodes from sample" $SAMPLE

bartender_extractor_com -f "$FILE" -o "${SAMPLE}_extracted" \
-p $BC_PATTERN \
-m $MISMATCHES

bartender_single_com -f "$SAMPLE"_extracted_barcode.txt -o "${SAMPLE}" \
-d $MIN_DISTANCE \
-c $MIN_COUNT \
-z -1 \
-t $THREADS 

mv "${SAMPLE}"*_cluster.csv "bartender_outputs/read_count_tables"
mv "${SAMPLE}"*.csv "bartender_outputs/other_files"
mv "${SAMPLE}"*.txt "bartender_outputs/other_files"
done


