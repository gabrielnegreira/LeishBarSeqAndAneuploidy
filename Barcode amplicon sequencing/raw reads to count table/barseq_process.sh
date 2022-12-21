#this code will take PE fastq files from amplicon sequencing, clean the reads and merge both reads in a single file by overlapping them. 
#it uses the fastp software. For more information visit https://github.com/OpenGene/fastp 
#OBS: put all files that must be analyzed in a clean directory. Set this directory in the INPUTS_PATH variable. 

#variables
WORK_DIR="/Volumes/LaCie/sequencing_data/MPU_data/lineage_tracing/per_experiment/PAT_FS_BPK282_barcoded_dorien_robin"
INPUTS_PATH="/Volumes/LaCie/sequencing_data/MPU_data/lineage_tracing/per_experiment/PAT_FS_BPK282_barcoded_dorien_robin"
#SPLIT_RIGHT is used to get samples name according to the names of the fastq files. This variable is used to split the name of the files and use the left size to determine samples name.
SPLIT_RIGHT="_R"  
#tells the code how to look for the fastq file of read 1 
READ1_ID="_R1" 
#tells the code how to look for the fastq file of read 2
READ2_ID="_R2" 
THREDS=4

#if you want the output files to be compressed set OUT_FORMAT to "fastq.gz". For uncompressed use "fastq"
#notice that the bartender pipeline used in the barseq_count script does not accept compressed files.
OUT_FORMAT="fastq"

#parameters: TO DO: add more options. 
AVG_MIN_Q=10
MAX_LEN=75
UMI_LOCATION=per_read
UMI_LENGTH=8
SUBSET=1000000
#OBS: UMI analysis not implemented yet

#moving to working directory
cd "$WORK_DIR"

#making output directories
mkdir fastp_outputs
mkdir fastp_outputs/processed_reads
mkdir fastp_outputs/unpaired_reads
mkdir fastp_outputs/fastp_reports

#getting samples
SAMPLES=($(find "$INPUTS_PATH" -maxdepth 1 -iname '*.fastq*' -exec basename {} \;| \
sed "s/$SPLIT_RIGHT/+/" | \
cut -d "+" -f 1 | \
sort | \
uniq))

echo "found '${#SAMPLES[@]}' samples"

#finding the R1 and R2 files of each sample
for SAMPLE in ${SAMPLES[@]}
do
FILE_R1=$(find "$INPUTS_PATH" -maxdepth 1 -name "*$SAMPLE*$READ1_ID*fastq*")
FILE_R2=$(find "$INPUTS_PATH" -maxdepth 1 -name "*$SAMPLE*$READ2_ID*fastq*")

echo "Processing files of sample" $SAMPLE

fastp --in1 $FILE_R1 --in2 $FILE_R2 \
--out1 "fastp_outputs/unpaired_reads/${SAMPLE}_unpaired_R1.fastq.gz" --out2 "fastp_outputs/unpaired_reads/${SAMPLE}_unpaired_R2.fastq.gz" \
--overrepresentation_analysis \
--average_qual $AVG_MIN_Q \
--max_len1 $MAX_LEN \
--max_len2 $MAX_LEN \
--reads_to_process $SUBSET \
--merge --merged_out "fastp_outputs/processed_reads/${SAMPLE}_merged.${OUT_FORMAT}" \
--overlap_diff_limit  10 \
--html "fastp_outputs/fastp_reports/${SAMPLE}_fastp_report.html" \
--json "fastp_outputs/fastp_reports/${SAMPLE}_fastp_report.json"
done
