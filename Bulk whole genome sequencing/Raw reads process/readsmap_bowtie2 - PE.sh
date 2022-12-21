#specifying requirements
#########################################################
#PBS -l walltime=23:59:00
#PBS -L tasks=1:lprocs=28

threads=4

#Description
#This code will take pair-end fastq files, map them to a reference genome, return an indexed and sorted bam file as well as a tab delimited table with binned read count and gc content for downstream analysis.  

#first specify the variables
WORK_DIR="/Volumes/LaCie/sequencing_data/MPU_data/bulk_wga/temp" #specifies in which directory this code will work. 
GENOME_PATH="/Users/gnegreira/OneDrive - ITG/ITM/PhD/IGV Files/Leishmania donovani/BPK282/Genomes/LdPBQ7G3I2I8.fasta" #points to the reference genome fasta file
INPUTS_PATH="/Volumes/LaCie/sequencing_data/MPU_data/bulk_wga/temp" #the directory where the fastq files are
INDEX_NAME="BPK282_index" #when the code builds the indexes for the reference genome, it will use this variable to name the indexes files. Don't use spaces.
SAMPLE_NAME_SPLIT="_R" #this is used to get samples name according to the names of the fastq files. This variable is used to split the name of the files and use the left size to determine samples name.  
READ1_ID="R1.fastq.gz" #tells the code how to look for the fastq file of read 1
READ2_ID="R3.fastq.gz" #tells the code how to look for the fastq file of read 2
BIN_SIZE=2500 #defines the bin size for the read count file. 

#going to the working directory
cd "$WORK_DIR"

#Loading modules
module load bowtie2
module load samtools
module load deeptools
module load bedtools

#making output directory
mkdir outputs

#finding fastqfiles and determining number an name of samples
#The name of the samples are based on the fastq files names. It asumes read1 and read2 from the same sample will start with the same name.

#first it will get the name of the folder where inputs are. This is useful for substringing the samples name later. 
INPUTS_FOLDER_NAME=$(echo "$INPUTS_PATH" | awk -F / '{print $NF}')

#then it uses the find function to find all fastq files. The files path are removed based on the name of the inputs folder,
#and the name is substringed based on the left side of the string defined in $SAMPLE_NAME_SPLIT, usually "_L00".
#So in summary, anything before the "_L00" string in the file name will define the samples name.
#So this assumes that for each sample, there will be two files, one for read 1 and another for read 2, and that both will start with the same name
#until the "L00" string.

SAMPLE_NAMES=$(find "$INPUTS_PATH" -iname '*.fastq*' | \
sed "s/\/$INPUTS_FOLDER_NAME\//+/" | \
sed "s/$SAMPLE_NAME_SPLIT/+/" | \
cut -d "+" -f 2 | \
uniq -d)

#creating index for reference genome
echo "Creating index for reference genome..."
bowtie2-build "$GENOME_PATH" $INDEX_NAME

#Once samples are defined, it will now find, for each sample, the read1 and read2 files
for SAMPLE in $SAMPLE_NAMES
    do
    READ1=$(find "$INPUTS_PATH" -iname $SAMPLE*$READ1_ID)
    READ2=$(find "$INPUTS_PATH" -iname $SAMPLE*$READ2_ID)

    echo "Mapping reads for sample $SAMPLE..."

    #####Actual mapping####################
    #creating a folder for each sample
    OUTPUT_FOLDER="./outputs/$SAMPLE"
    rm -r $OUTPUT_FOLDER
    mkdir $OUTPUT_FOLDER

    #mapping the reads
    bowtie2 -p $threads \
        -x $INDEX_NAME \
        -1 "$READ1" \
        -2 "$READ2" \
        -S $OUTPUT_FOLDER/$SAMPLE.sam

    echo "Converting sam file to bam..."
    #creating a bam file from the sam file
    samtools view -S -b --threads $((threads-1)) $OUTPUT_FOLDER/$SAMPLE.sam > $OUTPUT_FOLDER/$SAMPLE.bam
    #sorting the bam file
    samtools sort $OUTPUT_FOLDER/$SAMPLE.bam -o $OUTPUT_FOLDER/$SAMPLE.sorted.bam
    #indexing the bam file
    samtools index $OUTPUT_FOLDER/$SAMPLE.sorted.bam

    echo "Creating the table with the binned read count using bins of size $BIN_SIZE bp..."
    #creating a csv file with the read cout per bin using bamCoverage

    bamCoverage -b $OUTPUT_FOLDER/$SAMPLE.sorted.bam \
        -o $OUTPUT_FOLDER/$SAMPLE.reads_per_bin.bed \
        -of bedgraph \
        -bs $BIN_SIZE \
    -p $threads 

    #using bedtools to append gc content information
    echo "Appending gc content..."
    bedtools nuc -fi "$GENOME_PATH" \
        -bed $OUTPUT_FOLDER/$SAMPLE.reads_per_bin.bed \
        | cut -f 1-4,6 > $OUTPUT_FOLDER/$SAMPLE.binned_depth_gc.csv

    echo "Done"
done

#To do: add a step to check that there is only 1 READ1 and 1 READ2 files.
#       concatenate reads from multiple lanes. 
#       Skip indexing the genome if already done. 
#       Add proper column names to the csv file.
#       Add a table with summary metrics for alignment 