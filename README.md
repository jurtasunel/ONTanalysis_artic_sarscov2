#!bin/bash
############################################################################################################################################################
### This script runs the artic pipeline to produce a medaka consensus from ONT reads, then calls plot_depth_ONT.R and plot_coverage_ONT.R.               ###
### Plot_depth.R takes a depth_txt.file and a quality_report.csv as inputs. Plot_converage.R takes the medaka consensus and the quality_report as input. ###
### The scripts should be on the same directory because the Rscripts use previous output files as input arguments.                                       ###
### The barcode folders containing their ONT fastqsdirectory must be on the same directory. The folders must be named with the barcode as name.          ###
############################################################################################################################################################
##########################################################################################################
#### This script requires the artic environment to be activated with "source activate artic-ncov2019" ####
##########################################################################################################

### Documentation:
# Artic instalation: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
# Primer shcemes: https://github.com/phac-nml/primer-schemes
# Primer schemes documentation: https://psy-fer.github.io/interARTIC/primers/
# Artic version of muscle: https://ubuntu.pkgs.org/20.04/ubuntu-universe-amd64/muscle_3.8.1551-2build1_amd64.deb.html
# Muscle conda error: https://github.com/artic-network/artic-ncov2019/issues/93
# Artic common issues: https://github.com/artic-network/artic-ncov2019/issues

# Define variables for the barcode directories path, the prefix identifier of the analysis and the scripts directory:
barcode_path="/home/josemari/Desktop/Jose/Tests/ONt_test"
prefix="NetoVIR_run2"
Scripts="/home/josemari/Desktop/Jose/Projects/ONT_analysis/Scripts"

# Create a variable named barcode and store the names of the barcode directories. 
files=(${barcode_path}/*) # files will store the full path of the directories.
barcode=() # barcode is an empty array.
for i in "${!files[@]}"; do # loop throught the names in files.
    filename="$(basename "${files[i]}")" # extract the dir name from the path with basename.
    barcode+=("$filename") # Append it to the barcode array.
done
# Print the barcode names
echo -e "Barcodes: ${barcode[@]}\n"

# Create empty quality report csv file with headers that will be appended during the following iterations.
echo -e "Creating empty quality_report.csv file...\n"
echo -e "prefix_barcode,average_depth,perc_over_30X" > quality_report.csv

# Loop through the barcodes:
for i in ${barcode[@]}; do

# Quality check for one barcode with artic and concatenate all fasqs in one single file.
echo -e "Running artic guppyplex for ${i} \n"
artic guppyplex --min-length 300 --max-length 1400 --directory ${barcode_path}/${i} --prefix ${prefix}

# Run medaka with artic primers to produce consensus sequence.
echo -e "Running artic medaka for ${i} \n"
artic minion --medaka --medaka-model r941_min_high_g360 --normalise 200 --threads 8 --scheme-directory ~/artic-ncov2019/primer_schemes/midnight/ --read-file ${prefix}_${i}.fastq nCoV-2019/V1 ${prefix}_${i}

# Get depth text file from primmer trimmed sorted bam file produced by medaka with samtools.
echo -e "Getting depth from ${i} bam file...\n"
samtools depth -a -H ${prefix}_${i}.primertrimmed.rg.sorted.bam -o ${prefix}_${i}_depth.txt

# Call the Rscript to generate the depth plot. It will also append the quality report with the current depth and produce a tsv file with the low depth positions.
echo -e "Calling plot_depth.R...\n" 
Rscript plot_depth_ONT.R `pwd`/${prefix}_${i}_depth.txt `pwd`/quality_report.csv

# Make directory to store results.
mkdir ${i}_result
mv ${prefix}_${i}_depth.pdf ${prefix}_${i}_low_depth_positions.csv ${prefix}_${i}.sorted.bam ${prefix}_${i}.merged.vcf ${i}_result
# Copy the consensus for each barcode instead of moving it from directory because is needed for final coverage plot.
cp ${prefix}_${i}.consensus.fasta ${i}_result

done

# Get all consensus into one fasta.
echo -e "Merging consensus into one fasta\n"
cat *.consensus.fasta > ${prefix}_consensus_genomes.fasta

# Rscript to plot coverage. It will also append the quality report with the coverage of all samples.
echo -e "Calling plot_coverage.R...\n" 
Rscript plot_coverage_ONT.R `pwd`/${prefix}_consensus_genomes.fasta `pwd`/quality_report.csv

# Rename qualiy report and coverage barplot.
mv quality_report.csv ${prefix}_quality_report.csv
mv Coverage_barplot ${prefix}_coverage_barplot.png

# Make directory to store all barcodes and results and remove intermediate files.
mkdir ${prefix}_results
mv ${prefix}_quality_report.csv ${prefix}_consensus_genomes.fasta ${prefix}_coverage_barplot.png ${prefix}_results
rm ${prefix}*









