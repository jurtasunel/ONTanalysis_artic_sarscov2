# ONTanalysis_artic_sarscov2

This pipeline runs the artic pipeline to maps ONT reads from tile-amplicon sarscov2 to the reference and produces consensus fasta files, depth plots and a quality report csv.
This pipeline requires the artic environment to be activated with "source activate artic-ncov2019". 

All barcode directories must be on the same data directory specified on the ONT_analysis_pipeline.sh script.

ONT_analysis_pipeline.sh is the parent script to run. It produces a result directory for each of the barcodes and a final results directory with the quality_report.csv and coverage_barplot.png for all the samples.

### Documentation:
### Artic instalation: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
### Primer shcemes: https://github.com/phac-nml/primer-schemes
### Primer schemes documentation: https://psy-fer.github.io/interARTIC/primers/
### Artic version of muscle: https://ubuntu.pkgs.org/20.04/ubuntu-universe-amd64/muscle_3.8.1551-2build1_amd64.deb.html
### Muscle conda error: https://github.com/artic-network/artic-ncov2019/issues/93
### Artic common issues: https://github.com/artic-network/artic-ncov2019/issues


