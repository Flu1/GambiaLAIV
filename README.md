# GambiaLAIV
Code for analyzing next-gen sequencing data from LAIV from the Gambia.
This code is associated with the paper "The Consequences of Egg Adaptation in the H3N2 Component to the Immunogenicity of Live Attenuated Influenza Vaccine."

The first part of the code will demonstrate how to extract barcodes from a single sample.  The second part will show how all the samples were analyzed.

1. Sample sequences can be downloaded from https://www.ebi.ac.uk/ena (project number PRJEB34129.)
   The code here will give the example for D414H Day 2 (ERS3673683
   SAMEA5885229 is the accession number.) Download the sequences and put them in a new folder in your working directory.
2. Initial sample preparation - see initial_nextgen.R
3. Run extract_gambia_barcodes.R to extract the barcodes.

4. To make Figure 2c from the paper, you can download the cleaned sequence data which followed the steps above- see Figure2c.zip.  Then run gambia_figure() in a folder with the unzipped files.
