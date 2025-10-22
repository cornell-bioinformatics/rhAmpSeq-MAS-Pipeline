# rhAmpSeq-MAS-Pipeline
MAS Pipeline for Vitisgen3 project. It uses the rhAmpseq genotyping data to predict the presence/absence of desirable alleles of major QTLs.  

#### Required python packages:

- seaborn

- xlsxwriter

- multiqc

- scipy==1.10.1

  

#### Procedure:

1. Setup a config.yaml file. A template config.yaml file is provided. The the instructions below how to modify the file. 

   ```
   ## contents of the config file ##
   #
   # FAMILY_PATH: ./
   # # pathway for each family - keep the yaml in the same directory as the hap_genotype file generated from the amplicon.py script.
   #
   # MARKER_TRAIT_PATH: all_trait_key   ## this is the directory with the prediction rules for each QTL
   #
   # MARKER_TRAIT_FILE: # list of traits for MAS - naming here needs to match file names in all_trait_key
   #     - Msex
   #     - Fsex
   #     - 5OGT
   #     - Mus
   #     - Run1
   #     - Ren1
   #     - Ren2
   #     - Ren3and9
   #     - Ren4
   #     - Ren6
   #     - Ren7
   #     - Ren10
   #     - Ren11
   #     - Ren12
   #     - Rpv3-1_2023
   #     - Rpv3-2_2023
   #     - Rpv3-3_2023
   #     - Rpv10L
   #     - Rpv12
   #     - Rpv33
   #     - Seedless
   #     - Phylloxera
   #     - AcyAnthocyanin
   #
   # RM_IDV: ~
   # # can take a list of individuals to "remove" (marks them as contaminants in excel report)
   # # this file needs to be single column, one taxa per line, names matched to hapgeno
   # # I usually use this to mark the MDS outliers in the final output
   #
   # IDV_MISS: 0.5
   # sets a missingness threshold
   #
   # # names of parental individuals (if known) so they're marked in the MDS
   # # names matched to hapgeno and listed like traits above
   # MOTHER_NAMES: ~
   # FATHER_NAMES: ~
   ## end of config file ##
   ```

   

2. Run MAS Pipeline in QC Mode

   ```
   python mas_pipeline.py -c config.yaml -o test1 -s 01
   2
   # the "-s" option tells the pipeline what steps to run; the 0 means QC mode
   # I believe Cheng told me this checks the pipeline can execute and that markers
   # in trait predictions appear linked
   
   find test1/03.marker_QC/ -not -name '*.*' -exec sed -ir 's/#//g' {} \;
   # this comments out some markers
   
   multiqc test1/05.report -c test1/05.report/multiqc_config.yaml -f -n FamilyID
   # produces the multiqc report from the QC mode
   
   ## I take a look at the multqc report and see if anything looks unusual.
   ## I look at the MDS too and if there are outliers/selfs, I pull the pcs
   ## from /test1/02.family_QC/IDV_fil.recode.vcf.mds, find the IDs and add
   ## them to a file "rm_indv" (described in config) and modify config file
   ## with path
   
   ```

   

3. Run MAS Pipeline and Produce Excel Report

   ```
   python mas_pipeline.py -c onfig.yaml -o test1 -s 12
   # rerun the pipeline again, this time without "0"; it HAS to be run in QC mode f
   irst
   
   multiqc test1/05.report -c test1/05.report/multiqc_config.yaml -f -n FamilyID
   # recreate the multiqc report
   
   cp test1/05.report/Excel_report.xlsx FamilyID_Excel_report.xlsx
   # pulls the final excel report out and renames it
   
   ```

   
