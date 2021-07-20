#!/bin/bash

data/...
    #counts and fqn inputs for scripts
    #fqn are outputs of the normalisation script
    #TODO li_crc has both gene and ensemble.
        #avg_expr output/gene sigs are gene_ensemble
        #what do other algorithms use as input?
            #TODO check adobo and sccatch

program/runs_adobo.py #run adobo predictions and saves fqn normalised data
    (time python program/runs_adobo.py ./data/chung_breast/counts.tsv.xz ./data/chung_breast/chung_breast_clusters.csv ./program/output/cb_adobo.tsv) 2> ./times/cb_adobo.txt
    #TODO change these paths to the other datasets
    (time python program/runs_adobo.py ./data/chung_breast/counts.tsv.xz ./data/chung_breast/chung_breast_clusters.csv ./program/output/cb_adobo.tsv) 2> ./times/cb_adobo.txt
    (time python program/runs_adobo.py ./data/chung_breast/counts.tsv.xz ./data/chung_breast/chung_breast_clusters.csv ./program/output/cb_adobo.tsv) 2> ./times/cb_adobo.txt
    (time python program/runs_adobo.py ./data/chung_breast/counts.tsv.xz ./data/chung_breast/chung_breast_clusters.csv ./program/output/cb_adobo.tsv) 2> ./times/cb_adobo.txt
    (time python program/runs_adobo.py ./data/chung_breast/counts.tsv.xz ./data/chung_breast/chung_breast_clusters.csv ./program/output/cb_adobo.tsv) 2> ./times/cb_adobo.txt
    (time python program/runs_adobo.py ./data/chung_breast/counts.tsv.xz ./data/chung_breast/chung_breast_clusters.csv ./program/output/cb_adobo.tsv) 2> ./times/cb_adobo.txt
    (time python program/runs_adobo.py ./data/chung_breast/counts.tsv.xz ./data/chung_breast/chung_breast_clusters.csv ./program/output/cb_adobo.tsv) 2> ./times/cb_adobo.txt
    (time python program/runs_adobo.py ./data/chung_breast/counts.tsv.xz ./data/chung_breast/chung_breast_clusters.csv ./program/output/cb_adobo.tsv) 2> ./times/cb_adobo.txt
    #NOTE adobo sometimes throws an error when predicting cell types on macOS 
        # installations. Recommend using linux.
    #need a script that saves the filtered counts? or just provide filtered and prepped counts maybe

program/tme_inputs/avg_expr.ipynb #get avg expression from counts/fqn
        #counts.tsv.xz and fqn.tsv.xz for each dataset
        #clusters file for each dataset
        #the van galen output file gets filtered for cells that we aren't interested in
            #some cluster are also combined, make sure this is done
    #TODO some of these were run on the cluster because of memory issues
    #TODO maybe make this into a python script instead of ipynb

program/scripts/*.sh #run cluster labelling algorithms
    #TODO there are two van galen scripts, only publish one
    #Cibersort and GSEA have licenses, i can't distribute them
        #cibersort: https://cibersort.stanford.edu/download.php
        #gsea: http://software.broadinstitute.org/gsea/downloads.jsp
    #these scripts just need the compute canada stuff at the top removed
        #need r/3.5.2, perl/5.22.2, java/1.8
            #optparse, vioplot, GSA, data.table, precrec, ROCR, Seurat, dplyr, Rserve, e1071, colorRamps, stats
            #bioconductor: preprocessCore, GSVA, qvalue
            #perl_modules from javier's repo needs to be in the PERL5LIB env
                #can possibly do this in a main wrapper script?
        #could clone javier's git repo in the bash script
            https://github.com/jdime/scRNAseq_cell_cluster_labeling.git
        #need the tme_inputs folder for adobo, sccatch, and javier's stuff.
    
    program/sccatch/run_sccatch.R #sccatch
        program/sccatch/*.sh
        #TODO merge bash scripts into one
        #this'll be an r script with a bash wrapper for run time

    #TODO i don't think i can time these and get the right output file at the same time
        # had to remove the timing code from the scripts because of a weird perl error on\
        # compute canada. Might be able to fix that locally, but it'll give different results

cell_based_program/* #run cell based labelling algorithms
    #need to get ALL of ping's scripts and get them to output into the right folder
    #TODO change his scripts to read in genes x cells and transpose
    #TODO see what format his labels take in (column names, etc)
    CV.R #Generate cross validation folds for R based methods
        #TODO make sure the first column in all of my clusters files is the 
            # cell ID for cross validation
                #darmanis probably needs a change
        #arg1 is path to dataset (not actually needed)
        #arg2 is path to output folder
        #arg3 is path to the labels folder
    CV_r2py.py #convert .Rdata folds to .pkl for python
        #arg1 is the path of the folder containing the Rdata file 
            #pkl is saved in the same directory
    R_methods.R #run R methods for cell type prediction
        #arg1 is the dataset to run
        #arg2 is the relative path to ./data/...
            #TODO right now the counts files are in a subfolder
            #might need to change this
    Python_methods.py #run most python methods
        #args1 is path to data
        #args2 is the dataset
            #TODO right now the counts files are in a subfolder
            #might need to change this
            #TODO need to make sure this reads the labels properly
    run_LAmbDA.py
        #args1 is path to data
        #args2 is the dataset
    run_scVItool.py
        #args1 is path to data
        #args2 is the dataset
    other_scripts/results_table.ipynb #combine predictions into a single file
        #TODO this only gathers a single dataset right now. 
        #either make it take a command line arg or do a loop
    #TODO need to update output directory of the time_sim.py and subsequent scripts
    other_scripts/time_sim.py #combine prediction time files
        #TODO might need to update dataset names
    other_scripts/time_bar_plot.ipynb #make the supplemental bar plots
        #this reads in my runtimes as an npy file
        #TODO make sure that the cluster labelling timing scripts make an npy file
            #TODO ping sent me a script, gotta append it to my time gatherer
    other_scripts/time_coefficient_plot.ipynb #make the coef heatmap
        #TODO this needs to only output the dataframe


result_gathering.ipynb #gather predictions from cluster labelling outputs
    # iirc i was adding adobo and sccatch predictions manually... must fix
    #TODO this looks like it makes the right files but they need to be put on disk
    #just add the cluster mapping and cell-based matching to this python script
    #TODO #map cluster predictions to cell predictions
    #TODO #merge cell-based predictions with per-cell cluster predictions


#TODO #bootstrap the results
    #everything but f-measure is optional, but required for supplemental plots
    #TODO I had to change some paths in the f-measure script, make sure to do the same for everything else
    #predictions/*_predictions.tsv

#TODO timing scripts
    #will have to get scripts that read everything and generate the two plotting scripts from ping
    #my timing wrappers are separate from the wrappers that actually create the output
        #i would essentially need to run everything twice to get the time results

sensitivity_plots.ipynb #score the algorithms
    #predictions/*_predictions.tsv

main_figures.ipynb #generate MOST of the main figures
    #data_sizes.tsv #THIS GETS UPLOADED
    #Rdata/F-Measure-Bootstrap-Ensemble.tsv
    #times/df_for_heatmap.tsv #FROM PING
    #times/df_coef.tsv #FROM PING
    #performance/seurat/bigdf.tsv

subsampling_all_cells/looking_to_automate_singletons.ipynb #generate figures for imbalanced experiment
    #DATADIR IS performance/seurat
        #DATADIR/*classification_report.tsv
    #data/Lambrechts_LC_800.tsv
    #data/Peng_PC_800.tsv
    #data/vanGalan_AML_800.tsv
    #data/Darmanis_GBM_800.tsv
    #data/JA_Melanoma_800.tsv
    #data/Tirosh_Melanoma_800.tsv
    #ALL ARE FROM PING


rare_cell_types.ipynb #make the rare cell types plot
    #performance/seurat/bigdf.tsv

patient_data/predictions_results/score_patients.ipynb #generate the patient plots
    #pancreatic/Peng_patient_test.tsv
    #pancreatic/Peng_PC.tsv
    #pancreatic/Peng_PC_og_nocell.tsv
    #../pancreatic/pancreatic_patients.tsv

    #aml/vanGalan_patient_test.tsv
    #aml/vanGalan_AML.tsv
    #../../predictions/vg_predictions.tsv
    #../aml/aml_patients.tsv
    
    #metastatic_melanoma/Tirosh_patient_test.tsv
    #metastatic_melanoma/Tirosh_metastatic_melanoma.tsv
    #../../predictions/tm_predictions.tsv
    #../metastatic_melanoma/metastatic_melanoma_patients.tsv
    
    #melanoma/JA_patient_test.tsv
    #melanoma/JA_melanoma.tsv
    #../../predictions/jam_predictions.tsv
    #../melanoma/melanoma_patients.tsv

    #lung/patient_test.tsv
    #../../predictions/llc_predictions.tsv
    #../lung/lung_patients_unique.tsv 
    #../lung/lung_patient_counts_unique.tsv

supplementary_figures.ipynb #OPTIONAL, generate supplemental figures
    #performance/seurat/bigdf.tsv
    #Rdata/F-Measure-Bootstrap-Ensemble.tsv
    #subsampling_all_cells/performance/*.tsv
    #other heatmaps
        #Rdata_seurat/Homogeneity_bootstrap.tsv
        #Rdata_seurat/ARI_bootstrap.tsv
        #Rdata_seurat/percentage_correctly_assigned_bootstrap.tsv
        #Rdata_seurat/precision_bootstrap.tsv
        #Rdata_seurat/recall_bootstrap.tsv

#TODO do we need to run with both the paper clusters and the seurat clusters?
    # i think most of our analysis is with the seurat clusters but we may need a supplemental
    #TODO purge mention of the paper clusters in all scripts
    # file for the paper clusters if we mention it in the paper much
#TODO this is gonna take a long time to run. Provide an example with one small dataset that reviewers can run in order to test the full pipeline
    # chung breast will work but we don't include it in the patient analysis or subsampling
