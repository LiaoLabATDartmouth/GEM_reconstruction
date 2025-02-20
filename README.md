# Overview
A **GEnome-scale Metabolic (GEM)** model is a computational representation of an organism’s metabolism. GEMs are widely used to predict growth, optimize metabolic pathways, and explore gene functions through constraint-based modeling approaches.

The script **run_gem_recon.py** reconstructs strain-specific GEM models from whole-genome sequencing data. It begins with a pan-genome (universe) model, which includes all metabolic reactions present in an organism. For each query genome, we perform a BLAST search of their protein sequences against those in the pan-genome model to determine the presence or absence of genes and reactions in the query GEM. Three gap-filling algorithms are implemented:
	1.	Minimizing the total number of reactions added (gapfilling_method="minrxnnum").
	2.	Maximizing the average score of added reactions (gapfilling_method="maxavescore"). This score is derived from the protein identity score in the BLAST analysis of the constituting genes.
	3.	Removing all minimal cut sets, which are the smallest sets of reactions whose removal would lead to a “no growth” phenotype. For each minimal cut set, only the reaction with the highest reaction score will be included in the reconstructed GEM.

Please install IBM CPLEX solver, which is available for free for academics. Please also set up the following parameters in the configuration file **config.yaml**:

## Required Parameters
	•	pan_gem_model_file: Genome-scale metabolic model of the universe (pan-genome) model.
	•	pan_gem_gene_id_mapping_file: A two-column table that maps gene IDs in the pan-genome model (PanGemGeneID) to those in the query genomes (QueryGeneID).
	•	pan_gem_protein_sequence_file: Protein FASTA file of the pan-genome model.
	•	query_protein_sequence_files: Protein FASTA files of the query genomes (one FASTA file per query genome).

## Optional Parameters
	•	output_folder (default: “Output”): Data output directory (each query genome has its own subfolder under this directory).
	•	protein_identity_score_cutoff (default: 95): Cutoff value used to determine whether a protein in the query genome is present in the pan-genome model.
	•	gapfilling_method (default: “mcs”): Available options include “minrxnnum,” “maxavescore,” and “mcs.”
	•	n_cpus (default: 1): Number of CPUs to use.

After setting up the **config.yaml**, you can run the script by typing **python3.x run_gem_recon.py**. Please choose the python version compatible with IBM CPLEX solver.


