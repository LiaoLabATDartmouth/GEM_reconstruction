# GSMM_reconstruction

The script "run_gsmm_recon.py" reconstruct strain-specific genome-scale metabolic models from their whole-genome sequencing data. Before running the script, please set up the following parameters in the configuration file "config.yaml":

### required parameters
* pan_gem_model_file # genome-scale metabolic model of the universe (pan-genome) model
* pan_gem_gene_id_mapping_file # a two-column table that maps gene IDs in the pan-genome model (PanGemGeneID) to those in the query genomes (QueryGeneID)
* pan_gem_protein_sequence_file # protein fasta file of the pan-genome model
* query_protein_sequence_files # protein fasta files of the query genomes (one fasta file per query genome)

### optinal parameters
* output_folder (default: "Output") # data output directory (each query genome has its own subfolder under this directory)
* protein_identity_score_cutoff (default: 95) # cutoff value used to determine whether a protein in the query genome is present in the pan-genome model
* gapfilling_method (default: "mcs") # available options include "minrxnnum", "maxavescore", and "mcs"
* n_cpus (default: 1): number of cpus

