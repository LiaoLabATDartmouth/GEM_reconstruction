import os
import random
import re
import pandas as pd
import numpy as np
import cobra
from cobra.flux_analysis import flux_variability_analysis
from cobra.manipulation import remove_genes
from cobra import Model, Reaction
from optlang.interface import OPTIMAL
from optlang.symbolics import add
import straindesign as sd
import yaml


#######################
# Load parameter values
#######################

# default values for optional parameters
DEFAULT_CONFIG = {
    "output_folder": "Output",
    "protein_identity_score_cutoff": 95,
    "gapfilling_method": "mcs",  # Available options: "minrxnnum", "maxavescore", "mcs"
    "n_cpus": 1
}

# check if config.yaml exists
config_file = "config.yaml"
if not os.path.exists(config_file):
    raise FileNotFoundError(f"Configuration file '{config_file}' not found. Please create it with the required parameters.")
else:
    # load YAML file
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    # essential parameters that must be present
    required_keys = ["pan_gem_protein_sequence_file", "pan_gem_gene_id_mapping_file", "pan_gem_model_file", "query_protein_sequence_files"]

    # check for missing required parameters
    missing_keys = [key for key in required_keys if key not in config]
    if missing_keys:
        raise ValueError(f"Missing required configuration parameters: {', '.join(missing_keys)}")

    # load required parameters
    pan_gem_protein_sequence_file = config["pan_gem_protein_sequence_file"]
    pan_gem_gene_id_mapping_file = config["pan_gem_gene_id_mapping_file"]
    pan_gem_model_file = config["pan_gem_model_file"]
    query_protein_sequence_files = config["query_protein_sequence_files"]

    # load optional parameters with defaults
    output_folder = config.get("output_folder", DEFAULT_CONFIG["output_folder"])
    protein_identity_score_cutoff = config.get("protein_identity_score_cutoff", DEFAULT_CONFIG["protein_identity_score_cutoff"])
    gapfilling_method = config.get("gapfilling_method", DEFAULT_CONFIG["gapfilling_method"])
    n_cpus = config.get("n_cpus", DEFAULT_CONFIG["n_cpus"])


################################
# Functions defined in this file
################################


# set environment to glucose minimal media
def set_glucose_minimal_media(_model):
    # reset flux bounds to be between -1000 and 1000
    for rxn in _model.reactions:
        if rxn.lower_bound < -1000.0:
            rxn.lower_bound = -1000.0
        if rxn.upper_bound > 1000.0:
            rxn.upper_bound = 1000.0

    # synthetic minimal media recipe
    # 20 g/L glucose
    # 2.7 g/L ammonium sulphate
    # 0.05 g/L magnesium sulphate
    # 2 g/L potassium dihydrogen phosphate
    # 0.5 g/L calcium chloride
    # 100  g/L biotin (Sigma).

    # open medium component influx
    for reaction in _model.exchanges:
        if reaction.id.startswith('EX_'):
            if reaction.id in [
                'EX_C00031__extr',    # D-glucose
                'EX_C00014__extr',    # ammonium
                'EX_Oxygen__extr',    # O2
                'EX_C00009__extr',    # Phosphate
                'EX_C00059__extr',    # Sulfate
                'EX_C00120__extr'     # biotin
            ]:
                if reaction.id == 'EX_C00031__extr':
                    reaction.lower_bound = -2.0
                else:
                    reaction.lower_bound = -1000.0
            else:
                reaction.lower_bound = 0.0
    return _model


# get reaction-level confidence score by replacing genes in gene reaction rules with their protein similarity scores
def or2max_and2min(_gene_reaction_rule, _protein_similarity_scores):
    string2parse = '(' + _gene_reaction_rule + ')'
    parens = [m for m in re.finditer(r"(\([^(]*?\))", string2parse)]  # find the innermost parenthesis
    while len(parens) > 0:  # as long as a parenthesis can be found
        for paren in parens:
            extract = paren.group()
            content = extract[1:-1]  # remove "(" and ")"
            if 'and' in content:  # find the gene with minimum expression value
                assert 'or' not in content
                genes = [g.strip() for g in content.split('and')]
                scores = [_protein_similarity_scores[g] for g in genes]
                substitute = genes[np.argmin(scores)]
            elif 'or' in content:  # find the gene with maximum expression value
                assert 'and' not in content
                genes = [g.strip() for g in content.split('or')]
                scores = [_protein_similarity_scores[g] for g in genes]
                substitute = genes[np.argmax(scores)]
            else:
                substitute = content
            string2parse = string2parse.replace(extract, substitute)  # find the innermost parenthesis
        parens = [m for m in re.finditer(r"(\([^(]*?\))", string2parse)]
    return string2parse, _protein_similarity_scores[string2parse]


# solve flux balance model
def find_optimal_solution(_model, _solver_type, _mip_tol_int):
    # set up solver parameters
    # increase tolerance if returned primal are not integers
    _model.solver.problem.parameters.mip.tolerances.integrality.set(_mip_tol_int)

    # available methods implemented by CPLEX
    if _solver_type == 'lp':
        available_methods = ["auto", "primal", "dual", "network", "barrier", "sifting", "concurrent"]
    elif _solver_type == "qp":
        available_methods = ["auto", "primal", "dual", "network", "barrier"]  # sifting and concurrent are not valid qp_method
    else:
        print("Unknown optimization solver_type: %s" % _solver_type)
        return None, 0

    # Note that model.optimize will still raise error if solver status is None or does not allow retrieve primal values
    # the following solver status has primal values: NUMERIC, FEASIBLE, INFEASIBLE, SUBOPTIMAL, ITERATION_LIMIT, TIME_LIMIT
    # i.e., if the solver status is check_original_solver_status, it will throw out error regardless
    # CPXMIP_OPTIMAL_INFEAS is a CPLEX status that indicates that a preprocessed problem was solved to optimality
    # but the solver was unable to recover a feasible solution to the original problem.
    solution_found = 0
    solution = None
    for method in available_methods:
        if _solver_type == 'lp':
            _model.solver.configuration.lp_method = method
        elif _solver_type == 'qp':
            _model.solver.configuration.qp_method = method

        try:
            solution = _model.optimize()
            if _model.solver.status != OPTIMAL:
                # print("%s fails: solver status = %s"%(method, model.solver.problem.solution.get_status()))
                solution_found = 0
            else:
                solution_found = 1
                break
        except:
            # print("%s fails: solver status = %s"%(method, model.solver.problem.solution.get_status()))
            solution_found = 0

    if solution_found == 1:
        return solution, solution_found
    else:
        return None, solution_found


# prune model by removing reactions
def prune_model(_model_id, _reaction_ids_to_keep, _protein_besthit_matrix=None):
    pruned_model = Model(_model_id)
    pruned_model.name = pruned_model.id
    pruned_model.solver = 'cplex'

    # add reactions from pan model
    reactions_to_keep = []
    pan_gem_model = cobra.io.read_sbml_model(pan_gem_model_file)
    for reaction_id in list(_reaction_ids_to_keep):
        if reaction_id in pan_gem_model.reactions:
            reactions_to_keep.append(pan_gem_model.reactions.get_by_id(reaction_id))
    pruned_model.add_reactions(reactions_to_keep)
    fix_compartment_information(pruned_model, pan_gem_model)

    # remove genes not associated with any reaction
    genes_to_remove = []
    for gene in pruned_model.genes:
        if len(gene.reactions)==0:
            genes_to_remove.append(gene.id)
    remove_genes(pruned_model, genes_to_remove)

    # update gene ID, name, and GPR
    if _protein_besthit_matrix is not None:
        # update GPR
        _protein_besthit_dict = dict(_protein_besthit_matrix.loc[_protein_besthit_matrix[_model_id] != '', _model_id])
        keys = (re.escape(k) for k in _protein_besthit_dict.keys())
        pattern = re.compile(r'\b(' + '|'.join(keys) + r')\b')
        for rxn in pruned_model.reactions:
            gpr = rxn.gene_reaction_rule
            rxn.gene_reaction_rule = pattern.sub(lambda x: _protein_besthit_dict[x.group()], gpr)

        # update genes
        genes_to_remove = []
        for gene in pruned_model.genes:
            if len(gene.reactions)==0:
                genes_to_remove.append(gene.id)
        remove_genes(pruned_model, genes_to_remove)

        # reassign gene name
        for gene in pruned_model.genes:
            gene.name = gene.id

    # simulate growth
    pruned_model = set_glucose_minimal_media(pruned_model)
    biomass_reaction = pruned_model.reactions.get_by_id('e_Biomass__cyto')
    biomass_reaction.bounds = (0, 1000.0)
    pruned_model.objective = 'e_Biomass__cyto'
    return pruned_model


# gapfilling model by adding minimum number of reactions
def gapfill_minrxnnum(_model, _candidate_reactions, _reaction_confidence_scores):
    # add candidate reactions to the model
    _model.add_reactions(list(_candidate_reactions.values()))

    # set minimal growth
    biomass_reaction = _model.reactions.get_by_id('e_Biomass__cyto')
    biomass_reaction.bounds = (0.01, 0.01)

    # set constraints: flux is unconstrained when indicator variable is 1, otherwise 0
    indicator_vars = []
    indicator_vars_weighted = []
    for reaction_id, reaction in _candidate_reactions.items():
        var = _model.problem.Variable('indicator_var__' + reaction_id, lb=0, ub=1, type='binary')
        indicator_vars.append(var)
        indicator_vars_weighted.append(var * _reaction_confidence_scores[reaction_id])
        constrain_lower_bound = _model.problem.Constraint(
            (reaction.flux_expression - reaction.lower_bound * var).expand(),
            name='constrain_reaction_lower_bound__' + reaction_id,
            lb=0
        )
        constrain_upper_bound = _model.problem.Constraint(
            (reaction.flux_expression - reaction.upper_bound * var).expand(),
            name='constrain_reaction_upper_bound__' + reaction_id,
            ub=0
        )
        _model.add_cons_vars([var, constrain_lower_bound, constrain_upper_bound])

    _model.objective = add(*indicator_vars)
    _model.objective.direction = "min"
    _model.solver.update()

    # minimize number of reactions needed for gapfilling
    sol_min, success_found_min = find_optimal_solution(_model, _solver_type='lp', _mip_tol_int=1e-9)
    if success_found_min == 0:
        biomass_reaction.bounds = (0, 1000.0)
        return []
    else:
        minimum_reactions_added = 0
        gapfilled_reaction_ids = []
        for var in _model.variables:
            if var.name.startswith('indicator_var__') and var.primal == 1:
                minimum_reactions_added += 1
                gapfilled_reaction_ids.append(var.name.split('__')[1])
        assert np.round(sol_min.objective_value) == minimum_reactions_added

        # fix the number of active reactions and maximize the total reaction scores
        constrain_min_num_reactions = _model.problem.Constraint(
            add(*indicator_vars),
            name='constrain_total_reaction_number_added',
            ub=minimum_reactions_added
        )
        _model.add_cons_vars([constrain_min_num_reactions])
        _model.objective = add(*indicator_vars_weighted)
        _model.objective.direction = "max"
        _model.solver.update()

        sol_minmax, solution_found_minmax = find_optimal_solution(_model, _solver_type='lp', _mip_tol_int=1e-9)
        if solution_found_minmax == 0:
            biomass_reaction.bounds = (0, 1000.0)
            return gapfilled_reaction_ids
        else:
            # rewrite gapfilled_reaction_ids using new solution
            gapfilled_reaction_ids = []
            for var in _model.variables:
                if var.name.startswith('indicator_var__') and var.primal == 1:
                    gapfilled_reaction_ids.append(var.name.split('__')[1])
            assert len(gapfilled_reaction_ids) == minimum_reactions_added
            biomass_reaction.bounds = (0, 1000.0)
            return gapfilled_reaction_ids


# gapfilling model by adding reactions with maximum average score
def gapfill_maxavescore(_model, _candidate_reactions, _reaction_confidence_scores):
    # add candidate reactions to the model
    _model.add_reactions(list(_candidate_reactions.values()))

    # set minimal growth
    biomass_reaction = _model.reactions.get_by_id('e_Biomass__cyto')
    biomass_reaction.bounds = (0.01, 0.01)

    # set constraints: average flux = sum_i(reaction_i_flux * indicator_variable_i)/sum(indicator_variable_i)
    # see how to linearize the problem here
    # https://math.stackexchange.com/questions/2752558/how-to-linearize-a-weighted-average-with-a-decision-variable
    wx_vars = []
    z_vars = []
    M = max(_reaction_confidence_scores.values())
    y = _model.problem.Variable('y', lb=0, ub=M, type='continuous')
    for reaction_id, reaction in _candidate_reactions.items():
        xi = _model.problem.Variable('x__' + reaction_id, lb=0, ub=1, type='binary')
        constrain_reaction_lower_bound = _model.problem.Constraint(
            (reaction.flux_expression - reaction.lower_bound * xi).expand(),
            name='constrain_reaction_lower_bound__' + reaction_id,
            lb=0
        )
        constrain_reaction_upper_bound = _model.problem.Constraint(
            (reaction.flux_expression - reaction.upper_bound * xi).expand(),
            name='constrain_reaction_upper_bound__' + reaction_id,
            ub=0
        )
        zi = _model.problem.Variable('z__' + reaction_id, lb=0, ub=M, type='continuous')
        # zi <= xi*max(wi)
        constrain_zi_le_xiM = _model.problem.Constraint(zi - xi * M, name='constrain_zi_le_xiM__' + reaction_id, ub=0)
        # zi <= y
        constrain_zi_le_y = _model.problem.Constraint(zi - y, name='constrain_zi_le_y__' + reaction_id, ub=0)
        # zi >= y-(1-xi)*max(wi)
        constrain_zi_ge_y_minus_M_plus_xiM = _model.problem.Constraint(
            zi - y + M - M * xi,
            name='constrain_zi_ge_y_minus_M_plus_xiM__' + reaction_id,
            lb=0
        )
        wx_vars.append(xi * _reaction_confidence_scores[reaction_id])
        z_vars.append(zi)
        _model.add_cons_vars(
            [xi,
             constrain_reaction_lower_bound,
             constrain_reaction_upper_bound,
             zi,
             constrain_zi_le_xiM,
             constrain_zi_le_y,
             constrain_zi_ge_y_minus_M_plus_xiM
             ]
        )
    # sum(zi) = sum(wi*xi)
    constrain_sum_zi_eq_sum_wi_xi = _model.problem.Constraint(
        add(*z_vars) - add(*wx_vars),
        name='constrain_sum_zi_eq_sum_wi_xi',
        lb=0,
        ub=0
    )
    _model.add_cons_vars([constrain_sum_zi_eq_sum_wi_xi])

    _model.objective = y
    _model.objective.direction = "max"
    _model.solver.update()

    solution, solution_found = find_optimal_solution(_model, _solver_type='lp', _mip_tol_int=1e-9)
    if solution_found == 0:
        biomass_reaction.bounds = (0, 1000.0)
        return []
    else:
        gapfilled_reaction_ids = []
        for var in _model.variables:
            if var.name.startswith('x__') and var.primal == 1:
                gapfilled_reaction_ids.append(var.name.split('__')[1])
        total_reaction_score = sum([_reaction_confidence_scores[reaction_id] for reaction_id in gapfilled_reaction_ids])
        average_reaction_score = total_reaction_score / len(gapfilled_reaction_ids)
        assert abs(solution.objective_value - average_reaction_score) < 1e-10
        biomass_reaction.bounds = (0, 1000.0)
        return gapfilled_reaction_ids


# assign missing compartment to metabolites in a model
def fix_compartment_information(_model, _original_model):
    _model.compartments = _original_model.compartments
    for metabolite in _model.metabolites:
        if metabolite.compartment is None:
            metabolite.compartment = _original_model.metabolites.get_by_id(metabolite.id).compartment


if __name__ == '__main__':
    ################################################################
    # Blast pan gem protein sequences against input genome sequences
    ################################################################
    print("Blast start...")
    for seq_file in os.listdir(query_protein_sequence_files):
        if seq_file.endswith(".fasta") or seq_file.endswith(".faa") or seq_file.endswith(".fa"):
            # get sample name and suffix of file name
            suffix = seq_file.split('.')[-1]
            sample_name = seq_file.replace("." + suffix, "")

            # create an output folder for each genome
            if not os.path.exists("%s/%s" % (output_folder, sample_name)):
                os.system("mkdir -p %s/%s" % (output_folder, sample_name))

            # create a subfolder for blast database and run blastp
            if not os.path.exists("%s/%s/blast_db" % (output_folder, sample_name)):
                os.system("mkdir -p %s/%s/blast_db" % (output_folder, sample_name))
                os.system("makeblastdb -in %s/%s -title %s -dbtype prot -out %s/%s/blast_db/%s" % (
                    query_protein_sequence_files, seq_file, sample_name, output_folder, sample_name, sample_name
                ))

                print("Blast pan gem protein sequences against %s..." % sample_name)
                os.system("blastp -db %s/%s/blast_db/%s -query %s -out %s/%s/%s -outfmt \"6 qseqid qlen sseqid slen evalue bitscore length nident pident\" -num_threads %d" % (
                    output_folder, sample_name, sample_name,
                    pan_gem_protein_sequence_file,
                    output_folder, sample_name, sample_name + ".blast",
                    n_cpus
                ))
    print("Blast done.")

    #######################################################
    # Parse blast output and generate identity score matrix
    #######################################################
    protein_identity_matrix = pd.read_csv(pan_gem_gene_id_mapping_file)[['PanGemGeneID', 'QueryGeneID']]
    protein_besthit_matrix = pd.read_csv(pan_gem_gene_id_mapping_file)[['PanGemGeneID', 'QueryGeneID']]
    for sample_name in os.listdir(output_folder):
        # only iterate over existing folders
        if os.path.isdir(os.path.join(output_folder, sample_name)):
            df_blast = pd.read_csv("%s/%s/%s.blast" % (output_folder, sample_name, sample_name), sep="\t", header=None)
            df_blast.columns = ["qseqid", "qlen", "sseqid", "slen", "evalue", "bitscore", "length", "nident", "pident"]
            df_blast["nident"] = df_blast["nident"].astype(float)
            df_blast["score"] = df_blast["nident"] / df_blast["qlen"] * 100

            # protein identity matrix
            df_blast_iden = df_blast[["qseqid", "score"]].sort_values(["qseqid", "score"], ascending=False)
            df_blast_iden = df_blast_iden.drop_duplicates(subset="qseqid", keep="first")
            df_blast_iden = df_blast_iden.rename({"qseqid": "QueryGeneID", "score": sample_name}, axis=1)
            df_blast_iden['QueryGeneID'] = [_id.split('.')[0] for _id in df_blast_iden['QueryGeneID']]
            protein_identity_matrix = pd.merge(
                protein_identity_matrix,
                df_blast_iden,
                left_on="QueryGeneID",
                right_on="QueryGeneID",
                how="left"
            ).fillna(0)

            # protein besthit matrix
            df_blast_bh = df_blast[["qseqid", "sseqid", "score"]].sort_values(["qseqid", "score"], ascending=False)
            df_blast_bh = df_blast_bh.drop_duplicates(subset="qseqid", keep="first").drop('score', axis=1)
            df_blast_bh = df_blast_bh.rename({"qseqid": "QueryGeneID", "sseqid": sample_name}, axis=1)
            df_blast_bh['QueryGeneID'] = [_id.split('.')[0] for _id in df_blast_bh['QueryGeneID']]
            protein_besthit_matrix = pd.merge(
                protein_besthit_matrix,
                df_blast_bh,
                left_on="QueryGeneID",
                right_on="QueryGeneID",
                how="left"
            ).fillna('')
    protein_identity_matrix = protein_identity_matrix.drop("QueryGeneID", axis=1).set_index('PanGemGeneID')
    protein_identity_matrix = protein_identity_matrix.round(3)
    protein_identity_matrix.to_csv("%s/protein_identity_score_matrix.csv" % output_folder)

    # binarize identity matrix by using identity score cutoff
    binary_protein_identity_matrix = protein_identity_matrix.copy()
    binary_protein_identity_matrix[protein_identity_matrix >= protein_identity_score_cutoff] = 1
    binary_protein_identity_matrix[protein_identity_matrix < protein_identity_score_cutoff] = 0
    binary_protein_identity_matrix = binary_protein_identity_matrix.astype(int)
    binary_protein_identity_matrix.to_csv("%s/binary_protein_identity_score_matrix.csv" % output_folder)

    print("Construct protein identity matrix. Done.")

    protein_besthit_matrix = protein_besthit_matrix.drop("QueryGeneID", axis=1).set_index('PanGemGeneID')
    for sample_name in binary_protein_identity_matrix.columns:
        protein_besthit_matrix.loc[list(binary_protein_identity_matrix[binary_protein_identity_matrix[sample_name] == 0].index), sample_name] = ''
    protein_besthit_matrix.to_csv("%s/protein_besthit_sequence_matrix.csv" % output_folder)

    ####################
    # Load Pan GEM model
    ####################
    pan_gem_model = cobra.io.read_sbml_model(pan_gem_model_file)
    pan_gem_model.solver = 'cplex'
    pan_gem_model = set_glucose_minimal_media(pan_gem_model)
    pan_gem_model.objective = 'e_Biomass__cyto'  # set objective function (biomass)
    assert pan_gem_model.slim_optimize() >= 0.01
    print("Load Pan GEM model. Done.")

    ########################################################################
    # Convert gene presence/absence table to reaction presence/absence table
    ########################################################################
    reaction_presence_absence_matrix = pd.DataFrame(
        index=[rxn.id for rxn in pan_gem_model.reactions],
        columns=protein_identity_matrix.columns
    )
    reaction_confidence_matrix = pd.DataFrame(
        index=[rxn.id for rxn in pan_gem_model.reactions],
        columns=protein_identity_matrix.columns
    )
    for rid in reaction_presence_absence_matrix.index:
        rxn = pan_gem_model.reactions.get_by_id(rid)
        grr = str(rxn.gene_reaction_rule)
        if grr == '':
            # non genome-associated reactions are included by default
            reaction_presence_absence_matrix.loc[rid, :] = 1
            reaction_confidence_matrix.loc[rid, :] = 100.0
        else:
            for sample_name in reaction_presence_absence_matrix.columns:
                # reaction presence/absence
                grr2 = grr
                gene_values = {}
                while 'CPAR2' in grr2:
                    CPAR2_loc = grr2.find('CPAR2')
                    selected_gene = grr2[CPAR2_loc:CPAR2_loc + 12]
                    selected_gene_value = protein_identity_matrix.loc[selected_gene, sample_name]

                    if isinstance(selected_gene_value, float):
                        gene_values[selected_gene] = protein_identity_matrix.loc[selected_gene, sample_name]
                    else:
                        raise Exception("Wrong data type (%s): a float is expected. Check gene ID mapping file for duplicate KEGG IDs."%(type(selected_gene_value)))

                    binary_gene_value = binary_protein_identity_matrix.loc[selected_gene, sample_name]
                    if isinstance(binary_gene_value, (int, np.int64)):
                        grr2 = grr2.replace(selected_gene, str(binary_gene_value))
                    else:
                        raise Exception("Wrong data type (%s): an integer is expected. Check gene ID mapping file for duplicate KEGG IDs."%(type(binary_gene_value)))
                reaction_presence_absence_matrix.loc[rid, sample_name] = eval(grr2)

                # reaction confidence score
                _, reaction_confidence_matrix.loc[rid, sample_name] = or2max_and2min(grr, gene_values)

    reaction_presence_absence_matrix.to_csv("%s/reaction_presence_absence_matrix.csv" % output_folder)
    reaction_confidence_matrix.to_csv("%s/reaction_confidence_matrix.csv" % output_folder)

    print("Calculate reaction confidence score. Done.")

    #############################################
    # Find essential reactions and the active set
    #############################################
    if not os.path.exists("%s/essential_and_active_reaction_set.csv" % output_folder):
        biomass_reaction = pan_gem_model.reactions.get_by_id('e_Biomass__cyto')
        biomass_reaction.bounds = (0.01, 0.01)  # allow for slow growth
        fva = flux_variability_analysis(
            model=pan_gem_model,
            reaction_list=pan_gem_model.reactions,
            processes=n_cpus
        )

        #-----------------------------------------------------------------------------------------
        # the minimum cutoff flux used in FVA (current value: 1e-6) might matter
        #
        # if it is set too high, essential reactions may be filtered out.
        # this will lead to no growth of a model made of reactions solely from the active set.
        #
        # if it is set too low, inactive reactions may be unintentionally captured.
        # this does not do too much harm as long as the number of such false positives remain low.
        #-----------------------------------------------------------------------------------------

        # essential reactions (no growth if these reactions are deleted)
        essential_set = set(fva.query('minimum*maximum>0.000000000001').index)
        print("Found %d reactions in the essential set."%len(essential_set))
        # reactions that are not dead end (can carry flux in some condition)
        active_set = set(fva.query('abs(minimum)>0.000001 or abs(maximum)>0.000001').index)
        print("Found %d reactions in the active set."%len(active_set))
        pan_gem_model_active = None  # a model with only active reactions (will be updated when needed)
        biomass_reaction.bounds = (0, 1000.0)

        df_essential_set = pd.DataFrame([list(essential_set), [1]*len(essential_set)], index=["Reaction","Essential"]).T
        df_active_set = pd.DataFrame([list(active_set), [1]*len(active_set)], index=["Reaction","Active"]).T
        df_outjoin_set = pd.merge(df_essential_set, df_active_set, left_on="Reaction", right_on="Reaction", how="outer").fillna(0)
        df_outjoin_set.to_csv("%s/essential_and_active_reaction_set.csv"% output_folder, index=False)
    else:
        df_outjoin_set = pd.read_csv("%s/essential_and_active_reaction_set.csv"% output_folder)
        essential_set = set(df_outjoin_set[df_outjoin_set.Essential==1].Reaction)
        active_set = set(df_outjoin_set[df_outjoin_set.Active==1].Reaction)
        pan_gem_model_active = None  # a model with only active reactions (will be updated when needed)

    # the essential reactions will be added to the model regardless of their enzyme identity scores
    for reaction_id in essential_set:
        reaction_presence_absence_matrix.loc[reaction_id, :] = 1
    print("Search for essential and active reactions. Done.")

    ##############################################
    # Construct strain-specific GEM and gapfilling
    ##############################################
    gapfilling_summary = []
    for sample_name in protein_identity_matrix.columns:
        reaction_ids_to_keep = list(
            reaction_presence_absence_matrix[reaction_presence_absence_matrix[sample_name] == 1].index
        )
        draft_model = prune_model(
            sample_name,
            reaction_ids_to_keep,
            protein_besthit_matrix
        )
        if draft_model.slim_optimize() >= 0.01:
            # write to file
            cobra.io.write_sbml_model(draft_model, "%s/%s/%s.xml" % (output_folder, sample_name, sample_name))
            print("Generate a draft model for %s. No gapfilling needed. Save to file." % sample_name)
        else:
            # gapfilling starts
            if gapfilling_method == "minrxnnum" or gapfilling_method == "maxavescore":
                print("Gapfilling (method: %s) for %s start..." % (gapfilling_method, sample_name))
                # random seed in cplex can affect model convergence
                # try different random seed until it fails 20 times
                count_of_failure = 0
                while count_of_failure < 20:
                    # print(count_of_failure)
                    reaction_ids_to_keep = list(
                        reaction_presence_absence_matrix[reaction_presence_absence_matrix[sample_name] == 1].index
                    )
                    print("number of reactions to keep = %d" % len(reaction_ids_to_keep))
                    draft_model = prune_model(
                        sample_name,
                        reaction_ids_to_keep,
                        protein_besthit_matrix
                    )
                    draft_model.solver.problem.parameters.randomseed.set(random.randint(0, 99))
                    candidate_reaction_ids = list(
                        active_set.intersection(
                            set(reaction_presence_absence_matrix[reaction_presence_absence_matrix[sample_name] == 0].index)
                        )
                    )
                    print("number of reactions in candidate pool = %d" % len(candidate_reaction_ids))
                    candidate_reactions = {
                        reaction_id: pan_gem_model.reactions.get_by_id(reaction_id) for reaction_id in candidate_reaction_ids
                    }
                    reaction_confidence_scores = {
                        reaction_id: reaction_confidence_matrix.loc[reaction_id, sample_name] for reaction_id in candidate_reaction_ids
                    }

                    gapfilled_reaction_ids = []
                    if gapfilling_method == "minrxnnum":
                        gapfilled_reaction_ids = gapfill_minrxnnum(
                            draft_model, candidate_reactions, reaction_confidence_scores
                        )
                    elif gapfilling_method == "maxavescore":
                        gapfilled_reaction_ids = gapfill_maxavescore(
                            draft_model, candidate_reactions, reaction_confidence_scores
                        )
                    if len(gapfilled_reaction_ids) > 0:
                        num_reactions_added = len(gapfilled_reaction_ids)
                        average_confidence_score = sum(
                            [reaction_confidence_scores[reaction_id] for reaction_id in gapfilled_reaction_ids]
                        ) / num_reactions_added
                        gapfilling_summary.append([
                            sample_name,
                            gapfilling_method,
                            protein_identity_score_cutoff,
                            len(gapfilled_reaction_ids),
                            average_confidence_score,
                            ','.join(gapfilled_reaction_ids)
                        ])
                        gapfilled_model = prune_model(
                            sample_name,
                            reaction_ids_to_keep + gapfilled_reaction_ids,
                            protein_besthit_matrix
                        )
                        assert gapfilled_model.slim_optimize() >= 0.01
                        cobra.io.write_sbml_model(
                            gapfilled_model,
                            "%s/%s/%s_gf_%s.xml" % (output_folder, sample_name, sample_name, gapfilling_method)
                        )
                        print("Gapfilling done. %d reactions added. Save to file." % num_reactions_added)
                        break
                    else:
                        count_of_failure += 1

                if count_of_failure == 20:
                    print("Gapfilling failed.")
            elif gapfilling_method == "mcs":
                print("Gapfilling (method: %s) for %s start..." % (sample_name, gapfilling_method))

                # A pruned model with only active reactions is used for finding minimal cut set
                # running straindesign with original model might be computationally expensive
                if pan_gem_model_active is None:
                    pan_gem_model_active = prune_model(
                        "pan_gem_model_active",
                        list(active_set)
                    )
                    assert pan_gem_model_active.slim_optimize() >= 0.01
                reaction_ids_to_keep = list(
                    reaction_presence_absence_matrix[reaction_presence_absence_matrix[sample_name] == 1].index
                )
                candidate_reaction_ids = list(
                    active_set.intersection(
                        set(reaction_presence_absence_matrix[
                                reaction_presence_absence_matrix[sample_name] == 0].index)
                    )
                )
                reaction_confidence_scores = {
                    reaction_id: reaction_confidence_matrix.loc[reaction_id, sample_name] for reaction_id in
                    candidate_reaction_ids
                }
                ko_cost = {reaction_id: 1 for reaction_id in candidate_reaction_ids}
                modules = [sd.SDModule(pan_gem_model_active, sd.names.SUPPRESS, constraints='e_Biomass__cyto >= 0.01')]
                solutions = sd.compute_strain_designs(
                    pan_gem_model_active,
                    sd_modules=modules,
                    solution_approach='any',
                    ko_cost=ko_cost
                )
                # rescue one reaction (the one with highest confidence score) from each mcs
                gapfilled_reaction_ids = []
                for solution in solutions.reaction_sd:
                    cset = list(solution)
                    cset_scores = [reaction_confidence_matrix.loc[reaction_id, sample_name] for reaction_id in cset]
                    gapfilled_reaction_ids.append(cset[cset_scores.index(max(cset_scores))])
                gapfilled_reaction_ids = list(set(gapfilled_reaction_ids))

                if len(solutions.reaction_sd) > 0:
                    num_reactions_added = len(gapfilled_reaction_ids)
                    average_confidence_score = sum(
                        [reaction_confidence_scores[reaction_id] for reaction_id in gapfilled_reaction_ids]
                    ) / num_reactions_added
                    gapfilling_summary.append([
                        sample_name,
                        gapfilling_method,
                        protein_identity_score_cutoff,
                        len(gapfilled_reaction_ids),
                        average_confidence_score,
                        ','.join(gapfilled_reaction_ids)
                    ])
                    gapfilled_model = prune_model(
                        sample_name,
                        reaction_ids_to_keep + gapfilled_reaction_ids,
                        protein_besthit_matrix
                    )
                    assert gapfilled_model.slim_optimize() >= 0.01
                    cobra.io.write_sbml_model(
                        gapfilled_model,
                        "%s/%s/%s_gf_%s.xml" % (output_folder, sample_name, sample_name, gapfilling_method)
                    )
                    print("Gapfilling done. %d reactions added. Save to file." % num_reactions_added)
                else:
                    print("Gapfilling failed.")


    if len(gapfilling_summary) > 0:
        df_gapfilling_summary = pd.DataFrame(
            gapfilling_summary,
            columns=["Sample","Method","Cutoff","NumRxnAdded","AvgScore","GapfilledRxns"]
        )
        if not os.path.exists("%s/gapfilling_summary.csv" % output_folder):
            df_gapfilling_summary.to_csv("%s/gapfilling_summary.csv"% output_folder, index=False)
        else:
            df_summary_combine = pd.read_csv("%s/gapfilling_summary.csv" % output_folder)
            df_summary_combine = pd.concat([df_summary_combine, df_gapfilling_summary], ignore_index=True)
            df_summary_combine.to_csv("%s/gapfilling_summary.csv"% output_folder, index=False)

