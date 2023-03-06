# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 11:00:09 2023

@author: danie
"""

import os, sys
import numpy as np
import cobra
cobra_config = cobra.Configuration()
cobra_config.solver = 'glpk'

from pathlib import Path
path = Path.cwd()

from modelseedpy import MSBuilder, MSGenome
from modelseedpy.core import msmedia
from modelseedpy.core.rast_client import RastClient

sys.path.append(os.path.join(path.parents[1], 'scripts'))


import build_model, gapfill_function, NN_Predictor
from   reaction_class import Reaction
""" 
To install modelseedpy
Follow instructions from: https://github.com/ModelSEED/ModelSEEDpy 
using the git clone option
"""
from modelseedpy import MSBuilder, MSGenome
from modelseedpy.core import msmedia


"""
To annotate a genome with RAST (the following commands are to be executed on terminal): 
    1. Installation 
        Following the instructions here: https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html we get RAST locally:
        curl -O -L https://github.com/BV-BRC/BV-BRC-CLI/releases/download/1.040/bvbrc-cli-1.040.deb
        sudo dpkg -i bvbrc-cli-1.040.deb 
        sudo apt-get -f install
    
    2. Build a new genome object
    The genetic code used in translating to protein sequences. See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more information on genetic codes.
    rast-create-genome --scientific-name "Salmonella infantis" --genetic-code 11 --domain Bacteria --contigs *.fna > genome.gto
    3. Run the annotation step
    
    rast-process-genome < *.gto > annotated_genome.gto2
    4. Export as a .faa file 
    
    rast-export-genome protein_fasta < annotated_genome.gto2 > annotated_genome.faa
This script will use the annotated_genome.faa output of the RAST annotation as its input
"""

# Set the path to your genome
print("Build MSGenome object")
patric_genome = MSGenome.from_fasta(os.path.join(path.parents[1], 'files', 'genomes', 'bh_ensembl_cDNA.fa'), split = ' ')

rast = RastClient()
rast.annotate_genome(patric_genome)

print("Start building base model")
base_model = MSBuilder.build_metabolic_model(model_id = "Blautia hydrogenotrophica DSM 10507", 
                                             genome   = patric_genome, 
                                             index    = "0",
                                             classic_biomass = True, 
                                             gapfill_model   = False, 
                                             gapfill_media   = None, 
                                             annotate_with_rast = False,
                                             allow_all_non_grp_reactions = True
                                            )


model_name = "bh_ungapfilled_model.sbml"
cobra.io.write_sbml_model(cobra_model = base_model, filename = os.path.join(path.parents[0], 'models', model_name))

print("An ungapfilled modelSEED reconstructed model in now available.")

# Set paths to 
files_path  = os.path.join(path.parents[0])
models_path = os.path.join(files_path, 'models')
path_to_NN  = os.path.join(files_path, 'NN')
path_to_biochem  = os.path.join(files_path,  'biochemistry', 'reactions.tsv')

# Load NN predictor
trainedNN_MS = NN_Predictor.load_NN( path = os.path.join(path_to_NN, 'NN_MS.h5') )



draft_model = cobra.io.read_sbml_model(os.path.join(models_path, 'bh_ungapfilled_model.sbml')) 

# Build a Reaction object for the reactions present on your model
# draft_reaction = Reaction( model = os.path.join(models_path, 's_infantis_base_model.sbml') )
draft_reaction = Reaction( model = os.path.join(models_path, 'bh_ungapfilled_model.sbml') )

 
                 



# Here goes your medium if you want to have a defined one, otherwaise it is set to None
#def_med = {}

#def_med["EX_cpd00027_e0"] = {'lower_bound': -10., 'upper_bound': 10., 'metabolites': {"cpd00027_e0":-1.}}
#def_med["EX_cpd00001_e0"] = {'lower_bound': -10., 'upper_bound': 10., 'metabolites': {"cpd00001_e0":-1.}}



# Build a Reaction object for the exchange reactions; if you have a defined medium, set the fixed_bounds argument accordingly
exchange_reacs = Reaction( model        = os.path.join(models_path, 'exchangeReactions.sbml'), 
                           fixed_bounds = None)



# Build a Reaction object for the exchange reactions; if you have a defined medium, set the fixed_bounds argument accordingly

#if using a defined media, change to
# exchange_reacs = Reaction( model        = os.path.join(models_path, 'exchangeReactions.sbml'), 
#                            fixed_bounds = def_med)



# """ Consider adding this under the gapfill() function """
# for react in exchange_reacs.reactions:
#     if react not in def_med:
#         exchange_reacs.reactions[react]["lower_bound"] = 0





        
# Build a Reaction object for the reactions under the biochemistry folder 
db_reactions           = Reaction(biochem_input = path_to_biochem)
db_reactions.reactions = db_reactions.add_dict(exchange_reacs.reactions, db_reactions.reactions)

#Combine all reactions into one object
all_reactions           = Reaction(fixed_bounds = None)
#if using a defined media, change to
#all_reactions           = Reaction(fixed_bounds = def_media)

all_reactions.reactions = all_reactions.add_dict(draft_reaction.reactions, db_reactions.reactions)

draft_reaction_ids      = set(draft_reaction.reactions)


""" Consider adding this under the gapfill() function """
p = NN_Predictor.make_prediction( draft_reaction_ids, trainedNN = path_to_NN) 
weights = {}
for i in p:
    weights[i]  = np.round(1-p[i], 10)



model_NN_gf, obj, new_reacs = gapfill_function.gapfill( all_reactions, 
                                                        draft_reaction_ids, 
                                                        weights, 
                                                        'bio1', 
                                                        result_selection = 'min_reactions'
                                                    )

ref_model = build_model.refine_model(model_NN_gf, draft_model)

for reaction in ref_model.reactions:
     if not draft_model.reactions.has_id(reaction.id):
            print(reaction.id, "~~", reaction.build_reaction_string(use_metabolite_names = 1))
            
cobra.io.write_sbml_model(ref_model, filename = "bh_gapfilled_model.sbml")