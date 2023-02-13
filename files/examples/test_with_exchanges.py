

#Change working directory
import os
import sys
from pathlib import Path
path = Path.cwd()
print(path)
scripts_path = os.path.join(path.parents[1], 'scripts')
sys.path.append(scripts_path)

import NN_Predictor
import numpy as np
import pandas as pd
import json
import cobra
import build_model


from reaction_class import Reaction
import gapfill_function



file_path = os.path.join(path.parents[1], 'files')
NN_path = os.path.join(file_path, 'NN')

model_path = os.path.join(file_path, 'models')
NN_MS = NN_Predictor.load_NN(path= os.path.join(NN_path, 'NN_MS.h5'))


# Here goes our medium 
def_med = {}
def_med["EX_cpd00027_e0"] = {'lower_bound': -3., 'upper_bound': 0., 'metabolites': {"cpd00027_e0":-1.}}
def_med["EX_cpd00001_e0"] = {'lower_bound': -4., 'upper_bound': 0., 'metabolites': {"cpd00001_e0":-1.}}


# The user needs to provide his/her own model instead of our s_infantis_base_model.sbml file 
draft_model = cobra.io.read_sbml_model(os.path.join(model_path, 's_infantis_base_model.sbml')) 
draft_reaction = Reaction(model= os.path.join(model_path, 's_infantis_base_model.sbml') )

biochem = os.path.join(file_path,  'biochemistry', 'reactions.tsv')
db_reactions = Reaction(biochem_input=biochem)

exchange_reacs = Reaction(model = os.path.join(model_path, 'exchangeReactions.sbml'), fixed_bounds = def_med)
for react in exchange_reacs.reactions:
    if react not in def_med:
        exchange_reacs.reactions[react]["lower_bound"] = 0


db_reactions.reactions = db_reactions.add_dict(exchange_reacs.reactions, db_reactions.reactions)


#Combine all reactions into one object
all_reactions           = Reaction(fixed_bounds = def_med)
all_reactions.reactions = all_reactions.add_dict(draft_reaction.reactions, db_reactions.reactions) 
draft_reaction_ids = set(draft_reaction.reactions)

p = NN_Predictor.make_prediction(draft_reaction_ids, trainedNN = NN_path)


weights = {}
for i in p:
    weights[i]  = np.round(1-p[i])

# By default, medium is "complete"
NN_gf_model, obj, new_reacs = gapfill_function.gapfill(all_reactions, draft_reaction_ids, weights, 'bio1', medium = def_med, result_selection = 'min_reactions')



ref_model = build_model.refine_model(NN_gf_model, draft_model)


for reaction in ref_model.reactions:
     if not draft_model.reactions.has_id(reaction.id):
            print(reaction.id, "~~", reaction.build_reaction_string(use_metabolite_names = 1))



