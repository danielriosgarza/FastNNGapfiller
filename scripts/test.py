# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 18:40:11 2021

@author: u0139894
"""

#Avoid the error warnings about the sbml formating
import logging
logging.getLogger("cobra").setLevel(logging.ERROR)

from reaction_class import Reaction
import gapfill_function


import os
from pathlib import Path
path = Path.cwd()


# Example 1. Gapfilling a model

# Load the model

model_location = os.path.join(path.parent, 'files',  'models', 'all_bins-draft.xml')
draft_model = Reaction(model=model_location) 

#draft_model now contains names, bounds and stoichiometric information of all reactions in the model.
#These are all stored in a dictionary, called draft_model.reactions.

draft_reaction_ids = set(draft_model.reactions) #Set of reaction ids in draft model.


# Load the database with reactions used for gapfilling
biochem = os.path.join(path.parent, 'files',  'biochemistry', 'reactions.tsv')

db_reactions = Reaction(biochem_input=biochem)

#db_reactions now contains names, bounds and stoichiometric information of all reactions in the database.

#Combine all reactions into one object

all_reactions = Reaction()
all_reactions.reactions = all_reactions.add_dict(draft_model.reactions, db_reactions.reactions) 
#add_dict returns all items from both input dictionaries

# Gapfill

gapfilled_model_location = os.path.join(path.parent, 'files',  'models', 'all_bins-draft_gpfilled.xml')

model, obj, new_reacs = gapfill_function.gapfill(all_reactions, draft_reaction_ids, {}, 'bio1', latendresse_result_selection = 'min_reactions')

