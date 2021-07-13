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


# Example 1. Gapfilling a model

# Load the model

model_location = 'C:/Users/u0139894/Documents/Bram/files/models/all_bins-draft.xml'
draft_model = Reaction(model=model_location) 

#draft_model now contains names, bounds and stoichiometric information of all reactions in the model.
#These are all stored in a dictionary, called draft_model.reactions.

draft_reaction_ids = set(draft_model.reactions) #Set of reaction ids in draft model.


# Load the database with reactions used for gapfilling
biochem = 'C:/Users/u0139894/Documents/Bram/files/biochemistry/reactions.tsv'

db_reactions = Reaction(biochem_input=biochem)

#db_reactions now contains names, bounds and stoichiometric information of all reactions in the database.

#Combine all reactions into one object

all_reactions = Reaction()
all_reactions.reactions = all_reactions.add_dict(draft_model.reactions, db_reactions.reactions) 
#add_dict returns all items from both input dictionaries

# Gapfill

gapfilled_model_location = 'C:/Users/u0139894/Documents/Bram/files/models/gapfilled_model.xml'

model, obj, new_reacs = gapfill_function.gapfill(all_reactions, draft_reaction_ids, {}, 'bio1', latendresse_result_selection = 'min_reactions')

# #Certain reactions can be given a custom cost for gap filling, by adding them to the now empty dictionary {}. See example 2.
# #The biomass name must be specified, here it is 'bio1'.
# #write_sbml is off on default, but here it was given a location to write the gapfilled model to.

# #result[0] is the delta of the gapfill solution from the latendresse algorithm.
# #result[1] tells you if medium defining reactions are checked. 'checked_EX' means they are checked, and unusable ones are deleted.
# #result[3] lists all the reaction ids that are in the gap filled model.