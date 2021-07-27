# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 11:17:28 2021

@author: u0139894
"""
import cobra
from MSEED_compounds import Compounds
from MSEED_reactions import Reactions


#Use script from ModelSEED biochemistry to parse all metabolite/reaction info
compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
compounds_helper.saveCompounds(compounds_dict)
compounds_aliases_dict = compounds_helper.loadMSAliases()
compounds_helper.saveAliases(compounds_aliases_dict)
# compounds_names_dict = compounds_helper.loadNames()
# compounds_helper.saveNames(compounds_names_dict)


reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
reactions_aliases_dict = reactions_helper.loadMSAliases()
reactions_helper.saveAliases(reactions_aliases_dict)


def refine_model(gpfilledModel, draftModel=None):
    
    mc = gpfilledModel.copy()
    
    for metab in gpfilledModel.metabolites:
        
        
        
        
        metabName = metab.id.split('_')[0]
        metabSufix = metab.id.split('_')[1]
        
        if metabName in compounds_dict:
        
            # 1) Add metabolite info
            
            if metabName in compounds_aliases_dict:
                metab.annotation = compounds_aliases_dict[metabName].copy()
            else:
                metab.annotation = {}
                
            metab.charge = compounds_dict[metabName]['charge']
            
            metab.compartment = metabSufix
           
            if compounds_dict[metabName]['formula'] is None:
                compounds_dict[metabName]['formula']= ''
            
            metab.elements = compounds_helper.parseFormula(compounds_dict[metabName]['formula'])
            
            metab.formula = compounds_helper.buildFormula(compounds_helper.parseFormula(compounds_dict[metabName]['formula']))
            
            metab.name =  compounds_dict[metabName]['name'] + '_' + metabSufix
        else:
            #not in the modelSEED database (are specific to GapSeq)
            metab.annotation = {}
            metab.charge = None
            metab.compartment = metabSufix
            #metab.elements = ''
            metab.formula = ''
            metab.name = metab.id
    
    for reac in gpfilledModel.reactions:
        reacName = reac.id.split('_')[0]
        #reacSufix = reac.id.split('_')[1]
        
        if reacName in reactions_dict:
            
            if reacName in reactions_aliases_dict:
                reac.annotation = reactions_aliases_dict[reacName]
                
            else:
                reac.annotation = {}
        
            reac.name = reactions_dict[reacName]['name']
        
        else:#not in modelSEED
            reac.annotation = {}
            reac.name = reac.id
        
  
    if draftModel is not None:
        for reaction in draftModel.reactions:
            if len(reaction.genes)>0:
                if mc.reactions.has_id(reaction.id):
                    mc.reactions.get_by_id(reaction.id).gene_reaction_rule = reaction.gene_reaction_rule
                    
    #change internal fluxes to +/- 1000            
    for reaction in mc.reactions:
        
        if not reaction.boundary:
            reaction.lower_bound *=1000.
            reaction.upper_bound *= 1000.
            
        
            
    return mc
        
        


