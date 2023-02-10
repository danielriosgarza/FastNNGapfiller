#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 13:19:13 2019

@author: meine
"""
import numpy as np
import pandas as pd
import os
import sys
import tensorflow as tf
from tensorflow.keras.models import Sequential
from pathlib import Path
import cobra
path = Path.cwd()
sys.path.append(path)
NN_path = os.path.join(path.parent, 'files', 'NN')


#Function that loads the Neural network, maybe unneccesary; path is path to .h5 file
def load_NN(path=None):
    if path is None:
        print('Loading Default NN (ModelSEED)')
        path = os.path.join(NN_path, 'NN_MS.h5')
    else:
        print('Loading user provided network')
    return tf.keras.models.load_model(path, custom_objects={"custom_loss": 'binary_crossentropy'})

def load_ids(path=NN_path+'/rxn_ids_ModelSEED.npy'):
    ids = np.load(path, allow_pickle=True).astype('str')
    return ids
#Function that makes a prediction based on input_data using Neural Network (NN)
def make_prediction(input, NN=None, rxn_ids=None):
    if rxn_ids is None:
        rxn_ids = load_ids()
    if isinstance(input, cobra.core.model.Model):
        input = convert_reaction_list(input.reactions.list_attr('id'), rxn_ids)
    else:
        if (isinstance(input, pd.DataFrame)):
            input.reindex(rxn_ids)
            input = input.T
        else:
            if isinstance(input, dict):
                input = convert_reaction_list([i for i in input if input[i]==1], rxn_ids)
            else:
                if not np.isin(input, [0,1]).all():
                    print('Converting to binary array:')
                    try:
                        input = convert_reaction_list(input)
                    except:
                        raise Exception("Conversion failed")
                else:
                    input = np.asarray(input.T)

    if NN is None:
        NN = load_NN()
    if isinstance(NN, str):
            print('Loading network')
            NN = load_NN(NN_path)
    if not isinstance(NN, Sequential):
        raise Exception('Type: {} not supported'.format(type(NN)))


    single_input=False
    #test for single input (trips up NN)
    if np.ndim(input) == 1:
        single_input=True
        input = np.expand_dims(input,axis=0)
    if input.shape[-1] == NN.input_shape[-1]:
        prediction = np.zeros(input.shape)
        #   load network and make prediction
        t_result = NN.predict(input)
        prediction = np.asarray(t_result)

        if single_input:
            prediction = dict(zip(rxn_ids, np.squeeze(prediction)))
    else:
        raise Exception("data has wrong shape: ", input.shape, 'instead of ', NN.input_shape)
    if isinstance(input, pd.DataFrame):
        prediction = pd.DataFrame(index=rxn_ids, columns=input.columns, data=prediction.T)
    return prediction

#function that generates a binary input based on a list of reaction ids
def convert_reaction_list(reaction_set, NN_reaction_ids=None):
    if NN_reaction_ids is None:
        if(list(reaction_set)[0][:3] == 'rxn'):
            model_type = 'ModelSEED'
            NN_reaction_ids = np.load(NN_path+'/rxn_ids_ModelSEED.npy', allow_pickle=True).astype('str')
        else:
            model_type = 'BiGG'
            NN_reaction_ids = np.load(NN_path+'/rxn_ids_bigg.npy', allow_pickle=True).astype('str')
        print('Using {} ids'.format(model_type))
    else:
        print('Using user-provided ids')
    b_input = []
    if(list(reaction_set)[0][:3] == 'rxn'):
        reaction_list = [reaction[0:8] + "_c0" for reaction in reaction_set]
    else:
        reaction_list = list(reaction_set)
    for i in NN_reaction_ids:
        if i in reaction_list:
            b_input.append(1)
        else:
            b_input.append(0)
    cheat = set([i for i in set(reaction_set) if not 'EX' in i])
    print("#reactions not in NN_rxn: ", len(cheat.difference(NN_reaction_ids)))
    if(sum(b_input)==0):
        raise Exception("No reactions found")
    return np.array(b_input)
