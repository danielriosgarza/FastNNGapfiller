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
from pathlib import Path
path = Path.cwd()
sys.path.append(path)


#Function that loads the Neural network, maybe unneccesary; path is path to .h5 file
def load_NN(path):
    return tf.keras.models.load_model(path, custom_objects={"custom_loss": 'binary_crossentropy'})

#Function that makes a prediction based on input_data using Neural Network (NN)
def make_prediction(input_data, NN):

    single_input=False
    #test for single input (trips up NN)
    #workaround by just making the same prediction twice
    if np.ndim(input_data) == 1:
        single_input=True
        input_data = np.tile(input_data,[2,1])

    prediction = np.zeros(input_data.shape)

    #   load network and make prediction
    network = NN
    t_result = network.predict(input_data)
    prediction = np.asarray(t_result)

    #ugly i know. but do not have a different solution at this point
    if single_input:
        prediction  = prediction[0].T
    else:
        prediction = prediction.T

    return prediction

#function that generates a binary input based on a list of reaction ids
def convert_reaction_list(reaction_set):
    #load reaction ids
    file_path = os.path.join(path.parent, 'files', 'NN')
    reaction_ids = np.load(file_path+'/rxn_ids.npy').astype('str')
    b_input = []
    for i in reaction_ids:
        if i in reaction_set:
            b_input.append(1)
        else:
            b_input.append(0)

    return np.array(b_input)
