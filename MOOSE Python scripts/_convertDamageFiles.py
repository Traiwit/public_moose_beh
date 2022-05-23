#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 14:08:12 2022

@author: moose
"""


## you need damage_list.csv to run this script
## need to be in the folder that contain damage files

# this script removes the first column (node ID) from the damage file 

import pandas as pd
import csv
import numpy as np


damage_file_list = list(csv.reader(open('OLIMPIADA_damage_list.csv')))  #damage.csv

frame_num = len(damage_file_list)

for block in range(2, frame_num):

            damage_file = (damage_file_list[block][0])

            df = pd.read_csv(damage_file,header=None)
            first_column = df.columns[0]  
            df = df.drop([first_column], axis=1) #removing the element ID
                
            
            np.savetxt(damage_file, df ,fmt='%s',delimiter=",")  #save file
