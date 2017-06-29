# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 15:43:36 2017

@author: user
"""


def average(beamsize): 
    beam = 0.0
        
    for i in beamsize:
        beam += i
            
    avg = beam/len(beamsize)
#    bmaj = avg[0]
#    bmin = avg[1]
             
    return avg
