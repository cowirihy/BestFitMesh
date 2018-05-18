# -*- coding: utf-8 -*-
"""
Defines node class

@author: rihy
"""

import numpy as npy

import freedom

class Node():
    
    def __init__(self,ID,x,y):
        
        self.ID = ID
        """
        Integer ID for node
        """
        
        self.x = x
        """
        x position of node
        """
        
        self.y = y
        """
        y position of node"
        """
        
        self.elements = []
        """
        List of connected element objects
        """
        
        self.isActive = True
        """
        _Boolean_, denotes whether node is active or not
        """
        
        self.freedoms_list = [freedom.Freedom(4*ID+i,self,i) for i in range(4)]
        """
        List of freedom objects relating to this node
        """
        
        
    def connect_element(self,elementObj):
        
        self.elements.append(elementObj)
        
        
    def print_details(self):
        
        for key, val in vars(self).items():
            
            print("{0}:".format(key))
            print(val)
            print("")
        
        print("")
        
    def plot(self,ax,plot_ID:bool=False):
        """
        Visualise node by plotting onto given set of axes
        """
        
        if self.isActive:
            ms = '.'
            c = 'b'
        else:
            ms = 'x'
            c = 'r'
            
        ax.plot(self.x,self.y,color=c,marker=ms)
        
        if plot_ID:
            ax.text(self.x,self.y,"%d" % self.ID,color=c)
        

