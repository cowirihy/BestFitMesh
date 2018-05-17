# -*- coding: utf-8 -*-
"""
Defines element class

@author: rihy
"""

import numpy


class Element():
    
    def __init__(self,ID:int,connectedNodes,Le_x,Le_y):
        
        self.ID = ID
        """
        Integer ID for element
        """
        
        self.connectedNodes = connectedNodes
        """
        List of connected node objects
        """
        
        # Connect nodes
        for nodeObj in connectedNodes:
            nodeObj.connect_element(self)
            
        self.isActive = True
        """
        _Boolean_, denotes whether node is active or not
        """
        
        self.Le_x = Le_x
        """
        Length in x direction
        """
        
        self.Le_y = Le_y
        """
        Length in y direction
        """
        
    def print_details(self):
        """
        Print details of element to console
        """
        print("Element ID:\t%d" % self.ID)
        print("Nodes IDs:\t{0}".format([n.ID for n in self.connectedNodes]))
        
            
    def plot(self,ax,plot_ID=False):
        """
        Visualise element by plotting element onto given set of axes
        """
        
        x_vals = [node.x for node in self.connectedNodes]
        y_vals = [node.y for node in self.connectedNodes]
           
        if self.isActive:
            ax.plot(x_vals,y_vals,'k')
        else:
            ax.plot(x_vals,y_vals,'r')
        
        if plot_ID:
            ax.text(numpy.mean(x_vals),numpy.mean(y_vals),"%d" % self.ID)
            
