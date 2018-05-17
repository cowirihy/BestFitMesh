# -*- coding: utf-8 -*-
"""
Defines element class

@author: rihy
"""



class Element():
    
    def __init__(self,ID:int,connectedNodes):
        
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
