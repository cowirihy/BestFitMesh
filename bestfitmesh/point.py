# -*- coding: utf-8 -*-
"""
Defines point class

@author: whoever
"""

class Point():
    
    def __init__(self,ID:int,xyz_list):
        
        self.ID = ID
        
        self.x = xyz_list[0]
        self.y = xyz_list[1]
        self.z = xyz_list[2]
        
        self.element_obj = None
        """
        Element object, within the xy domain of which point lies
        """
        
        
        
    def print_details(self):
        
        for key, val in vars(self).items():
            
            print("{0}:".format(key))
            print(val)
        
        print("")

