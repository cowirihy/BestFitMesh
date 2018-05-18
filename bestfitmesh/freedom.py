# -*- coding: utf-8 -*-
"""
Defines freedom class

@author: rihy
"""



class Freedom():
    
    def __init__(self,node_obj,index:int):
                
        self.node_obj = node_obj
        """
        Node object to which freedom relates
        """
        
        self.index = index
        """
        _Integer_, denotes freedom type:
            
        * 0 = z
        
        * 1 = dz/dz
        
        * 2 = dz/dy
        
        * 3 = d2z/dxdy
        
        """
        
        self.value = None
        """
        _Float_, used to store freedom value, obtained by best fit analysis
        """

    def print_details(self):
        
        for key, val in vars(self).items():
            
            print("{0}:".format(key))
            print(val)
            print("")
        
        print("")