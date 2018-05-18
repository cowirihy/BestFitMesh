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
        
        self.points_list = []
        """
        List of points associated with element
        """
        
    def add_point(self,point_obj,check_domain=True,verbose=False):
        """
        Associates 'Point' object with element
        """
        
        if verbose:
            print("Point %d assigned to element %d" % (point_obj.ID,self.ID))
        
        self.check_point_within_domain(point_obj)
        
        self.points_list.append(point_obj)
        point_obj.element_obj = self
        
        
    def check_point_within_domain(self,point_obj) -> bool:
        """
        Tests whether point lies within the element domain
        """
        
        min_x = self.connectedNodes[0].x
        min_y = self.connectedNodes[0].y
        max_x = min_x + self.Le_x
        max_y = min_y + self.Le_y
        
        x = point_obj.x
        y = point_obj.y
        
        if x < min_x or x > max_x or y < min_y or y > max_y:
            raise ValueError("Point outside element domain!\n" + 
                             "Element ID: %d\n" % self.ID + 
                             "Domain: {0}".format([[min_x,max_x],[min_y,max_y]])+
                             "Point ID: %d\n" % point_obj.ID +
                             "Point x: %.2f\n" % point_obj.x + 
                             "Point y: %.2f\n" % point_obj.y)
            
        return True
        
        
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
            
