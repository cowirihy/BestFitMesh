# -*- coding: utf-8 -*-
"""
Defines element class

@author: rihy
"""

import numpy

from scipy.stats import kurtosis
from scipy.stats import skew


import shape_funcs


class Element():
    """
    Class to define rectangular element of dimensions [Le_x,Le_y]
    """
        
    # Define shape functions function applicable to this class of element
    # Hold as class attribute (common to all instances)
    shapeFuncs = shape_funcs.ShapeFuncs().evaluate
    
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
    
    
    def calc_stats(self,stat_name:str,use_residuals:bool=False):
        """
        Calculate statistics based on z- coordinates of points associated with 
        this element
        
        ***
        Accepted `stat_name` arguments:
            
        * 'num', returns number of points
        
        * 'mean', returns mean z value
        
        * 'std', returns standard deviation of z values
        
        * 'skewness', return skewness (3rd moment) of z values
        
        * 'kurtosis', return excess kurtosis (4th moment) of z values

    	 Also `ID` can be used to return element ID.
        
        ***
        Options:
            
        * `use_residuals`, can be set to express statistics based on z 
          residuals rather than z values themselves.
        
        """
        
        z_vals = [p.z for p in self.points_list]
        
        if use_residuals:
            z_vals = z_vals
            print("Not yet implemented!!!")

        if stat_name == 'ID':
            return self.ID
        
        elif stat_name == 'num':
            return len(z_vals)
        
        elif stat_name == 'mean':
            return numpy.mean(z_vals)
        
        elif stat_name == 'std' or stat_name == 'std dev':
            return numpy.std(z_vals)
        
        elif stat_name == 'skewness' or stat_name == 'skew':
            return skew(z_vals)
        
        elif stat_name == 'kurtosis' or stat_name == 'kurt':
            return kurtosis(z_vals,fisher=True) # 0.0 expected for normal dist
        
        else:
            raise ValueError("Unexpected `stat_name` passed")
        
        
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
            
            
    def get_point_pos(self,use_local=True):
        """
        Returns xy positions of all points associated with element
        
        If option `use_local` is set, positions are expressed in local 
        (eta,nu) coordinate system
        """
        
        xi = [p.x for p in self.points_list]
        yi = [p.y for p in self.points_list]
        zi = [p.z for p in self.points_list]
        
        if use_local:
            
            node0 = self.connectedNodes[0]
            xmin = node0.x
            ymin = node0.y
            
            eta_i = 2*(xi - xmin) / self.Le_x - 1
            nu_i = 2*(yi - ymin) / self.Le_y - 1
            
            return eta_i, nu_i
        
        else:
            
            return xi, yi, zi
        
        
    def calc_shape_funcs(self):
        """
        Calculates shape function matrix for all points related to this element
        """
            
        # Get point positions in local coordinates
        eta, nu = self.get_point_pos(use_local=True)
        
        # Evaluate shape functions at these positions
        shape_func_vals = self.shapeFuncs(eta,nu)
        
        return shape_func_vals
    
    
    def get_freedoms(self):
        """
        Get list of `Freedom` objects related to this element
        
        Given each element has 4 nodes, and each node has 4 freedoms,
        a list of 16no freedom objects will be obtained
        """
        freedoms_list = []
        
        for node_obj in self.connectedNodes:
            freedoms_list += node_obj.freedoms_list
            
        return freedoms_list
    
    
    def assemble_regression_matrices(self,nFreedoms):
        """
        Assembles matrices `A` and `b` as per canonical form of linear 
        least squares regression.
        
        Shapes:
            
        * `A`: [nPoints,nFreedoms]
        
        * `b`: [nPoints,1]
        
        where _nPoints_ is the number of points related to this element and
        nFreedoms is the number of freedoms in the global problem.
        
        Both matrices are returned as numpy matrices
        """
        
        nPoints = len(self.points_list)
        
        # Evalute shape functions at point positions
        shapefunc_vals = self.calc_shape_funcs()
            
        # Get z coordinates of all points in this element
        # This is the b vector to be used in least squares problem
        b = self.get_point_pos(use_local=False)[2]
        b = numpy.asmatrix(b).T
        
#        for row, sf_vals in zip(range(nPoints),shapefunc_vals):
#            
#            print(row)
#            print(sf_vals.shape)
#            print(b_val)
#        
#        rows = list(range(nPoints))
#        
#        freedom_indexs = [f.ID for f in elem_obj.get_freedoms()]
#        
#        A_submtrx = sparse.coo_matrix((sf_vals,(rows, freedom_indexs)),
#                                      shape=(nPoints, nFreedoms))
#        print(A_submtrx.shape)
            
            
    
# Test routines    
if __name__ == "__main__":
    
    testRoutine = 0
    
    if testRoutine == 0:
        
        pass
    
    elif testRoutine == 1:
    
        pass
            
        