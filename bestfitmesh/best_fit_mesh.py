# -*- coding: utf-8 -*-
"""
Created on Thu May 17 22:06:18 2018

@author: Richard Hollamby, RIHY, COWI UK Bridge
"""

import numpy
import matplotlib.pyplot as plt


import element
import node


class BestFitMesh():
    """
    Class to analyse point cloud data for best fit mesh
    """
    
    def __init__(self,
                 length_x:float, length_y:float,
                 Le_x:float, Le_y:float,
                 verbose=True):
        """
        Initialisation function
        
        ***
        Required:
            
        * `points_fName`, string to denote .csv file defining point cloud data
        
        * `length_x`, `length_y`; _floats_ to define dimensions of rectangular 
          domain to which point cloud data applies
        
        * `Le_x`, `Le_y`; _floats_ to denote target element size in x- 
          and y- directions respectively
        
        """

        # Define elements & nodes
        self.DefineMesh(length_x,length_y,Le_x,Le_y,verbose=verbose)

        

        
    def DefineMesh(self,
                   L_x,L_y,
                   Le_x_target,Le_y_target,
                   verbose=False):
        """
        Define a grid of elements and nodes
        
        Note: grid wraps in x- direction i.e. first and last columns of 
        elements connect to 1st column of nodes
        """
        
        self.L_x = L_x
        self.L_y = L_y
        
        # Calculate number of elements in each grid direction
        nElements_x = int(L_x / Le_x_target + 1)
        nElements_y = int(L_y / Le_y_target + 1)
        self.nElements_x = nElements_x
        self.nElements_y = nElements_y
        
        # Calculate actual length dimensions for elements
        Le_x = L_x / nElements_x
        Le_y = L_y / nElements_y
        self.Le_x = Le_x
        self.Le_y = Le_y
        
        # Define node objects
        self.node_dict = self.define_nodes(Le_x,Le_y,nElements_x,nElements_y)
        """
        _Dict_ of `Node` objects, defining mesh nodes
        """
        
        # Define element objects
        self.element_dict = self.define_elements(nElements_x,nElements_y)
        """
        _Dict_ of `Element` objects, defining elements
        """
        
        if verbose:
            self.print_mesh_details()
        
            
    def define_nodes(self,Le_x:float,Le_y:float,
                    nElements_x:int,nElements_y:int):
        """
        Define positions of nodes, arranged in rectangular grid which wraps in 
        x- direction
        
        Nodes are numbered row-wise in ascending order
        """

        # Determine number of nodes
        # Note nodes wrap-around in circumferential x- direction
        nNodes_x = nElements_x
        nNodes_y = nElements_y + 1
        self.nNodes_x = nNodes_x
        self.nNodes_y = nNodes_y
        
        # Determine xy coordinates of nodes on regular grid
        x_vals = numpy.arange(0,Le_x*nNodes_x,Le_x)
        y_vals = numpy.arange(0,Le_y*nNodes_y,Le_y)
        
        node_xy = numpy.array([[x,y] for y in y_vals for x in x_vals])
    
        # Define node objects
        node_dict = {}
        
        for n in range(node_xy.shape[0]):
            
            node_dict[n] = node.Node(ID=n,x=node_xy[n,0],y=node_xy[n,1])
            
        return node_dict
    
    def define_elements(self,n_x,n_y):
        """
        Defines grid of elements
        
        Elements are numbered row-wise in ascending order
        
        Each element is attached to 4no corner nodes. Topologically these are 
        assigned clockwise, starting from bottom-left
        """
        
        element_dict = {}
        
        # Loop over rows and columns  
        index_rc = [(r,c) for r in range(n_y) for c in range(n_x)]
        
        for e, (r,c) in enumerate(index_rc):
                        
            # Define IDs of connecting nodes
            node0 = n_x*r + c
            node1 = node0 + 1
            
            if c==n_x-1: # correct for wrap-around effect
                node1 = node1 + - n_x
                
            node2 = node1 + n_x
            node3 = node0 + n_x
            
            node_IDs=[node0,node1,node2,node3]
            
            # Get node objects, given IDs
            node_objs = [self.node_dict[x] for x in node_IDs]
            
            # Define element
            element_dict[e] = element.Element(ID=e,connectedNodes=node_objs)
        
        return element_dict
    
    
    def print_mesh_details(self):
        """
        Prints details of mesh
        """

        print("Length in x:\t%.2fm" % self.L_x)
        print("Length in y:\t%.2fm" % self.L_y)
        print("nElements_x:\t%d" % self.nElements_x)
        print("nElements_y:\t%d" % self.nElements_y)
        print("Le_x:\t\t%.2f" % self.Le_x)
        print("Le_y:\t\t%.2f" % self.Le_y)
        print("nNodes_x:\t%d" % self.nNodes_x)
        print("nNodes_y:\t%d" % self.nNodes_y)
        
        print("nElements:\t%d" % len(self.element_dict))
        print("nNodes:\t\t%d" % len(self.node_dict))
        
    def plot_mesh(self):
        
        fig, ax = plt.subplots(1)
        fig.suptitle("Mesh numbering")
           
        for elem_obj in self.element_dict.values():
            
            IDs = [node.ID for node in elem_obj.connectedNodes]
            x_vals = [node.x for node in elem_obj.connectedNodes]
            y_vals = [node.y for node in elem_obj.connectedNodes]
                
            ax.plot(x_vals,y_vals,'k',marker='.',markeredgecolor='r')
            
            for _x,_y,_ID in zip(x_vals,y_vals,IDs):
                ax.text(_x,_y,"%d" % _ID,color='r')
                
            ax.text(numpy.mean(x_vals),numpy.mean(y_vals),"%d" % elem_obj.ID)
            
        ax.set_xlim([0,self.L_x])
        ax.set_ylim([0,self.L_y])
        ax.set_xlabel("Distance in x")
        ax.set_ylabel("Distance in y")
        
                
# Test routine
if __name__ == "__main__":
    
    analysis = BestFitMesh(1.0,2.0,0.4,0.52)
    analysis.plot_mesh()

#
#        # -----
#
#        # Define mesh
#        global nodeCoords, nNodes, nDofs
#        nodeCoords, nNodes, nDofs = DefineNodes()
#
#        global elementConnectivity
#        elementConnectivity = DefineElements()
#
#        # -----
#
#        # Read in results data from file
#        resultsData = npy.asmatrix(npy.genfromtxt("Data/" + fName,delimiter=",",skip_header=0))
#
#        Xi = resultsData[:,0]
#        Yi = resultsData[:,1]
#        Zi = resultsData[:,2]
#
#        # Check data 
#        Xi,Yi,Zi = CheckInputData(Xi,Yi,Zi)
#
#        # -----
#
#        # Calculate least squares mesh   
#        A_matrix, x_est, b_est, residuals, colCounts, freedomsToUse, elementIndexs= SolveForLeastSquaresMesh(Xi,Yi,Zi,minPivotTol)
#
#        # -----
#
#        # Re-form full freedoms vector, using NaN to denote freedoms not obtained by least-squares solution
#        global x_full
#        x_full = npy.asmatrix(numpy.zeros((nDofs,1)))
#        x_full[freedomsToUse] = x_est
#
#        # Populate boolean vector to denote which nodes have their freedoms defined
#        global nodeIsValid
#        nodeIsValid = numpy.zeros((nNodes),dtype=bool)
#        nodeIsValid[freedomsToUse//4] = True
#        
#        # Count number of points within each element
#        elementIndexs = npy.ravel(elementIndexs)
#        global elementCounts 
#        elementCounts = npy.zeros((nElements,))
#
#        for e in elementIndexs:
#            elementCounts[e] = elementCounts[e] + 1
#
#        # Determine whether element is valid based on number of survey data points
#        global elementIsValid
#        elementIsValid = npy.zeros_like(elementCounts,dtype="bool")
#        elementIsValid = npy.where(elementCounts>elementCountTolerance,True,False)
#            
#        # -----
#
#        # Create plots to show results
#
#        zLimits = [-6,6]
#        
#        print("\nPlotting best fit surface, residuals and calculation results")
#
#        fig = figure(num=None, figsize=(20, 16), dpi=80)
#
#        gs = gridspec.GridSpec(3, 3, height_ratios=[3,1,1],width_ratios=[1, 1,1]) 
#
#        fig.suptitle("%s" % fName,fontsize=16)
#
#        ax1 = fig.add_subplot(gs[0], projection='3d')
#        ScatterPlot(Xi,Yi,Zi,zLimits,10,"Raw data in s-z-Dr coordinates",ax1)
#
#        ax2 = fig.add_subplot(gs[1], projection='3d')
#        ScatterPlot(Xi,Yi,b_est,zLimits,10,"Best fit surface: Dr (mm)",ax2)
#
#        ax3 = fig.add_subplot(gs[2], projection='3d')
#        ScatterPlot(Xi,Yi,residuals,zLimits,10,"Residuals (mm)",ax3)
#
#        ax4 = fig.add_subplot(gs[3])
#        PlotResidualsHistogram(residuals,zLimits,"Histogram of residuals\n",ax4)
#
#        ax5 = fig.add_subplot(gs[4])
#        PlotCounts(colCounts[::4].reshape((nElements_y+1,nElements_x)),"Count of survey points relating to each node",ax5)
#
#        ax6 = fig.add_subplot(gs[5])
#        PlotIsValid(nodeIsValid.reshape((nElements_y+1,nElements_x)),"Under-defined nodes",ax6)
#        
#        ax7 = fig.add_subplot(gs[7])
#        PlotCounts(elementCounts.reshape((nElements_y,nElements_x)),"Count of survey points relating to each element",ax7)
#        
#        ax8 = fig.add_subplot(gs[8])
#        PlotIsValid(elementIsValid.reshape((nElements_y,nElements_x)),"Under-defined elements",ax8)
#
#        savefig("Plots/" + fName_root +"_1.png",transparent=0, bbox_inches='tight')
#        plt.close(fig)
#
#        # -----
#
#        # ----- POST-PROCESSING TO GIVE GAUGE IMPERFECTIONS -----
#
#        # Define grid of gauge positions
#        gauge_x = npy.linspace(0,length_x,200)
#        gauge_y = npy.linspace(0,length_y,200)
#        gaugeX, gaugeY = npy.meshgrid(gauge_x,gauge_y)
#        gauge_x = npy.asmatrix(npy.ravel(gaugeX)).T
#        gauge_y = npy.asmatrix(npy.ravel(gaugeY)).T
#        gauge_xy = npy.hstack((gauge_x,gauge_y))
#        nGaugePos = gauge_x.shape[0]
#
#        # Calculate gauge lengths
#        lg_x_short, lg_x_long, lg_y = CalcGaugeLengths(r,t,L)
#
#        # Calculate codified imperfection limits
#        # Refer Table 8.4 of BS EN 1993-1-6
#        U0max = [0.006,0.010,0.016]
#        limit_w_x_short = npy.multiply(lg_x_short,U0max) * 1000
#        limit_w_x_long = npy.multiply(lg_x_long,U0max) * 1000
#        limit_w_y = npy.multiply(lg_y,U0max) * 1000
#
#        # Calculate imperfections by using best fit mesh
#        w_x_short, gauge_IsValid_w_x_short = CalcGaugeImperfection(gauge_xy,npy.asmatrix([1,0]),lg_x_short)
#        w_x_short = w_x_short.reshape(gaugeX.shape)
#        gauge_IsValid_w_x_short = gauge_IsValid_w_x_short.reshape(gaugeX.shape)
#        w_x_short_count = npy.count_nonzero(gauge_IsValid_w_x_short) / nGaugePos
#
#        w_x_long, gauge_IsValid_w_x_long = CalcGaugeImperfection(gauge_xy,npy.asmatrix([1,0]),lg_x_long)
#        w_x_long = w_x_long.reshape(gaugeX.shape)
#        gauge_IsValid_w_x_long = gauge_IsValid_w_x_long.reshape(gaugeX.shape)
#        w_x_long_count = npy.count_nonzero(gauge_IsValid_w_x_long) / nGaugePos
#
#        w_y, gauge_IsValid_w_y = CalcGaugeImperfection(gauge_xy,npy.asmatrix([0,1]),lg_y)
#        w_y = w_y.reshape(gaugeX.shape)
#        gauge_IsValid_w_y = gauge_IsValid_w_y.reshape(gaugeX.shape)
#        w_y_count = npy.count_nonzero(gauge_IsValid_w_y) / nGaugePos
#
#        # Classify imperfection results
#        classification_w_x_short = ClassifyData(w_x_short,limit_w_x_short).reshape(gaugeX.shape)
#        classification_w_x_long = ClassifyData(w_x_long,limit_w_x_long).reshape(gaugeX.shape)
#        classification_w_y = ClassifyData(w_y,limit_w_y).reshape(gaugeX.shape)
#
#        # Plot results
#        fig = figure(num=None, figsize=(20, 16), dpi=80)
#
#        fig.suptitle("%s" % fName,fontsize=16)
#
#        print("Plotting gauge imperfection results")
#        
#        ax1 = fig.add_subplot(231)
#        titleStr = "Meridional gauge: $L_g$ = %.0fmm\n" % (1000*lg_y)
#        xlabelStr = "Surveyed area: %.0f %%" % (100*w_y_count)
#        PlotImperfectionResults(w_y,limit_w_y[2],titleStr,xlabelStr,ax1)
#        
#        ax2 = fig.add_subplot(232)
#        titleStr = "Short circumferential gauge: $L_g$ = %.0fmm\n" % (1000*lg_x_short)
#        xlabelStr = "Surveyed area: %.0f %%" % (100*w_x_short_count)
#        PlotImperfectionResults(w_x_short,limit_w_x_short[2],titleStr,xlabelStr,ax2)
#        
#        ax3 = fig.add_subplot(233)
#        titleStr = "Long circumferential gauge: $L_g$ = %.0fmm\n" % (1000*lg_x_long)
#        xlabelStr = "Surveyed area: %.0f %%" % (100*w_x_long_count)
#        PlotImperfectionResults(w_x_long,limit_w_x_long[2],titleStr,xlabelStr,ax3)
#
#        ax4 = fig.add_subplot(234)
#        CategoryPixelPlot(classification_w_y,ax4)
#        
#        ax5 = fig.add_subplot(235)
#        CategoryPixelPlot(classification_w_x_short,ax5)
#        
#        ax6 = fig.add_subplot(236)
#        CategoryPixelPlot(classification_w_x_long,ax6)
#        
#        savefig("Plots/" + fName_root +"_2.png",transparent=0, dpi=80, bbox_inches='tight')
#        plt.close(fig)
#        
#    else:
#        print("File %d skipped, as defined by input" % f)
#    
#    print("\n******\n")
#    
#print("SCRIPT TERMINATED SUCCESSFULLY!")

