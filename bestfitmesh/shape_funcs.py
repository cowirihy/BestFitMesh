# -*- coding: utf-8 -*-
"""
Created on Fri May 18 14:08:37 2018

@author: rihy
"""

import numpy

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable




class ShapeFuncs():
    """
    Class for defining shape functions for Bogner-Fox-Schmit element
    """
    
    def __init__(self,applyScaling:bool=True,verbose=False,makePlots=False):
        """
        Initialisation function
        
        Power series coefficients are determined numerically, to give required 
        shape function behaviour
        """
        print("Initialising Bogner-Fox-Schmit shape function formulation...")
        
        self.polycoeffs = self._calc_polycoeffs(applyScaling=applyScaling)
        
        if verbose:
            print("Polynomial coefficients matrix:\n{0}"
                  .format(self.polycoeffs))
            
        if makePlots:
            self.plot_polycoeffs()
            self.plot()
            
        
    def evaluate(self,eta,nu):
        """
        Evaluate shape functions at local coordinates (eta,nu)
        """
        
        # Flatten arrays
        eta = numpy.ravel(eta)
        nu = numpy.ravel(nu)
        
        # Define polynomial terms vector
        polyterms = self._calc_poly(eta,nu)
        
        # Evaluate shape functions at specified (eta,nu)
        Nvector = polyterms * self.polycoeffs
    
        return Nvector
    
    
    def plot(self,ax=None,cmap='bwr'):
        """
        Produce 3D plots to illustrate shape functions 
        as defined by class methods
        """
        
        # Get figure
        fig = plt.figure()
        fig.suptitle("Shape functions for nodal freedoms",fontsize=16)
        fig.set_size_inches((8,8))

        # Define grid of (eta,nu) coordinates to evaluate shape functions at
        eta = numpy.arange(-1, 1, 0.01)
        nu = numpy.arange(-1, 1, 0.01)
        Eta, Nu = numpy.meshgrid(eta, nu)
    
        eta = numpy.ravel(Eta)
        nu = numpy.ravel(Nu)
    
        # Evaluate shape functions at (eta,nu) coordinates
        Zvals = self.evaluate(eta,nu)
        
        # Produce plots for each shape function as 4x4 grid
        for i in range(16):
        
            ax = fig.add_subplot(4,4, i+1, projection='3d')
    
            Z = Zvals[:,i].reshape(Eta.shape)
            
            absmax = numpy.max(numpy.abs(Z))
            
            ax.plot_surface(Eta, Nu, Z,
                            cmap=cmap,
                            linewidth=0,
                            vmin=-absmax, vmax=+absmax)
        
            ax.set_zlim([-1,1])
        
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
    
        
        return fig, ax
        
    
    def plot_polycoeffs(self,cmap='bwr'):
        """
        Produces pixel plot of polynomial coefficients matrix
        """
        
        fig, ax = plt.subplots()
        
        # Get axis for colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        
        # Make plot
        h = ax.imshow(self.polycoeffs,interpolation="none",cmap=cmap)
        fig.colorbar(h,cax=cax,orientation='vertical')
        
        ax.set_title("Polynomial coefficients matrix\n"+
                     "(used when evaluating shape functions)")
        
        return fig, ax


    
    def _calc_poly(self,eta,nu):
        """
        Calculates 3rd order 2D polynomial terms given local coords (eta,nu)
        """
        
        # Check shapes
        if eta.shape != nu.shape:
            raise ValueError("`eta` and `nu` to have same shape!")
    
        # Define polynomial terms matrix
        
        polyterms = numpy.zeros((16,eta.shape[0]))
        
        polyterms[0,:]=1
        
        polyterms[1,:]=eta
        polyterms[2,:]=nu
        
        polyterms[3,:]=eta**2
        polyterms[4,:]=eta*nu
        polyterms[5,:]=nu**2
        
        polyterms[6,:]=eta**3
        polyterms[7,:]=eta**2 * nu
        polyterms[8,:]=eta * nu**2
        polyterms[9,:]=nu**3
        
        polyterms[10,:]=eta**3 * nu
        polyterms[11,:]=eta**2 * nu**2
        polyterms[12,:]=eta * nu**3
        
        polyterms[13,:]=eta**3 * nu**2
        polyterms[14,:]=eta**2 * nu**3
        
        polyterms[15,:]=eta**3 * nu**3
        
        # Transpose and convert to matrix format
        polyterms = numpy.asmatrix(polyterms.T)
                
        return polyterms

    
    def _calc_poly_deriv_eta(self,eta,nu):
        """
        Calculates d/d_eta of shape function polynomial for given (eta,nu)
        """
        
        # Define polynomial terms matrix
        
        polyterms = numpy.zeros((16,eta.shape[0]))
        
        polyterms[0,:]=0
        
        polyterms[1,:]=1
        polyterms[2,:]=0
        
        polyterms[3,:]=2*eta
        polyterms[4,:]=nu
        polyterms[5,:]=0
        
        polyterms[6,:]=3*eta**2
        polyterms[7,:]=2*eta * nu
        polyterms[8,:]=nu**2
        polyterms[9,:]=0
        
        polyterms[10,:]=3*eta**2 * nu
        polyterms[11,:]=2*eta * nu**2
        polyterms[12,:]=nu**3
        
        polyterms[13,:]=3*eta**2 * nu**2
        polyterms[14,:]=2*eta * nu**3
        
        polyterms[15,:]=3*eta**2 * nu**3
        
        # Transpose and convert to matrix format
        polyterms = numpy.asmatrix(polyterms.T)
        
        return polyterms
    
    
    def _calc_poly_deriv_nu(self,eta,nu):
        """
        Calculate d/d_nu of shape function polynomial for given (eta,nu)
        """
                
        # Define polynomial terms matrix
        
        polyterms = numpy.zeros((16,eta.shape[0]))
        
        polyterms[0,:]=0
        
        polyterms[1,:]=0
        polyterms[2,:]=1
        
        polyterms[3,:]=0
        polyterms[4,:]=eta
        polyterms[5,:]=2*nu
        
        polyterms[6,:]=0
        polyterms[7,:]=eta**2
        polyterms[8,:]=eta * 2*nu
        polyterms[9,:]=3*nu**2
        
        polyterms[10,:]=eta**3
        polyterms[11,:]=eta**2 * 2*nu
        polyterms[12,:]=eta * 3*nu**2
        
        polyterms[13,:]=eta**3 * 2*nu
        polyterms[14,:]=eta**2 * 3*nu**2
        
        polyterms[15,:]=eta**3 * 3*nu**2
        
        # Transpose and convert to matrix format
        polyterms = numpy.asmatrix(polyterms.T)
        
        return polyterms
    
    
    def _calc_poly_deriv_etanu(self,eta,nu):
        """
        Calculate mixed derivative d2/d_eta_nu of shape function polynomial 
        for given (eta,nu)
        """
               
        # Define polynomial terms matrix
        
        polyterms = numpy.zeros((16,eta.shape[0]))
        
        polyterms[0,:]=0
        
        polyterms[1,:]=0
        polyterms[2,:]=0
        
        polyterms[3,:]=0
        polyterms[4,:]=1
        polyterms[5,:]=0
        
        polyterms[6,:]=0
        polyterms[7,:]=2*eta
        polyterms[8,:]=2*nu
        polyterms[9,:]=0
        
        polyterms[10,:]=3*eta**2
        polyterms[11,:]=2*eta * 2*nu
        polyterms[12,:]=3*nu**2
        
        polyterms[13,:]=3*eta**2 * 2*nu
        polyterms[14,:]=2*eta * 3*nu**2
        
        polyterms[15,:]=3*eta**2 * 3*nu**2
        
        # Transpose and convert to matrix format
        polyterms = numpy.asmatrix(polyterms.T)
        
        return polyterms
        
    
    def _calc_polycoeffs(self,applyScaling:bool=True):
        """
        Calculate polynomial coefficients matrix
        (which is required to define unit shape functions)
        """
        
        # Define polynomial terms matrix
        polyterms_matrix = numpy.asmatrix(numpy.zeros((16,16)))
        
        # Define corner coordinates for nodes 0-4
        # These are numbered clockwise from bottom-left corner of element
        eta = numpy.array([-1,+1,+1,-1])
        nu = numpy.array([-1,-1,+1,+1])
        
        # Evaluate shape function polynomial terms at corner coordinates
        u = self._calc_poly(eta,nu)
        u_deriv_eta = self._calc_poly_deriv_eta(eta,nu)
        u_deriv_nu = self._calc_poly_deriv_nu(eta,nu)
        u_deriv_etanu = self._calc_poly_deriv_etanu(eta,nu)
        
        # Assemble polynomial terms matrix
        polyterms_matrix[0::4,:]=u
        polyterms_matrix[1::4,:]=u_deriv_eta
        polyterms_matrix[2::4,:]=u_deriv_nu
        polyterms_matrix[3::4,:]=u_deriv_etanu
        
        # Polynomial coefficients matrix is inverse of this
        # such that identity matrix is obtained by multiplying matrices
        polycoeffs_matrix = numpy.linalg.inv(polyterms_matrix)
        
        if applyScaling:
            
            # Define scaling matrix
            # (used to scale the effect of freedoms, to improve numerical 
            # conditioning of least squares problem)
            sf = [1,3.38,3.38,11.39]
            S = numpy.diagflat([sf,sf,sf,sf])
            polycoeffs_matrix = polycoeffs_matrix * S
        
        return polycoeffs_matrix

    
# Test routines
if __name__ == "__main__":
    
    testRoutine = 0
    
    if testRoutine == 0:
        
        sf = ShapeFuncs()
        eta = numpy.linspace(-1,1,5)
        nu = numpy.linspace(-1,1,5)
        sf.evaluate(eta,nu)
    
    elif testRoutine == 1:
    
        pass
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
