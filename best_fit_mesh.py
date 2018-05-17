# -*- coding: utf-8 -*-
"""
Created on Thu May 17 17:41:41 2018

@author: rihy
"""

import numpy as npy
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import ListedColormap
% pylab

# Set up progress bar functionality
from ipywidgets import FloatProgress
from IPython.display import display

from scipy.stats import norm

# Import from other local modules
import node

#%%
# *** MESHING FUNCTIONS ***/

# Define node positions
def DefineNodes():

    nNodes = int(nElements_x*(nElements_y+1))   # note nodes wrap-around in circumferential x- direction
    print("nNodes=%d" % nNodes)

    nodeCoords = npy.asmatrix(npy.zeros((nNodes,2)))

    for n in arange(0,nNodes):
        r = n // nElements_x
        c = n % nElements_x
        nodeCoords[n,0] = c*Le_x
        nodeCoords[n,1] = r*Le_y

    nDofs = int(4*nNodes)
    print("nDofs=%d" % nDofs)

    return nodeCoords, nNodes, nDofs
    
def DefineElements():
    
    # Each element relates to 4 nodes, numbered 0 to 3 working clockwise from (1,1)
    # Define element connectivity matrix whereby columns correspond to nodes 0 to 8 and rows denote each element
    # Matrix index is index of node

    elementConnectivity = npy.asmatrix(npy.zeros((nElements,4)))   # create blank array

    nodesPerRow = nElements_x

    for e in arange(0,nElements):

        rowIndex = e // nElements_x
        colIndex = e % nElements_x

        if colIndex == nElements_x-1:   # special case of last column

            elementConnectivity[e,2] = rowIndex*nodesPerRow + colIndex
            elementConnectivity[e,1] = elementConnectivity[e,2] + nodesPerRow
            elementConnectivity[e,3] = rowIndex*nodesPerRow
            elementConnectivity[e,0] = rowIndex*nodesPerRow + nodesPerRow

        else:
            elementConnectivity[e,2] = rowIndex*nodesPerRow + colIndex
            elementConnectivity[e,1] = elementConnectivity[e,2] + nodesPerRow
            elementConnectivity[e,0] = elementConnectivity[e,1] + 1
            elementConnectivity[e,3] = elementConnectivity[e,2] + 1

    return elementConnectivity

#%%

# *** FUNCTIONS USED TO EVALUATE SHAPE FUNCTIONS ***

# Calculate polynomial terms vector for given (eta,nu) coords
def checkshape(eta,nu):
    
    # Expect eta,nu to be row vectors
    if eta.shape[0]!=1 or nu.shape[0]!=1:
        print("Error: expected row vector!")
    
    if eta.shape[1]!=nu.shape[1]:
        print("Error: expected row vectors of identical length")

def calc_poly(eta,nu):
    
    checkshape(eta,nu)
    
    # Convert to numpy array format (to allow element-wise operations)
    eta = npy.array(eta)
    nu = npy.array(nu)
    
    # Define polynomial terms matrix
    # Row vectors are for each (eta,nu) coordinate
    polyterms = npy.zeros((16,eta.shape[1]))
    
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
    
    # Transpose
    polyterms = npy.asmatrix(polyterms.T)
    
    return polyterms

# Calculate d/d_eta of shape function polynomial for given (eta,nu) coords
def calc_poly_deriv_eta(eta,nu):
    
    # Convert to numpy array format (to allow element-wise operations)
    eta = npy.array(eta)
    nu = npy.array(nu)
    
    # Define polynomial terms matrix
    # Row vectors are for each (eta,nu) coordinate
    polyterms = npy.zeros((16,eta.shape[1]))
    
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
    
    # Transpose
    polyterms = npy.asmatrix(polyterms.T)
    
    return polyterms

# Calculate d/d_nu of shape function polynomial for given (eta,nu) coords
def calc_poly_deriv_nu(eta,nu):
    
    # Convert to numpy array format (to allow element-wise operations)
    eta = npy.array(eta)
    nu = npy.array(nu)
    
    # Define polynomial terms matrix
    # Row vectors are for each (eta,nu) coordinate
    polyterms = npy.zeros((16,eta.shape[1]))
    
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
    
    # Transpose
    polyterms = npy.asmatrix(polyterms.T)
    
    return polyterms

# Calculate mixed derivative d2/d_eta_nu of shape function polynomial for given (eta,nu) coords
def calc_poly_deriv_etanu(eta,nu):
       
    # Convert to numpy array format (to allow element-wise operations)
    eta = npy.array(eta)
    nu = npy.array(nu)
    
    # Define polynomial terms matrix
    # Row vectors are for each (eta,nu) coordinate
    polyterms = npy.zeros((16,eta.shape[1]))
    
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
    
    # Transpose
    polyterms = npy.asmatrix(polyterms.T)
    
    return polyterms
    
# Calculate polynomial coefficients matrix (which is used to define unit shape functions)
def calc_polycoeffs():
    
    # Define polynomial terms matrix
    polyterms_matrix = npy.asmatrix(npy.zeros((16,16)))
    
    # Define corner coordinates
    eta = npy.asmatrix([1,-1,-1,1])
    nu = npy.asmatrix([1,1,-1,-1])
    
    # Evaluate shape function polynomial terms at corner coordinates
    u=calc_poly(eta,nu)
    u_deriv_eta=calc_poly_deriv_eta(eta,nu)
    u_deriv_nu=calc_poly_deriv_nu(eta,nu)
    u_deriv_etanu=calc_poly_deriv_etanu(eta,nu)
    
    # Assemble polynomial terms matrix
    polyterms_matrix[0::4,:]=u
    polyterms_matrix[1::4,:]=u_deriv_eta
    polyterms_matrix[2::4,:]=u_deriv_nu
    polyterms_matrix[3::4,:]=u_deriv_etanu
    
    # Define scaling matrix (used to scale the effect of freedoms, to ensure least squares problem is well-conditioned)
    sf = [1,3.38,3.38,11.39]
    S = numpy.diagflat([sf,sf,sf,sf])
            
    # Polynomial coefficients matrix is inverse of this, such that identify matrix is obtained by multiplying matrices
    polycoeffs_matrix = npy.linalg.inv(polyterms_matrix) * S
    
    return polycoeffs_matrix
   
# Calculate shape function matrix given (eta, nu) and coeffs already calculated
def shapefunc_eta_nu(eta,nu):
    
    # Define polynomial terms vector
    polyterms = calc_poly(eta,nu)
    
    # Evaluate shape functions at specified (eta,nu)
    Nvector = polyterms * polycoeffs
    
    return Nvector

polycoeffs = calc_polycoeffs()
#imshow(polycoeffs,interpolation="none",cmap="bwr")
#colorbar()

# Plot shape functions, calculated using the above functions
def PlotShapeFunctions():

    eta = np.arange(-1, 1, 0.01)
    nu = np.arange(-1, 1, 0.01)
    Eta, Nu = np.meshgrid(eta, nu)

    eta = npy.asmatrix(npy.ravel(Eta))  # flattened
    nu = npy.asmatrix(npy.ravel(Nu))    # flattened

    Zvals = shapefunc_eta_nu(eta,nu)
    
    fig = figure(num=None, figsize=(14, 10), dpi=80)
    fig.suptitle("Shape functions for nodal freedoms",fontsize=16)

    for i in arange(0,16):
    
        ax = fig.add_subplot(4,4, i+1, projection='3d')

        Z = Zvals[:,i].reshape(Eta.shape)
        
        ax.plot_surface(Eta, Nu, Z,cmap=cm.coolwarm,linewidth=0,vmin=npy.min(Z), vmax=npy.max(Z))
    
        ax.set_zlim([-1,1])
    
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

    plt.show()
    
#PlotShapeFunctions()

#%%

# *** LEAST SQUARES FUNCTIONS ***/

def Define_A_matrix(xi,yi):
    
    nPoints = xi.shape[1]
     
    # Determine which element data point is within, given (x,y) coords 
    e = GetElement(xi,yi)
    
    # Determine which nodal freedoms are relevant for this element
    rF = GetRelevantNodalFreedoms(e)
    
    # Calculate position of point within element, in (eta,nu) coords
    eta_i, nu_i = GetEtaNuCoords(e,xi,yi)    
    
    # Evaluate shape functions at (eta,nu) co-ordinates
    shapeFuncs = shapefunc_eta_nu(eta_i,nu_i)
    
    # Populate sparse A matrix using shape function ordinates and relevant freedoms already determined
    A_matrix = npy.asmatrix(npy.zeros((nPoints,nDofs)))
    
    colCounts = npy.zeros((nDofs,))  # used to count number of non-zero entries per column
        
    for p in arange(0,nPoints):
        A_matrix[p,rF[p,:]] = shapeFuncs[p,:]
        colCounts[rF[p,:]] = colCounts[rF[p,:]] + 1

    return A_matrix, colCounts, e
    
# Set up least squares solution
def SolveForLeastSquaresMesh(Xi,Yi,Zi,minPivotTol):
            
    nResults = Zi.shape[0]
    
    print("No data points: %d" % nResults)
    print("No mesh dofs: %d" % nDofs)
    if nResults < nDofs:
        print("********** Warning: least squares problem is under-defined (nResults<nDofs) *********")
    
    # Define A matrix for least squares problem
    A_matrix, colCounts, elementIndexs = Define_A_matrix(npy.asmatrix(Xi).T,npy.asmatrix(Yi).T)
            
    # Retain only the nodal freedoms which have have sufficient pivots
    freedomsToUse = npy.ravel(argwhere(colCounts>minPivotTol))
    A_matrix = A_matrix[:,freedomsToUse]
      
    print("\nA_matrix shape:")
    print(A_matrix.shape)
    print("Solving for best-fit mesh using least squares")
    
    # Define b vector
    b = Zi
    
    # Solve Ax=b using psuedo inverse
    A_pinv = npy.linalg.pinv(A_matrix)   # normal Ax=b solver
    x_est = A_pinv * b
    
    # Alternatively use SciPy's solver for sparse matrices
    #x_est = spla.lsqr(A_matrix, b,show=True)
    
    b_est = A_matrix * x_est
    residuals = b - b_est
    
    print("-> Finished!")
    
    return A_matrix, x_est, b_est, residuals, colCounts,freedomsToUse,elementIndexs

#%%

# *** MISC FUNCTIONS ***/

def GetElement(xi,yi):   # returns element index, given (x,y) coords
   
    # Determine which element data point is within, given (x,y) coords  
    colNo = (xi // Le_x).astype(int)
    rowNo = (yi // Le_y).astype(int)
    
    # Correct for nodes that are on the boundary of the rectangular domain  
    colNo = npy.where(colNo< nElements_x,colNo,nElements_x-1)
    rowNo = npy.where(rowNo< nElements_y,rowNo,nElements_y-1)
        
    # Get element index
    e = nElements_x*(rowNo)+colNo
    e = e.astype(int)
    
    return e

def GetConnectedNodes(e):
    
    # Determine which nodes relate to this element
    connectedNodes = elementConnectivity[e,:].astype(int)
    
    return connectedNodes

def GetEtaNuCoords(e,xi,yi):
    
    # Determine which nodes relate to this element
    connectedNodes = GetConnectedNodes(e)
    
    # Calculate position of point within element, in (eta,nu) coords
    node2indexs = connectedNodes[:,2].T
    x_min_e = nodeCoords[node2indexs,0]
    y_min_e = nodeCoords[node2indexs,1]   
    eta_i = 2*(xi - x_min_e) / Le_x - 1
    nu_i = 2*(yi - y_min_e) / Le_y - 1
    
    return eta_i, nu_i

def GetRelevantNodalFreedoms(e):
    
    # Determine which nodes relate to this element
    connectedNodes = GetConnectedNodes(e)
    
    # Determine relevant nodal freedoms
    rF0 = connectedNodes[:,0]*4 + arange(0,4)
    rF1 = connectedNodes[:,1]*4 + arange(0,4)
    rF2 = connectedNodes[:,2]*4 + arange(0,4)
    rF3 = connectedNodes[:,3]*4 + arange(0,4)
    
    rF = npy.hstack((rF0,rF1))
    rF = npy.hstack((rF,rF2))
    rF = npy.hstack((rF,rF3))
    
    rF = rF.astype(int)
    
    return rF

def TestForValidNodes(connectedNodes):
    
    isValid_node0 = nodeIsValid[connectedNodes[:,0]]
    isValid_node1 = nodeIsValid[connectedNodes[:,1]]
    isValid_node2 = nodeIsValid[connectedNodes[:,2]]
    isValid_node3 = nodeIsValid[connectedNodes[:,3]]
    
    isValid_element = np.logical_and(isValid_node0,isValid_node1)
    isValid_element = np.logical_and(isValid_element,isValid_node2)
    isValid_element = np.logical_and(isValid_element,isValid_node3)
    
    return isValid_element

def IsWithinValidElement(xi,yi):
    
    # Get elements based on (x,y) coords
    e = GetElement(xi,yi)
    
    # Get connected nodes
    connectedNodes = GetConnectedNodes(e)
    
    # Test whether all connected nodes are valid (i.e. freedoms included in the least squares solution)
    validElements = TestForValidNodes(connectedNodes).T
    
    # Reject if outside [0,L] domain
    validElements = npy.where(xi>length_x,False,validElements)
    validElements = npy.where(xi<0,False,validElements)
    validElements = npy.where(yi>length_y,False,validElements)
    validElements = npy.where(yi<0,False,validElements)
    
    # Reject if not within a valid element
    validElements = np.logical_and(validElements,elementIsValid[e])
    
    return validElements

def CheckInputData(Xi,Yi,Zi):
    
    # Check data lies within meshed domain
    minX = Xi.min()
    maxX = Xi.max()
    minY = Yi.min()
    maxY = Yi.max()
    print("Range of input (x,y) data:")
    print("X: [%.2f, %.2f]" % (minX,maxX))
    print("Y: [%.2f, %.2f]" % (minY,maxY))
    
    # Modify data as necessary to lie in [0,L] domain in both x and y directions
    isModified = False
    
    # Shift data within xy plane
    if minX < 0:
        Xi = Xi - minX
        isModified = True
    if minY < 0:
        Yi = Yi - minY
        isModified = True

    minX = Xi.min()
    maxX = Xi.max()
    minY = Yi.min()
    maxY = Yi.max()

    # Stretch data within xy plane
    tolerance = 0.001
    
    if maxX-minX > length_x:
        sf = (length_x-tolerance) / (maxX-minX)
        Xi = sf * Xi
        isModified = True
        
    if maxY-minY > length_y:
        sf = (length_y-tolerance) / (maxY-minY)
        Yi = sf * Yi
        isModified = True
        
    minX = Xi.min()
    maxX = Xi.max()
    minY = Yi.min()
    maxY = Yi.max()
        
    # Subtract mean from Zi
    Zi = Zi - mean(Zi)
    
    if isModified:
        print("Range of input (x,y) data after modification:")
        print("X: [%.2f, %.2f]" % (minX,maxX))
        print("Y: [%.2f, %.2f]" % (minY,maxY))
        
    return Xi,Yi,Zi

#%%

# *** PLOTTING FUNCTIONS ***/

def ScatterPlot(Xi,Yi,Zi,zRange,markerSize,titleStr,ax=None):
    
    if ax is None:
        ax = plt.gca()
    
    # Flatten data
    Xi = npy.ravel(Xi)
    Yi = npy.ravel(Yi)
    Zi = npy.ravel(Zi)
        
    p=ax.scatter(Xi,Yi,Zi,c=Zi,cmap="bwr",s=markerSize,edgecolors="none",vmin=zRange[0],vmax=zRange[1])
    
    fig = plt.gcf()
    cbar=fig.colorbar(p)
    cbar.ax.tick_params(labelsize=10) 
    
    ax.view_init(90,-90)
    
    ax.set_xlim([0,length_x])
    ax.set_ylim([0,length_y])
    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    
    ax.set_xlabel("Circumferential position (m)",fontsize=8)
    ax.set_ylabel("Vertical position (m)",fontsize=8)
    
    ax.set_title(titleStr,fontsize=10)
    
def WireframePlot(Xi,Yi,Zi,zRange,markerSize,titleStr,ax=None):
    
    if ax is None:
        ax = plt.gca()
           
    p=ax.plot_wireframe(Xi,Yi,Zi,linewidth=0.5, colors='k',alpha=0.5)
    
    ax.set_xlim([0,length_x])
    ax.set_ylim([0,length_y])
    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    
    ax.set_xlabel("Circumferential position (m)")
    ax.set_ylabel("Vertical position (m)")
    
    ax.set_title(titleStr)    

def PlotResidualsHistogram(residuals,xrange,title,ax=None):
    
    if ax is None:
        ax = plt.gca()

    bins = numpy.linspace(xrange[0], xrange[1], 20)
    ax.hist(residuals,bins,normed=True)
    
    ax.set_xlabel("Residual (mm)",fontsize=8)
    ax.set_title(title,fontsize=10)

    # Overlay best fit normal dist
    mean_residuals = npy.mean(residuals)
    std_residuals = npy.std(residuals)
    
    x = np.linspace(norm.ppf(0.005,loc=mean_residuals,scale=std_residuals),norm.ppf(0.995,loc=mean_residuals,scale=std_residuals), 100)
    ax.plot(x, norm.pdf(x,loc=mean_residuals,scale=std_residuals),'r--', lw=1, alpha=0.6, label='norm pdf')

    ax.set_xlim(xrange)
    ax.set_yticks([])
    ax.tick_params(axis='x', which='major', labelsize=8)
    
    ax.text(0.8, 0.8,("Std err = %.2fmm" % std_residuals), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    
def PlotIsValid(isValid_mtrx,title,ax=None):

    if ax is None:
        ax = plt.gca()
        
    ax.imshow(isValid_mtrx,interpolation="None",cmap="Blues",vmin=0,vmax=1,origin='lower')
    ax.set_title(title,fontsize=10)
    
    ax.set_xlabel("Circumferential direction",fontsize=8)
    ax.set_ylabel("Vertical direction",fontsize=8)
    
    ax.set_xticks([])
    ax.set_yticks([])
    
def PlotCounts(counts_mtrx,title,ax=None):  
    
    if ax is None:
        ax = plt.gca()
       
    p=ax.imshow(counts_mtrx,interpolation="None",cmap="Greens",vmin=0,origin='lower')
    
    fig = plt.gcf()
    cbar=fig.colorbar(p)
    cbar.ax.tick_params(labelsize=10) 
    
    ax.set_title(title,fontsize=10)
    ax.set_xlabel("Circumferential direction",fontsize=8)
    ax.set_ylabel("Vertical direction",fontsize=8)
    
    ax.set_xticks([])
    ax.set_yticks([])
    
def PlotMesh(title,ax=None):

    if ax is None:
        ax = plt.gca()
    
#%%
    
# *** IMPERFECTION GAUGE FUNCTIONS ***

# Calculate gauge lengths
def CalcGaugeLengths(r,t,L):

    lg_x_short = 4 * (r*t)**0.5

    lg_x_long = 2.3 * (L**2 * r * t)**0.25
    if lg_x_long > r:
        lg_x_long = r

    lg_y = lg_x_short

    print("\nlg_x_short = %.3f" % lg_x_short)
    print("lg_x_long = %.3f" % lg_x_long)
    print("lg_y = %.3f" % lg_y)
    
    return lg_x_short, lg_x_long, lg_y

# Evaluate gauge imperfections
def CalcGaugeImperfection(gaugePos_xy,gaugeDir,gaugeLength):
    
    nPos = gaugePos_xy.shape[0]
    
    # Expect gaugePos to be a matrix whose row vectors are x y coords
    if gaugePos_xy.shape[1]!=2:
        print("Error: unexpected shape to gaugePos_xy argument!")
        print(gaugePos_xy.shape)
    
    # Expect gaugeDir to be a row vector
    if gaugeDir.shape[0]!=1 and gaugeDir.shape[1]!=2:
        print("Error: unexpected shape to gaugeDir argument!")
        print(gaugeDir.shape)
        
    # Expect gaugeLength to be a scalar
    if not(numpy.isscalar(gaugeLength)):
        print("Error: unexpected gaugeLength to be scalar!")

    # Calculate coordinates of points on gauge
    pos0 = gaugePos_xy - gaugeLength/2 * gaugeDir
    pos1 = gaugePos_xy
    pos2 = gaugePos_xy + gaugeLength/2 * gaugeDir
    pos = npy.vstack((pos0,pos1))
    pos = npy.vstack((pos,pos2))
    
    # Correct for wrap-around effect
    pos[:,0] = pos[:,0] % length_x
    
    # Set up Ax=b problem for valid points and use A and x_est to solve in forward manner for b_est
    A, c, e = Define_A_matrix(pos[:,0].T,pos[:,1].T)
    pos_z = A * x_full
    
    # Unpack end z coordinates from stacked array
    pos0_z = pos_z[0:nPos]
    pos1_z = pos_z[nPos:2*nPos]
    pos2_z = pos_z[2*nPos:3*nPos]
    
    # Calculate imperfection at middle of gauge
    w_imp = pos1_z - 0.5*(pos0_z+pos2_z)
    
    # Test whether gauge points are valid
    pos0_IsValid = IsWithinValidElement(pos0[:,0].T,pos0[:,1].T)
    pos1_IsValid = IsWithinValidElement(pos1[:,0].T,pos1[:,1].T)
    pos2_IsValid = IsWithinValidElement(pos2[:,0].T,pos2[:,1].T)
    
    gauge_IsValid = np.logical_and(pos0_IsValid,pos1_IsValid)
    gauge_IsValid = np.logical_and(gauge_IsValid,pos2_IsValid)
    gauge_IsValid = gauge_IsValid.T
            
    # Filter imperfection results
    w_imp = npy.where(gauge_IsValid,w_imp,npy.nan)
        
    return w_imp, gauge_IsValid

def PlotImperfectionResults(w,w_tol,titleStr,xlabel,ax):
    
    if ax is None:
        ax = plt.gca()
    
    p=ax.imshow(w,interpolation="None",cmap="bwr",vmin=-w_tol,vmax=+w_tol,origin='lower')
    fig = plt.gcf()
    fig.colorbar(p)
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.set_title(titleStr)
    ax.set_xlabel(xlabel)
    
# Function to classify data according to pre-defined limits
def ClassifyData(matrix_w,w_limits):
    
    # Create similar array filled with NaN
    classification_w = npy.zeros_like(matrix_w)
    classification_w.fill(npy.nan)

    # Loop through array
    for i in arange(0,classification_w.shape[0]):
        for j in arange(0,classification_w.shape[1]):

            abs_val = npy.abs(matrix_w[i,j])

            if abs_val <= w_limits[0]:
                classification_w[i,j] = 0

            elif abs_val <= w_limits[1]:
                classification_w[i,j] = 1

            elif abs_val <= w_limits[2]:
                classification_w[i,j] = 2

            elif npy.isnan(abs_val):
                classification_w[i,j] = npy.nan

            else:
                classification_w[i,j] = 3
                
    return classification_w

def CategoryPixelPlot(classification_w,ax):
    
    if ax is None:
        ax = plt.gca()
    
    cMap = ListedColormap(['green', 'orange', 'red','purple'])

    p=ax.imshow(classification_w,interpolation=None,cmap=cMap,vmin=0, vmax=3,origin='lower')
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    fig=plt.gcf()
    cbar=fig.colorbar(p)
    cbar.ax.get_yaxis().set_ticks([])
    
    for j, lab in enumerate(['     $Class A$','     $Class B$','    $Class C$','   $>Class C$']):
        cbar.ax.text(.5, (2 * j + 1) / 8.0, lab, ha='left',va='center')

#%%
        
# *** MAIN ROUTINE ***/

# All units in m

# Definitions file
defs_fName = "Best fit mesh - Definitions file.csv"
defs_data = npy.genfromtxt(defs_fName,dtype=None,delimiter=",",names=True)
nFiles = defs_data.shape[0]

# Initialise progress bar
pBar = FloatProgress(min=0, max=nFiles)
pBar.description = "Processing point cloud data"
display(pBar)

for f in arange(0,nFiles):
    
    # Update progress bar
    pBar.value += 1
    
    # Read details from definitions file
    fName = defs_data["Filename"][f].decode('UTF-8')
    OD = defs_data["OD"][f]                         # outer diameter in m
    t = defs_data["t"][f]                           # wall thickness in m
    L = defs_data["L"][f]                           # segment length / spacing between ring flanges
    Le_x_target = defs_data["Le_x_target"][f]        # target length of elements in x direction
    Le_y_target = defs_data["Le_y_target"][f]        # target length of elements in y direction 
    minPivotTol = defs_data["minPivotTol"][f]          # minimum number of points relating to given mesh node for nodal freedoms to be included
    elementCountTolerance = defs_data["elementCountTol"][f]   # minimum number of points within a given element for it to be considered well-defined
    solutionReqd = bool(defs_data["solutionReqd"][f])
    
    if solutionReqd:

        # Derive other geometry properties, given inputs provided
        r = (OD-t)/2                                       # midwall radius

        fName_root = str(fName[:(len(fName)-4)])
        print(fName_root+"\n")

        # Define rectanglar domain to mesh

        global length_x
        length_x = npy.pi * OD

        global length_y
        length_y = L

        print("Length x : %.2fm" % length_x)
        print("Length y : %.2fm" % length_y)

        global nElements_x
        nElements_x = length_x // Le_x_target + 1
        nElements_x = int(nElements_x)

        global nElements_y
        nElements_y = length_y // Le_y_target + 1
        nElements_y = int(nElements_y)

        print("nElements_x=%d" % nElements_x)
        print("nElements_y=%d" % nElements_y)

        global nElements
        nElements = nElements_x * nElements_y
        print("nElements=%d" % nElements)

        global Le_x
        global Le_y
        Le_x = length_x / nElements_x
        Le_y = length_y / nElements_y
        print("Le_x: %.2f" % Le_x)
        print("Le_y: %.2f" % Le_y)

        # -----

        # Define mesh
        global nodeCoords, nNodes, nDofs
        nodeCoords, nNodes, nDofs = DefineNodes()

        global elementConnectivity
        elementConnectivity = DefineElements()

        # -----

        # Read in results data from file
        resultsData = npy.asmatrix(npy.genfromtxt("Data/" + fName,delimiter=",",skip_header=0))

        Xi = resultsData[:,0]
        Yi = resultsData[:,1]
        Zi = resultsData[:,2]

        # Check data 
        Xi,Yi,Zi = CheckInputData(Xi,Yi,Zi)

        # -----

        # Calculate least squares mesh   
        A_matrix, x_est, b_est, residuals, colCounts, freedomsToUse, elementIndexs= SolveForLeastSquaresMesh(Xi,Yi,Zi,minPivotTol)

        # -----

        # Re-form full freedoms vector, using NaN to denote freedoms not obtained by least-squares solution
        global x_full
        x_full = npy.asmatrix(numpy.zeros((nDofs,1)))
        x_full[freedomsToUse] = x_est

        # Populate boolean vector to denote which nodes have their freedoms defined
        global nodeIsValid
        nodeIsValid = numpy.zeros((nNodes),dtype=bool)
        nodeIsValid[freedomsToUse//4] = True
        
        # Count number of points within each element
        elementIndexs = npy.ravel(elementIndexs)
        global elementCounts 
        elementCounts = npy.zeros((nElements,))

        for e in elementIndexs:
            elementCounts[e] = elementCounts[e] + 1

        # Determine whether element is valid based on number of survey data points
        global elementIsValid
        elementIsValid = npy.zeros_like(elementCounts,dtype="bool")
        elementIsValid = npy.where(elementCounts>elementCountTolerance,True,False)
            
        # -----

        # Create plots to show results

        zLimits = [-6,6]
        
        print("\nPlotting best fit surface, residuals and calculation results")

        fig = figure(num=None, figsize=(20, 16), dpi=80)

        gs = gridspec.GridSpec(3, 3, height_ratios=[3,1,1],width_ratios=[1, 1,1]) 

        fig.suptitle("%s" % fName,fontsize=16)

        ax1 = fig.add_subplot(gs[0], projection='3d')
        ScatterPlot(Xi,Yi,Zi,zLimits,10,"Raw data in s-z-Dr coordinates",ax1)

        ax2 = fig.add_subplot(gs[1], projection='3d')
        ScatterPlot(Xi,Yi,b_est,zLimits,10,"Best fit surface: Dr (mm)",ax2)

        ax3 = fig.add_subplot(gs[2], projection='3d')
        ScatterPlot(Xi,Yi,residuals,zLimits,10,"Residuals (mm)",ax3)

        ax4 = fig.add_subplot(gs[3])
        PlotResidualsHistogram(residuals,zLimits,"Histogram of residuals\n",ax4)

        ax5 = fig.add_subplot(gs[4])
        PlotCounts(colCounts[::4].reshape((nElements_y+1,nElements_x)),"Count of survey points relating to each node",ax5)

        ax6 = fig.add_subplot(gs[5])
        PlotIsValid(nodeIsValid.reshape((nElements_y+1,nElements_x)),"Under-defined nodes",ax6)
        
        ax7 = fig.add_subplot(gs[7])
        PlotCounts(elementCounts.reshape((nElements_y,nElements_x)),"Count of survey points relating to each element",ax7)
        
        ax8 = fig.add_subplot(gs[8])
        PlotIsValid(elementIsValid.reshape((nElements_y,nElements_x)),"Under-defined elements",ax8)

        savefig("Plots/" + fName_root +"_1.png",transparent=0, bbox_inches='tight')
        plt.close(fig)

        # -----

        # ----- POST-PROCESSING TO GIVE GAUGE IMPERFECTIONS -----

        # Define grid of gauge positions
        gauge_x = npy.linspace(0,length_x,200)
        gauge_y = npy.linspace(0,length_y,200)
        gaugeX, gaugeY = npy.meshgrid(gauge_x,gauge_y)
        gauge_x = npy.asmatrix(npy.ravel(gaugeX)).T
        gauge_y = npy.asmatrix(npy.ravel(gaugeY)).T
        gauge_xy = npy.hstack((gauge_x,gauge_y))
        nGaugePos = gauge_x.shape[0]

        # Calculate gauge lengths
        lg_x_short, lg_x_long, lg_y = CalcGaugeLengths(r,t,L)

        # Calculate codified imperfection limits
        # Refer Table 8.4 of BS EN 1993-1-6
        U0max = [0.006,0.010,0.016]
        limit_w_x_short = npy.multiply(lg_x_short,U0max) * 1000
        limit_w_x_long = npy.multiply(lg_x_long,U0max) * 1000
        limit_w_y = npy.multiply(lg_y,U0max) * 1000

        # Calculate imperfections by using best fit mesh
        w_x_short, gauge_IsValid_w_x_short = CalcGaugeImperfection(gauge_xy,npy.asmatrix([1,0]),lg_x_short)
        w_x_short = w_x_short.reshape(gaugeX.shape)
        gauge_IsValid_w_x_short = gauge_IsValid_w_x_short.reshape(gaugeX.shape)
        w_x_short_count = npy.count_nonzero(gauge_IsValid_w_x_short) / nGaugePos

        w_x_long, gauge_IsValid_w_x_long = CalcGaugeImperfection(gauge_xy,npy.asmatrix([1,0]),lg_x_long)
        w_x_long = w_x_long.reshape(gaugeX.shape)
        gauge_IsValid_w_x_long = gauge_IsValid_w_x_long.reshape(gaugeX.shape)
        w_x_long_count = npy.count_nonzero(gauge_IsValid_w_x_long) / nGaugePos

        w_y, gauge_IsValid_w_y = CalcGaugeImperfection(gauge_xy,npy.asmatrix([0,1]),lg_y)
        w_y = w_y.reshape(gaugeX.shape)
        gauge_IsValid_w_y = gauge_IsValid_w_y.reshape(gaugeX.shape)
        w_y_count = npy.count_nonzero(gauge_IsValid_w_y) / nGaugePos

        # Classify imperfection results
        classification_w_x_short = ClassifyData(w_x_short,limit_w_x_short).reshape(gaugeX.shape)
        classification_w_x_long = ClassifyData(w_x_long,limit_w_x_long).reshape(gaugeX.shape)
        classification_w_y = ClassifyData(w_y,limit_w_y).reshape(gaugeX.shape)

        # Plot results
        fig = figure(num=None, figsize=(20, 16), dpi=80)

        fig.suptitle("%s" % fName,fontsize=16)

        print("Plotting gauge imperfection results")
        
        ax1 = fig.add_subplot(231)
        titleStr = "Meridional gauge: $L_g$ = %.0fmm\n" % (1000*lg_y)
        xlabelStr = "Surveyed area: %.0f %%" % (100*w_y_count)
        PlotImperfectionResults(w_y,limit_w_y[2],titleStr,xlabelStr,ax1)
        
        ax2 = fig.add_subplot(232)
        titleStr = "Short circumferential gauge: $L_g$ = %.0fmm\n" % (1000*lg_x_short)
        xlabelStr = "Surveyed area: %.0f %%" % (100*w_x_short_count)
        PlotImperfectionResults(w_x_short,limit_w_x_short[2],titleStr,xlabelStr,ax2)
        
        ax3 = fig.add_subplot(233)
        titleStr = "Long circumferential gauge: $L_g$ = %.0fmm\n" % (1000*lg_x_long)
        xlabelStr = "Surveyed area: %.0f %%" % (100*w_x_long_count)
        PlotImperfectionResults(w_x_long,limit_w_x_long[2],titleStr,xlabelStr,ax3)

        ax4 = fig.add_subplot(234)
        CategoryPixelPlot(classification_w_y,ax4)
        
        ax5 = fig.add_subplot(235)
        CategoryPixelPlot(classification_w_x_short,ax5)
        
        ax6 = fig.add_subplot(236)
        CategoryPixelPlot(classification_w_x_long,ax6)
        
        savefig("Plots/" + fName_root +"_2.png",transparent=0, dpi=80, bbox_inches='tight')
        plt.close(fig)
        
    else:
        print("File %d skipped, as defined by input" % f)
    
    print("\n******\n")
    
print("SCRIPT TERMINATED SUCCESSFULLY!")

