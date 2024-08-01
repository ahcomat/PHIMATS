from numba import njit, prange, int64, void, float64
import numpy as np

import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
#Set fonts
from matplotlib import rc, rcParams, cm
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Arial']})
#Use tex for math
rcParams['mathtext.fontset'] = 'cm'
rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"

def CPlot2D(Cx, dpi=100, cmap='RdBu_r'):
                
    fig, ax = plt.subplots(1, 1, figsize=(1*56.25/25.4, 1.5*56.25/25.4), dpi=dpi, tight_layout=True)
    
    c = ax.contourf(Cx, 256, cmap=cmap, vmin=Cx.min(), vmax=Cx.max())
    cbar = fig.colorbar(c, ax=ax, ticks=np.linspace(Cx.min(), Cx.max(), 4), orientation="horizontal", pad=0.05)

    cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
    cbar.ax.tick_params(labelsize=10)
    ax.set_aspect('equal', 'box')
    ax.set_xticks([])
    ax.set_yticks([])

    pass
    
def PlotLine(y = 2*56.25/25.4, x = 1.5*56.25/25.4, dpi=100):
    fig, ax = plt.subplots(figsize=(y, x), dpi=dpi, tight_layout=True)
    return ax   

@njit()
def CalcGradx_2D(Cx, dx):
    dCx = np.zeros_like(Cx)
    x,y = Cx.shape
    for i in prange(5, y-5):
        for j in range(5, x-5):
            dCx[i,j] =  - 0.00079365*Cx[i, j-5]   \
                        + 0.00992063*Cx[i, j-4]  \
                        - 0.05952380*Cx[i, j-3]  \
                        + 0.23809523*Cx[i, j-2]  \
                        - 0.83333333*Cx[i, j-1]  \
                        + 0.83333333*Cx[i, j+1]  \
                        - 0.23809523*Cx[i, j+2]  \
                        + 0.05952380*Cx[i, j+3]  \
                        - 0.00992063*Cx[i, j+4]  \
                        + 0.00079365*Cx[i, j+5]
    return dCx/dx

@njit()
def CalcGrady_2D(Cx, dx):
    dCx = np.zeros_like(Cx)
    x,y = Cx.shape
    for i in prange(5, y-5):
        for j in range(5, x-5):
            dCx[i,j] =  - 0.00079365*Cx[i-5, j]   \
                        + 0.00992063*Cx[i-4, j]  \
                        - 0.05952380*Cx[i-3, j]  \
                        + 0.23809523*Cx[i-2, j]  \
                        - 0.83333333*Cx[i-1, j]  \
                        + 0.83333333*Cx[i+1, j]  \
                        - 0.23809523*Cx[i+2, j]  \
                        + 0.05952380*Cx[i+3, j]  \
                        - 0.00992063*Cx[i+4, j]  \
                        + 0.00079365*Cx[i+5, j]  
    return dCx/dx

@njit(nogil=True)
def CalcLaplacian_2D(Cx, dx):
    
    laplacian = np.zeros_like(Cx)
    x,y = Cx.shape
    
    for i in prange(5, y-5):
        for j in range(5, x-5):
            laplacian[i,j] =  0.00031746*Cx[i-5, j] \
                            - 0.00496031*Cx[i-4, j] \
                            + 0.03968253*Cx[i-3, j] \
                            - 0.23809523*Cx[i-2, j] \
                            + 1.66666666*Cx[i-1, j] \
                            + 0.00031746*Cx[i, j-5] \
                            - 0.00496031*Cx[i, j-4] \
                            + 0.03968253*Cx[i, j-3] \
                            - 0.23809523*Cx[i, j-2] \
                            + 1.66666666*Cx[i, j-1] \
                            - 5.85444444*Cx[i, j]   \
                            + 1.66666666*Cx[i, j+1] \
                            - 0.23809523*Cx[i, j+2] \
                            + 0.03968253*Cx[i, j+3] \
                            - 0.00496031*Cx[i, j+4] \
                            + 0.00031746*Cx[i, j+5] \
                            + 1.66666666*Cx[i+1, j] \
                            - 0.23809523*Cx[i+2, j] \
                            + 0.03968253*Cx[i+3, j] \
                            - 0.00496031*Cx[i+4, j] \
                            + 0.00031746*Cx[i+5, j]
                            
    return laplacian/(dx**2)

@njit()
def CalcGradx_3D(Cx, dx):
    
    dCx = np.zeros_like(Cx)
    x,y,z = Cx.shape
    
    for i in prange(5, y-5):
        for j in range(5, x-5):
            for k in range(5, z-5):
                
                dCx[i,j,k] = - 0.00079365*Cx[i, j-5, k]  \
                             + 0.00992063*Cx[i, j-4, k]  \
                             - 0.05952380*Cx[i, j-3, k]  \
                             + 0.23809523*Cx[i, j-2, k]  \
                             - 0.83333333*Cx[i, j-1, k]  \
                             + 0.83333333*Cx[i, j+1, k]  \
                             - 0.23809523*Cx[i, j+2, k]  \
                             + 0.05952380*Cx[i, j+3, k]  \
                             - 0.00992063*Cx[i, j+4, k]  \
                             + 0.00079365*Cx[i, j+5, k]

    return dCx/dx


@njit()
def CalcGrady_3D(Cx, dx):
    
    dCx = np.zeros_like(Cx)
    x,y,z = Cx.shape
    
    for i in prange(5, y-5):
        for j in range(5, x-5):
            for k in range(5, z-5):
                
                dCx[i,j,k] = - 0.00079365*Cx[i-5, j, k]  \
                             + 0.00992063*Cx[i-4, j, k]  \
                             - 0.05952380*Cx[i-3, j, k]  \
                             + 0.23809523*Cx[i-2, j, k]  \
                             - 0.83333333*Cx[i-1, j, k]  \
                             + 0.83333333*Cx[i+1, j, k]  \
                             - 0.23809523*Cx[i+2, j, k]  \
                             + 0.05952380*Cx[i+3, j, k]  \
                             - 0.00992063*Cx[i+4, j, k]  \
                             + 0.00079365*Cx[i+5, j, k]

    return dCx/dx

@njit()
def CalcGradz_3D(Cx, dx):
    
    dCx = np.zeros_like(Cx)
    x,y,z = Cx.shape
    
    for i in prange(5, y-5):
        for j in range(5, x-5):
            for k in range(5, z-5):
                
                dCx[i,j,k] = - 0.00079365*Cx[i, j, k-5]  \
                             + 0.00992063*Cx[i, j, k-4]  \
                             - 0.05952380*Cx[i, j, k-3]  \
                             + 0.23809523*Cx[i, j, k-2]  \
                             - 0.83333333*Cx[i, j, k-1]  \
                             + 0.83333333*Cx[i, j, k+1]  \
                             - 0.23809523*Cx[i, j, k+2]  \
                             + 0.05952380*Cx[i, j, k+3]  \
                             - 0.00992063*Cx[i, j, k+4]  \
                             + 0.00079365*Cx[i, j, k+5]

    return dCx/dx

@njit(nogil=True)
def CalcLaplacian_3D(Cx, dx):
    
    laplacian = np.zeros_like(Cx)
    x,y,z = Cx.shape
    
    for i in prange(5, y-5):
        for j in range(5, x-5):
            for k in range(5,z-5):
                
                laplacian[i,j,k] =    0.00031746*Cx[i-5, j, k] \
                                    - 0.00496031*Cx[i-4, j, k] \
                                    + 0.03968253*Cx[i-3, j, k] \
                                    - 0.23809523*Cx[i-2, j, k] \
                                    + 1.66666666*Cx[i-1, j, k] \
                                    + 0.00031746*Cx[i, j-5, k] \
                                    - 0.00496031*Cx[i, j-4, k] \
                                    + 0.03968253*Cx[i, j-3, k] \
                                    - 0.23809523*Cx[i, j-2, k] \
                                    + 1.66666666*Cx[i, j-1, k] \
                                    + 0.00031746*Cx[i, j, k-5] \
                                    - 0.00496031*Cx[i, j, k-4] \
                                    + 0.03968253*Cx[i, j, k-3] \
                                    - 0.23809523*Cx[i, j, k-2] \
                                    + 1.66666666*Cx[i, j, k-1] \
                                    - 8.78166665*Cx[i, j, k]   \
                                    + 1.66666666*Cx[i, j, k+1] \
                                    - 0.23809523*Cx[i, j, k+2] \
                                    + 0.03968253*Cx[i, j, k+3] \
                                    - 0.00496031*Cx[i, j, k+4] \
                                    + 0.00031746*Cx[i, j, k+5] \
                                    + 1.66666666*Cx[i, j+1, k] \
                                    - 0.23809523*Cx[i, j+2, k] \
                                    + 0.03968253*Cx[i, j+3, k] \
                                    - 0.00496031*Cx[i, j+4, k] \
                                    + 0.00031746*Cx[i, j+5, k] \
                                    + 1.66666666*Cx[i+1, j, k] \
                                    - 0.23809523*Cx[i+2, j, k] \
                                    + 0.03968253*Cx[i+3, j, k] \
                                    - 0.00496031*Cx[i+4, j, k] \
                                    + 0.00031746*Cx[i+5, j, k]
                            
    return laplacian/(dx**2)

def CalcFluxGB_2D(theta, gPhi, dx, zeta_gb, DL, T):
    
    R = 8.31446261815324
    
    thetaGradx = CalcGradx_2D(theta, dx)
    thetaGrady = CalcGrady_2D(theta, dx)
    
    gPhiGradx = CalcGradx_2D(gPhi, dx)
    gPhiGrady = CalcGrady_2D(gPhi, dx)
    
    JxTotal_2D = 0
    JyTotal_2D = 0
                
    # case int(1):
    #     if varDiffusivity:
    #         JxTotal_2D = -DField*(thetaGradx - zeta_gb/(R*T)*theta*gPhiGradx)
    #         JyTotal_2D = -DField*(thetaGrady - zeta_gb/(R*T)*theta*gPhiGrady)
    #     else:
                    
    JxTotal_2D = -DL*(thetaGradx - zeta_gb/(R*T)*theta*gPhiGradx)
    JyTotal_2D = -DL*(thetaGrady - zeta_gb/(R*T)*theta*gPhiGrady)
                    
    # case int(2):
        
    #     JxTotal_2D = -DField*(thetaGradx - zeta_ff/(R*T)*theta*gPhi_ffGradx_2D \
    #                                                     - zeta_fM/(R*T)*theta*gPhi_fMGradx_2D \
    #                                                     - zeta_MM/(R*T)*theta*gPhi_MMGradx_2D \
    #                                                     - zeta_M/(R*T)*theta*martensiteGradx_2D)
        
    #     JyTotal_2D = -DField*(thetaGrady - zeta_ff/(R*T)*theta*gPhi_ffGrady_2D \
    #                                                     - zeta_fM/(R*T)*theta*gPhi_fMGrady_2D \
    #                                                     - zeta_MM/(R*T)*theta*gPhi_MMGrady_2D \
    #                                                     - zeta_M/(R*T)*theta*martensiteGrady_2D)
        
    # case int(3):
        
    #     JxTotal_2D = - DField*(6/Vm*thetaGradx \
    #                 - zeta_HAGB/(R*T)*theta*gPhi_HAGBGradx_2D \
    #                 - zeta_LAGB/(R*T)*theta*gPhi_LAGBGradx_2D)

        
    #     JyTotal_2D = - DField*(6/Vm*thetaGrady \
    #                 - zeta_HAGB/(R*T)*theta*gPhi_HAGBGrady_2D \
    #                 - zeta_LAGB/(R*T)*theta*gPhi_LAGBGrady_2D)
    
    # case int(4):
        
    #     JxTotal_2D = -DField*(thetaGradx - zeta_gb/(R*T)*theta*gPhiGradx_2D \
    #                                           - Vh/(R*T)*theta*sigma_hGradx_2D \
    #                                           - zeta_dis/(R*T)*theta*epslGradx_2D)
        
    #     JyTotal_2D = -DField*(thetaGrady - zeta_gb/(R*T)*theta*gPhiGrady_2D \
    #                                           - Vh/(R*T)*theta*sigma_hGrady_2D \
    #                                           - zeta_dis/(R*T)*theta*epslGrady_2D)

    return [JxTotal_2D, JyTotal_2D]