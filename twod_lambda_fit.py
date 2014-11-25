# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 11:56:21 2014

@author: jholt
"""

import pylab as pl
import scipy 
import scipy.optimize 

from mpl_toolkits.mplot3d import Axes3D 
import numpy as np

##****************************************************************************** 
##******************************************************************************
##------------------------------------------------------------------------------ 
def __residual(params, f, x, y): 
    """ 
    Define fit function; 
    Return residual error. 
    """ 
    a0, a1, a2, a3, a4, a5 = params 
    return np.ravel(a0 + a1*x + a2*x**2 + a3*y + a4*x*y + a5*(x**2)*y - f)

def __residual_linear(params, f, x, y):
    a0, a1, a2 = params 
    return np.ravel(a0 + a1*x + a2*y - f)

def twodfit(dataX, dataY, dataZ, logger, lower_len_points=10., sigma_max=0.5):
#    dataX = np.array([47, 347, 561, 565, 174, 240, 339, 360, 454, 479, 517, 648, 
#                      713, 892, 946, 81, 264, 771, 1010])
#    dataY0 = np.array([37, 37, 37, 37, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 
#                       35, 35, 35, 35])

    #dataY = 1./np.array(dataY0)
#    dataZ = np.array([20412.305, 20499.357, 20562.953, 20564.143, 21012.189, 
#                      21033.09, 21062.689, 21067.729, 21096.316, 21107.104, 
#                      21115.652, 21155.926, 21176.387, 21232.354, 21249.48,21580.498, 
#                      21636.973, 21802.111, 21873.348])
##    

    testing = False

    dataZ_orig = dataZ

    dataXX, dataYY = scipy.meshgrid(dataX,dataY) 
    nx,ny=100,100
   # dataXX,dataYY = np.meshgrid(np.linspace(dataX.min()-5, dataX.max()+5, nx), np.linspace(dataY.min()-0.005, dataY.max()+0.005, ny))
    
    ## guess initial values for parameters 
    p0 = [137.9, 0., 1./36, 750000, 10, 0.] 
    bad_points=[]
    #print __residual(p0, dataZZ, dataXX, dataYY) 
    sigma=100.
    p1, pcov, infodict, errmsg, success = scipy.optimize.leastsq(__residual, x0=p0, args=(dataZ, dataX, dataY), full_output=1)
    k=0

    if testing:
        pl.figure(14,figsize=(15,8))
        pl.clf()
        ax1=pl.subplot(411)
        pl.title("2d fitting")
                #pl.plot(sciorder.dx_2dfit,sciorder.cont,'r',label='avg of central rows')
                #pl.xlim([sciorder.dx_2dfit[0],sciorder.dx_2dfit[-1]])
        ax2=pl.subplot(412)
        ax3=pl.subplot(413)
        ax4=pl.subplot(414)

        datax_forplot=dataX
        datay_forplot=dataY
        dataz_forplot=dataZ
        points = ['r.','g.','c.','k.','m.','b.','y.',
                  'rx','gx','cx','kx','mx','bx','yx',
                  'r*','g*','c*','k*','m*','b*','y*','r.','g.','c.','k.','m.','b.','y.',
                  'rx','gx','cx','kx','mx','bx','yx',
                  'r*','g*','c*','k*','m*','b*','y*']

        lines = ['r-.','g.-','c-.','k-.','m-.','b-.','y-.',
                  'r--','g--','c--','k--','m--','b--','y--',
                  'r-','g-','c-','k-','m-','b-','y-','r-.','g.-','c-.','k-.','m-.','b-.','y-.',
                  'r--','g--','c--','k--','m--','b--','y--',
                  'r-','g-','c-','k-','m-','b-','y-']

    while len(dataZ) > lower_len_points-1. and sigma > sigma_max:
        
        p1, pcov, infodict, errmsg, success = scipy.optimize.leastsq(__residual, x0=p0, args=(dataZ, dataX, dataY), full_output=1) 
        
        newoh=np.ravel(p1[0] + p1[1]*dataX + p1[2]*dataX**2 + p1[3]*dataY + p1[4]*dataX*dataY + p1[5]*(dataX**2)*dataY ) 

        if (len(dataZ) > len(p0)) and pcov is not None:
                s_sq = (__residual(p1,dataZ,dataX,dataY)**2).sum() / (len(dataY) - len(p0))
                
                if testing:

                    ax1.plot(newoh, dataZ, points[k], newoh, dataZ, lines[k], label='fit')
                    ax2.plot(__residual(p1,dataz_forplot, datax_forplot, datay_forplot),
                            points[k],__residual(p1, dataz_forplot, datax_forplot, datay_forplot), lines[k], label=str(k)+' fit')
                    #ax2.ylabel('Residuals')
                    ax3.plot(dataX, dataZ,points[k],dataX, dataZ, lines[k], label=str(k)+' fit')
                    ax4.plot(1./dataY, dataZ,points[k],1./dataY, dataZ, lines[k], label=str(k)+' fit')

                var = ((__residual(p1,dataZ,dataX,dataY)**2).sum()) / (len(dataZ)-1)

                sigma = np.sqrt(var)

                residual=__residual(p1,dataZ,dataX,dataY)
                
                reduced_chi_sq=(residual**2).sum() / (len(dataY) - len(p0))

                covar = pcov*reduced_chi_sq
                std_error = np.array([np.sqrt(covar[i,i]) for i in range(len(p0))])
                
        else:
                pcov = 1000000.
          
        error =[]
        for j in range(len(p1)):
            try:
                error.append(np.absolute(pcov[j][j])**0.5)
            except:
                error.append(0.00)

        # new arrays made for second pass
        if sigma > sigma_max and len(dataZ) > lower_len_points:
            residuals=np.absolute(newoh-dataZ)
            bad_points.append(residuals.argmax())
            logger.info('stddev='+str(sigma))
            #logger.info('removed matched oh line from order '+str(1./dataY[residuals.argmax()]))
            #logger.info('removed matched oh line from fit '+str(dataX[residuals.argmax()]))
            #logger.info('removed matched oh line from fit '+str(dataZ[residuals.argmax()]))
            #print 'index = ',residuals.argmax(),' residual = ',residuals[residuals.argmax()]
            if testing:
                dataz_forplot[residuals.argmax()] = dataz_forplot[residuals.argmin()]
                datax_forplot[residuals.argmax()] = datax_forplot[residuals.argmin()]
                datay_forplot[residuals.argmax()] = datay_forplot[residuals.argmin()]

            dataZ = np.delete(dataZ, residuals.argmax())
            dataX = np.delete(dataX, residuals.argmax())
            dataY = np.delete(dataY, residuals.argmax())

        elif sigma > sigma_max:
            logger.info('cannot remove more points! we must live with sigma='+str(sigma))
            break
        
        k+=1
     
    #logger.info('sigma='+str(sigma))
    logger.info('std_error for each parameter='+str(std_error))

    dataZZ = p1[0] + p1[1]*dataXX + p1[2]*dataXX**2 + p1[3]*dataYY + p1[4]*dataXX*dataYY + p1[5]*(dataXX**2)*dataYY 
    if testing:
        pl.legend()
        pl.figure(15)
        pl.clf()
        pl.plot(newoh, dataZ,'b.',newoh,dataZ,'b')

    ## plot data 
#    pylab.close('all') 
#    fig = pylab.figure() 
#    ax = Axes3D(fig) 
#    ax.plot_wireframe(dataXX, dataYY, dataZZ) 
#    ax.scatter(dataX,dataY,dataZ,color='red')
#    ax.set_xlabel('x') 
#    ax.set_ylabel('y') 
#    ax.set_zlabel('z') 
#    pylab.show()
    
    return p1