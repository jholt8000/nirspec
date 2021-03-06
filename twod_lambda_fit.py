# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 11:56:21 2014

@author: jholt
"""

import pylab as pl
import scipy
import scipy.optimize

import numpy as np

import statsmodels.api as smapi
from statsmodels.formula.api import ols
import logging
import nirspec_util

# #******************************************************************************
# #******************************************************************************
# #------------------------------------------------------------------------------
def __residual(params, f, x, y):
    """ 
    Define fit function; 
    Return residual error. 
    """
    a0, a1, a2, a3, a4, a5 = params
    return np.ravel(a0 + a1 * x + a2 * x ** 2 + a3 * y + a4 * x * y + a5 * (x ** 2) * y - f)


def twodfit(dataX, dataY, dataZ, logger='', lower_len_points=10., sigma_max=0.5):
    """

    :param dataX: First independent variable
    :param dataY: Second independent variable
    :param dataZ: Dependent variable
    :param logger: logger instance
    :param lower_len_points: lowest
    :param sigma_max:
    :return:
    """
    if len(dataX) < 6:
        return [], [], []

    testing = True

    if not isinstance(logger, logging.Logger):
        logger = nirspec_util.NirspecBookkeeping.setup_logger('two_d_fit.log', '.', verbose=False)

    newoh = 9999

    dataXX, dataYY = scipy.meshgrid(dataX, dataY)

    # # guess initial values for parameters
    p0 = [137.9, 0., 1. / 36, 750000, 10, 0.]
    bad_points = []
    # print __residual(p0, dataZZ, dataXX, dataYY)
    sigma = 100.

    # This call is just to set up the plots
    p1, pcov, infodict, errmsg, success = scipy.optimize.leastsq(__residual, x0=p0, args=(dataZ, dataX, dataY),
                                                                 full_output=1)


    k = 0

    if testing:
        pl.figure(14, figsize=(15, 8))
        pl.clf()
        ax1 = pl.subplot(211)
        pl.title("2d fitting")
        ax2 = pl.subplot(212)

        points = ['r.', 'g.', 'c.', 'k.', 'm.', 'b.', 'y.',
                  'rx', 'gx', 'cx', 'kx', 'mx', 'bx', 'yx',
                  'r*', 'g*', 'c*', 'k*', 'm*', 'b*', 'y*', 'r.', 'g.', 'c.', 'k.', 'm.', 'b.', 'y.',
                  'rx', 'gx', 'cx', 'kx', 'mx', 'bx', 'yx',
                  'r*', 'g*', 'c*', 'k*', 'm*', 'b*', 'y*']

        lines = ['r-.', 'g.-', 'c-.', 'k-.', 'm-.', 'b-.', 'y-.',
                 'r--', 'g--', 'c--', 'k--', 'm--', 'b--', 'y--',
                 'r-', 'g-', 'c-', 'k-', 'm-', 'b-', 'y-', 'r-.', 'g.-', 'c-.', 'k-.', 'm-.', 'b-.', 'y-.',
                 'r--', 'g--', 'c--', 'k--', 'm--', 'b--', 'y--',
                 'r-', 'g-', 'c-', 'k-', 'm-', 'b-', 'y-']

        ax2.plot(__residual(p1, dataZ, dataX, dataY),
                         points[k], __residual(p1, dataZ, dataX, dataY), lines[k],
                         label=str(k) + ' fit')

    dataZ_new=np.copy(dataZ)
    dataY_new=np.copy(dataY)
    dataX_new=np.copy(dataX)

    residual = __residual(p1, dataZ, dataX, dataY)
    x_res = np.arange(len(residual))
    regression = ols("data ~ x_res", data=dict(data=residual, x=x_res)).fit()
    test = regression.outlier_test()
    outliers = ((x_res[i], residual[i]) for i,t in enumerate(test.icol(2)) if t < 0.9)
    #print 'outliers=',list(outliers)
    x=list(outliers)
    logger.info('residual outliers='+str(x))
    xhap=0

    for j in range(len(x)):
        logger.info('deleting outlier = '+str(dataX_new[x[j][0]-xhap])+' order = ' +str(1/dataY_new[x[j][0]-xhap])+' wavelength = '+str(dataZ_new[x[j][0]-xhap]))
        dataZ_new = np.delete(dataZ_new, x[j][0]-xhap)
        dataX_new = np.delete(dataX_new, x[j][0]-xhap)
        dataY_new = np.delete(dataY_new, x[j][0]-xhap)

        xhap+=1


    happened=0
    while len(dataZ_new) > lower_len_points - 1. and sigma > sigma_max:

        p1, pcov, infodict, errmsg, success = scipy.optimize.leastsq(__residual, x0=p0, args=(dataZ_new, dataX_new, dataY_new),
                                                                     full_output=1)


        residual = __residual(p1, dataZ_new, dataX_new, dataY_new)
        x_res = np.arange(len(residual))
        regression = ols("data ~ x_res", data=dict(data=residual, x=x_res)).fit()
        test = regression.outlier_test()
        outliers = ((x_res[i], residual[i]) for i,t in enumerate(test.icol(2)) if t < 0.9)
        #print 'outliers=',list(outliers)
        x=list(outliers)
        xhap=0

        for j in range(len(x)):
            logger.info('deleting outlier = '+str(dataX_new[x[j][0]-xhap])+' order = ' +str(1/dataY_new[x[j][0]-xhap])+' wavelength = '+str(dataZ_new[x[j][0]-xhap]))
            dataZ_new = np.delete(dataZ_new, x[j][0]-xhap)
            dataX_new = np.delete(dataX_new, x[j][0]-xhap)
            dataY_new = np.delete(dataY_new, x[j][0]-xhap)

            xhap+=1

        dataX_new_forplot = np.copy(dataX_new)
        dataY_new_forplot = np.copy(dataY_new)
        dataZ_new_forplot = np.copy(dataZ_new)

        newoh = np.ravel(p1[0] + p1[1] * dataX_new + p1[2] * dataX_new ** 2 + p1[3] * dataY_new + p1[4] * dataX_new * dataY_new + p1[5] * (
            dataX_new ** 2) * dataY_new)

        if (len(dataZ_new) > len(p0)) and pcov is not None:

            if testing:
                ax1.plot(newoh, dataZ_new, points[k], newoh, dataZ_new, lines[k], label='fit')
                ax2.plot(__residual(p1, dataZ_new_forplot, dataX_new_forplot, dataY_new_forplot),
                         points[k], __residual(p1, dataZ_new_forplot, dataX_new_forplot, dataY_new_forplot), lines[k],
                         label=str(k) + ' fit')


            residual = np.abs(__residual(p1, dataZ_new, dataX_new, dataY_new))
            var = ((residual ** 2).sum()) / (len(dataZ_new) - 1)

            sigma = np.sqrt(var)


        # new arrays made for second pass
        if sigma > sigma_max and len(dataZ_new) > lower_len_points:

            bad_points.append(residual.argmax())
            logger.info('sigma=' + str(sigma))

            if testing:
                dataZ_new_forplot[residual.argmax()-happened ] = dataZ_new_forplot[residual.argmin()]
                dataX_new_forplot[residual.argmax()-happened ] = dataX_new_forplot[residual.argmin()]
                dataY_new_forplot[residual.argmax()-happened ] = dataY_new_forplot[residual.argmin()]


            #logger.info('removing residual val='+str(residual[residual.argmax()])+' index = '+str(residual.argmax()))
            #logger.info(' removing datax='+str(dataX_new[residual.argmax()])+' datay='+str(1/dataY_new[residual.argmax()])+' dataz=',dataZ_new[residual.argmax()])

            dataZ_new = np.delete(dataZ_new, residual.argmax())
            dataX_new = np.delete(dataX_new, residual.argmax())
            dataY_new = np.delete(dataY_new, residual.argmax())
            happened +=1


        elif sigma > sigma_max:
            logger.info('cannot remove more points! we must live with sigma=' + str(sigma))
            break

        k += 1

    logger.info('sigma='+str(sigma))

    logger.info('final wavelength solution = '+str(p1[0])+' + '+str(p1[1])+' * pixel + '+str(p1[2])+'* pixel ** 2 +'+\
                str(p1[3])+' * (1/order) + '+str(p1[4])+" * pixel * (1/order) * " + str(p1[5])+" * (pixel ** 2) * (1/order)")

    dataZZ = p1[0] + p1[1] * dataXX + p1[2] * dataXX ** 2 + p1[3] * dataYY + p1[4] * dataXX * dataYY + p1[5] * (
        dataXX ** 2) * dataYY

    if testing:
        ax2.legend()
        pl.show()
    return p1, newoh, dataZZ


def applySolution(order_object, p1):
    if len(p1) > 0:
        newdx = np.arange(1024)
        newy = 1. / order_object.sciorder.order_num
        newoh = np.ravel(
            p1[0] + p1[1] * newdx + p1[2] * newdx ** 2 + p1[3] * newy + p1[4] * newdx * newy + p1[5] * (
                newdx ** 2) * newy)

        # setattr(reduced_order_object.sciorder, "dx_2dfit", astro_math.conv_ang_to_mu(newoh))
        # reduced_order_object.sciorder.dx = astro_math.conv_ang_to_mu(reduced_order_object.sciorder.dx)
        # reduced_order_object.lineobj.matchesohx = astro_math.conv_ang_to_mu(reduced_order_object.lineobj.matchesohx)
        # reduced_order_object.lineobj.bigohx = astro_math.conv_ang_to_mu(reduced_order_object.lineobj.bigohx)
        return newoh
    else:
        return []