# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 09:52:38 2014

@author: jholt
"""
import itertools

import numpy as np
from numpy import fft
import scipy as sp
from scipy import optimize

import fits
import robust


try:
    from scipy.signal import argrelextrema
except:
    print 'need to update scipy to get argrelextrema'

# collection of array math utilities that do not need self.data as in am classes
def conv_ang_to_mu(dx):
    mu_dx = dx / 10000.
    return mu_dx


def median_comb(cslist):
    """ stack up all arrays or fits files in cslist and take the median
    :param cslist:
    Parameters:
    --------------
    cslist : data list or array containing string FITS file names
    Returns:
    -------------
    median combined numpy array     
    """
    data_array = []
    for name in cslist:
        data = fits.Handle_fits.ensure_array(name)
        data_array.append(data)

    # images is a comma separated list of files to be median combined
    all_images = np.array(data_array)
    med_comb = np.median(all_images, axis=0)

    return med_comb


def fit_gauss(array, sigma=0.2):
    """ fit a gaussian to a 1D numpy array 
    Parameters:
    -----------------------
    array: numpy array
    sigma: float
        starting sigma to use in gaussian equation least squares fit
        
    Returns:
    gaussian fit parameters : p[0] * e** (-(x-p[1]**2 / (2*p[2]**2)) / p[2]*sqrt(2pi)
    fit_function: the fit parameters applied to an array of length paramter array
    
    """
    x = np.arange(len(array))
    if not x:
        return '', ''
    fit_func = lambda p, x: (p[0] / (p[2] * np.sqrt(2 * np.pi))) * np.exp(-(x - p[1]) ** 2 / 2 * p[2] ** 2)
    errfunc = lambda p, x, y: fit_func(p, x) - y
    p0 = [array.max(), array.argmax(), sigma]
    fit = sp.optimize.leastsq(errfunc, p0, args=(x, array))

    return fit, fit_func(fit[0], x)


def cheb_fit(matchesohx, matchesdx, dx):
    """
          not used
        """
    if len(matchesohx) > 4:
        chebord = 3
    else:
        chebord = 2

    # c is Legendre series coefficients, chebord is degree of the fitting serious
    c, [residuals, rank, singular_values, rcond] = np.polynomial.chebyshev.chebfit(matchesdx, matchesohx, chebord,
                                                                                   full=True)

    # find Chebyshev fit with no Rank error and lowest order
    while singular_values.min() < 0.000000000000001 and chebord > 0:
        chebord -= 1
        c, [residuals, rank, singular_values, rcond] = np.polynomial.chebyshev.chebfit(matchesdx, matchesohx, chebord,
                                                                                       full=True)

        # apply coefficients to make a chebyshev function
    dispcx = np.polynomial.chebyshev.chebval(matchesdx, c)

    # fits a line to the chebyshev solution, this is only used to find and remove bad fits
    # disp_sol,disp,chebtoline_residual=linearize(dispcx)

    # if float(disp) < 0.1 or float(disp) > 0.9 or abs(chebtoline_residual) > 400000: 
    # return 'nofit'

    return dispcx


def linearize(dispcx):
    """ called from cheb fit, which is not used"""

    p, residuals, rank, singular_values, rcond = np.polyfit(np.arange(1024), dispcx, 1, full=True)
    disp_sol = p[0] * np.arange(1024) + p[1]

    return disp_sol, p[0], residuals


def arg_max_corr(a, b):
    """ Find the maximum of the cross-correlation - includes upsampling
    """
    if len(a.shape) > 1:
        raise ValueError('Needs a 1-dimensional array.')
    length = len(a)
    if not length % 2 == 0:
        raise ValueError('Needs an even length array.')

    if not a.shape == b.shape:
        raise ValueError('The 2 arrays need to be the same shape')

    # Start by finding the coarse discretised arg_max
    coarse_max = np.argmax(np.correlate(a, b, mode='full')) - length + 1

    omega = np.zeros(length)
    omega[0:length / 2] = (2 * np.pi * np.arange(length / 2)) / length
    omega[length / 2 + 1:] = (2 * np.pi *
                              (np.arange(length / 2 + 1, length) - length)) / length

    fft_a = fft.fft(a)

    def correlate_point(tau):
        rotate_vec = np.exp(1j * tau * omega)
        rotate_vec[length / 2] = np.cos(np.pi * tau)

        return np.sum((fft.ifft(fft_a * rotate_vec)).real * b)

    start_arg, end_arg = (float(coarse_max) - 1, float(coarse_max) + 1)

    max_arg = optimize.fminbound(lambda tau: -correlate_point(tau),
                                 start_arg, end_arg)
    # print 'coarse_max=',coarse_max,' max_arg=',max_arg

    return max_arg


def fit_poly(cm, xes='default', deg=4):
    """
    fit a polynomial of degree=degree to array cm (usually output from spectroid)
    
    """
    p0 = np.polyfit(np.arange(len(cm) - 10), cm[:-10], deg=deg)  # end always drops off
    if xes == 'default':
        xes = np.arange(len(cm))
    if deg == 4:
        cmfit = p0[0] * xes ** 4 + p0[1] * xes ** 3 + p0[2] * xes ** 2 + p0[3] * xes + p0[4]
    elif deg == 3:
        cmfit = p0[0] * xes ** 3 + p0[1] * xes ** 2 + p0[2] * xes + p0[3]
    elif deg == 2:
        cmfit = p0[0] * xes ** 2 + p0[1] * xes + p0[2]
    elif deg == 1:
        cmfit = p0[0] * xes + p0[1]
    else:
        return p0
    return cmfit, p0


def fit_plane_ltsq(X, y, z):
    # Fits a plane to a point cloud, 
    # Where z = aX + bY + c        ----Eqn #1
    # Rearanging Eqn1: aX + bY -z +c =0
    # Gives normal (a,b,-1)
    # Normal = (a,b,-1)
    # [rows,cols] = XYZ.shape
    rows = len(X)
    cols = 3
    # print"X=",X
    # print"Y=",Y
    # print"z=",z
    # G = np.ones((rows,3))
    # G[:,0] = XYZ[:,0]  #X
    # G[:,1] = XYZ[:,1]  #Y
    # z = XYZ[:,2]
    # (a,b,c),resid,rank,s = np.linalg.lstsq(G,z)
    G = np.ones((rows, cols))
    G[:, 0] = X
    G[:, 1] = y

    (a, b, c), resid, rank, s = np.linalg.lstsq(G, z)

    # print'plane fit coeffs,resid,rank,s=',(a,b,c),resid,rank,s
    return a, b, c


def plotplane(coeffs, X, Y, Z):
    """

    :param coeffs:
    :param X:
    :param Y:
    :param Z:
    """
    import matplotlib.pyplot as plt
    # from matplotlib import cm    

    axs = np.array(X)
    xsurf = np.arange(axs.min() - 1, axs.max() + 1, 5)

    orders = np.array(Y)
    ysurf = np.arange(orders.min() - 0.005, orders.max() + 0.005, 1)
    x_surf, y_surf = np.meshgrid(xsurf, ysurf)
    z_surf = coeffs[0] * x_surf + coeffs[1] * y_surf + coeffs[2]

    fig = plt.figure()
    plt.clf()
    plt.hold(True)
    ax = fig.gca(projection='3d')
    ax.plot_wireframe(x_surf, y_surf, z_surf)
    ax.scatter(X, Y, Z, color='red')
    ax.set_xlabel('pixel ')
    ax.set_ylabel('order')
    ax.set_zlabel('wavlength')

    plt.show()


def find_peaks(skys, indices):
    """ not used for now, finds the peaks in the array 'skys' around the indices
    in indices"""
    new_indices = []
    for idx in indices:
        locpeak = argrelextrema(skys[idx - 2:idx + 2], np.greater)
        newloc = locpeak + idx - 2

        new_indices.append(newloc)

    return new_indices


def fit3dline(XYZ):
    """
    fits a 3d line to input file XYZ (order, disp, zeropoint)
    
    """
    datamean = XYZ.mean(axis=0)
    # print'datamean=',datamean

    uu, dd, vv = np.linalg.svd(XYZ - datamean)
    # vv[0] is the unit direction vector of the best fit line

    xs = []
    ys = []
    zs = []
    for point1 in XYZ:
        xs.append(point1[0])
        ys.append(point1[1])
        zs.append(point1[2])

    # parametric equations of line
    t = np.arange(-1000, 1000)
    x = datamean[0] + t * vv[0][0]
    y = datamean[1] + t * vv[0][1]
    z = datamean[2] + t * vv[0][2]

    return vv[0], datamean, [x, y, z]


def newpoint(order, three_d_slope, datamean):
    """
    finds a newpoint using the linear 2d wavelength solution found in fit3dline
    """

    # using parametric equations of line

    newx = order
    newtee = (newx - datamean[0]) / three_d_slope[0]
    newy = datamean[1] + newtee * three_d_slope[1]
    newz = datamean[2] + newtee * three_d_slope[2]

    newpoint = [newx, newy, newz]
    return newpoint


def plot3dline(XYZ, newpoint):
    datamean = XYZ.mean(axis=0)
    uu, dd, vv = np.linalg.svd(XYZ - datamean)
    # vv[0] is the unit direction vector of the best fit line

    xs = []
    ys = []
    zs = []
    for point1 in XYZ:
        xs.append(point1[0])
        ys.append(point1[1])
        zs.append(point1[2])

    # parametric equations of line
    t = np.arange(-1000, 1000)
    x = datamean[0] + t * vv[0][0]
    y = datamean[1] + t * vv[0][1]
    z = datamean[2] + t * vv[0][2]

    import matplotlib.pyplot as plt

    fig = plt.figure()
    plt.clf()
    plt.hold(True)
    ax = fig.gca(projection='3d')

    nxs = np.array(xs)
    nys = np.array(ys)
    nzs = np.array(zs)

    [newx, newy, newz] = newpoint

    ax.scatter(xs, ys, zs, color='red')
    ax.scatter(newx, newy, newz, color='blue')
    ax.plot(x, y, z, color='black')
    ax.set_xlabel('order')
    ax.set_ylabel('dispersion')
    ax.set_zlabel('zeropoint')
    ax.set_xlim(nxs.min() - 1, nxs.max() + 2)
    ax.set_ylim(nys.min() - 0.5, nys.max() + 0.5)
    ax.set_zlim(nzs.min() - 100, nzs.max() + 100)
    plt.show()


def plot_all_disp(XYZ):
    import pylab as pl

    xs = []
    ys = []
    zs = []
    for point1 in XYZ:
        xs.append(point1[0])
        ys.append(point1[1])
        zs.append(point1[2])

    pl.figure(1)
    pl.clf()
    pl.plot(zs, xs)
    pl.plot(zs, xs, 'r.')
    pl.title('order vs zeropoint')

    pl.figure(2)
    pl.clf()
    pl.plot(ys, xs)
    pl.plot(ys, xs, 'r.')
    pl.title('order vs disp')

    pl.show()


def fit_curve(X, Y, Z):
    from scipy.optimize import curve_fit

    # ensure X,Y,Z are numpy arrays
    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)

    if len(Z) > 4:
        def func(w, a, b, c, d, e):
            return a * w[0] ** 2 + b * w[0] + c * w[1] ** 2 + d * w[1] + e
    else:
        def func(w, a, b, c):
            return a * w[0] + b * w[1] + c

            # 
            # xdata=np.array(xs)
            # ydata=np.array(ys)    
            # zdata=np.array(zs)
    xys = np.array([X, Y])

    popt, pcov = curve_fit(func, xys, Z)

    if len(Z) > 4:
        disp_set_fit = popt[0] * X ** 2 + popt[1] * X + popt[2] * Y ** 2 + popt[3] * Y + popt[4]
    else:
        disp_set_fit = popt[0] * X + popt[1] * Y + popt[2]

    # print 'popt=',popt
    # printX,Y,Z,disp_set_fit
    import pylab as pl

    pl.figure(1)
    pl.clf()
    pl.plot(X, Z, 'r*')
    pl.plot(X, disp_set_fit, 'r')

    pl.figure(2)
    pl.clf()
    pl.plot(Y, Z, 'r*')
    pl.plot(Y, disp_set_fit, 'r')
    pl.show()
    return popt, X, Y, Z, disp_set_fit


def newpoint2(popt, new_order):
    new_zeropoint = extend(XYZ, new_order)
    new_disp = popt[0] * new_order ** 2 + popt[1] * new_order + popt[2] * new_zeropoint ** 2 + \
               popt[3] * new_zeropoint + popt[4]

    return new_disp, new_zeropoint


def extend(XYZ, neworder):
    import pylab as pl
    from scipy.interpolate import InterpolatedUnivariateSpline

    xs = []
    ys = []
    zs = []
    for point1 in XYZ:
        xs.append(point1[0])
        ys.append(point1[1])
        zs.append(point1[2])

    xdata = np.array(xs)
    ydata = np.array(ys)
    zdata = np.array(zs)

    xdata = xdata[::-1]
    ydata = ydata[::-1]
    zdata = zdata[::-1]
    xvals = np.arange(xdata.min() - 2, xdata.max() + 2, 0.1)

    # k=1 is linear spline 2 quadratic 3 cubic
    s = InterpolatedUnivariateSpline(xdata, zdata, k=1)
    z = s(xvals)

    # k=1 is linear spline 2 quadratic 3 cubic
    s2 = InterpolatedUnivariateSpline(xdata, ydata, k=1)
    y = s2(xvals)

    newz = s(neworder)
    # print'new_zeropoint=',newz
    newy = s2(neworder)
    # print'new_disp=',newy
    # print'new_order=',neworder

    pl.figure(2)
    pl.clf()
    pl.plot(xdata, zdata, 'bo')
    pl.plot(xvals, z, '-x')
    pl.plot(neworder, newz, 'ro')
    pl.title('extrapolated zeropoints vs order')
    pl.show()

    pl.figure(3)
    pl.clf()
    pl.plot(xdata, ydata, 'bo')
    pl.plot(xvals, y, '-x')
    pl.plot(neworder, newy, 'ro')
    pl.title('extrapolated dispersion vs order')
    pl.show()

    return newy, newz


def actual_to_theory(loc1, loc2, threshold='40'):
    if abs(loc1 - loc2) < threshold and loc1 > 1:
        return True
    else:
        return False


def order_wavelength_solution(matchesdx, matchesohx, dx, order):
    # print'in wavelength solution matchesdx=',matchesdx
    # print'matches ohx = ',matchesohx
    p4 = np.poly1d(np.polyfit(matchesdx, matchesohx, 4))
    p3 = np.poly1d(np.polyfit(matchesdx, matchesohx, 3))
    p2 = np.poly1d(np.polyfit(matchesdx, matchesohx, 2))
    p1 = np.poly1d(np.polyfit(matchesdx, matchesohx, 1))
    p20 = np.poly1d(np.polyfit(matchesdx, matchesohx, 20))
    # c1 = cheb_fit(matchesohx, matchesdx, dx)

    # v4=np.polyval(np.polyfit(matchesdx, matchesohx, 4),matchesdx)
    v3 = np.polyval(np.polyfit(matchesdx, matchesohx, 3), matchesdx)
    v2 = np.polyval(np.polyfit(matchesdx, matchesohx, 2), matchesdx)
    v1 = np.polyval(np.polyfit(matchesdx, matchesohx, 1), matchesdx)
    v20 = np.polyval(np.polyfit(matchesdx, matchesohx, 20), matchesdx)

    cc2 = robust.polyfit(matchesdx, matchesohx, 2, iterMax=25)
    pc2 = np.poly1d(cc2)
    vc2 = np.polyval(cc2, matchesdx)

    f1, residuals1, rank1, singular_values1, rcond1 = np.polyfit(matchesdx, matchesohx, 1, full=True)
    f2, residuals2, rank2, singular_values2, rcond2 = np.polyfit(matchesdx, matchesohx, 2, full=True)
    f3, residuals3, rank3, singular_values3, rcond3 = np.polyfit(matchesdx, matchesohx, 3, full=True)
    f4, residuals4, rank4, singular_values4, rcond4 = np.polyfit(matchesdx, matchesohx, 4, full=True)
    f20, residuals20, rank20, singular_values20, rcond20 = np.polyfit(matchesdx, matchesohx, 20, full=True)

    # print'1 results=',f1, residuals1, rank1, singular_values1, rcond1
    # print'2 results=',f2, residuals2, rank2, singular_values2, rcond2 
    # print'3 results=',f3, residuals3, rank3, singular_values3, rcond3
    # print'4 results=',f4, residuals4, rank4, singular_values4, rcond4 
    # print'20 results=', residuals20, rank20, singular_values20, rcond20

    import pylab as pl

    pl.figure(1)
    pl.clf()
    pl.plot(matchesdx, matchesohx, 'k*')
    pl.plot(matchesdx, p4(matchesdx), 'r-', label='4th order chisq=' + str(residuals4))
    pl.plot(matchesdx, p1(matchesdx), 'g-', label='1st order chisq=' + str(residuals1))
    pl.plot(matchesdx, p2(matchesdx), 'y-', label='2nd order chisq=' + str(residuals2))
    pl.plot(matchesdx, p3(matchesdx), 'c-', label='3rd order chisq=' + str(residuals3))
    pl.plot(matchesdx, p20(matchesdx), 'm-', label='20th order chisq=' + str(residuals20))
    pl.plot(matchesdx, pc2(matchesdx), 'b-', label='cleaned 2 order fit')
    pl.xlabel('actual sky line locations')
    pl.ylabel('sky line OH list values')
    pl.legend(loc=4)
    pl.show()

    chi1 = np.sum((np.polyval(f1, matchesdx) - matchesohx) ** 2)
    chi2 = np.sum((np.polyval(f2, matchesdx) - matchesohx) ** 2)
    chi3 = np.sum((np.polyval(f3, matchesdx) - matchesohx) ** 2)
    chi4 = np.sum((np.polyval(f4, matchesdx) - matchesohx) ** 2)
    chic2 = np.sum((np.polyval(cc2, matchesdx) - matchesohx) ** 2)
    chi20 = np.sum((np.polyval(f20, matchesdx) - matchesohx) ** 2)
    # print'chis=',chi1,chi2,chi3,chi4,chi20,chic2

    pl.figure(2)
    pl.clf()
    pl.plot(matchesdx, v1 - matchesohx, 'g-', label='1st order  chisq=' + str(chi1))
    pl.plot(matchesdx, v2 - matchesohx, 'y-', label='2nd order  chisq=' + str(chi2))
    pl.plot(matchesdx, v3 - matchesohx, 'c-', label='3rd order  chisq=' + str(chi3))
    pl.plot(matchesdx, vc2 - matchesohx, 'k-', label='2 order cleaned  chisq=' + str(chic2))

    pl.plot(matchesdx, v20 - matchesohx, 'm-', label='20th chisq=' + str(chi20))
    pl.title('residuals')
    pl.legend(loc=4)

    # this does nothing later as we apply the 2d solution to np.arange(1024)
    dx = p3(dx)

    # find the real dispersion solution, with at x=(1,1024)
    n20 = np.poly1d(np.polyfit(np.arange(1024), dx, 20))
    n3 = np.poly1d(np.polyfit(np.arange(1024), dx, 3))
    n2 = np.poly1d(np.polyfit(np.arange(1024), dx, 2))
    n1 = np.poly1d(np.polyfit(np.arange(1024), dx, 1))

    nf1, nresiduals1, nrank1, nsingular_values1, nrcond1 = np.polyfit(np.arange(1024), dx, 1, full=True)
    nf2, nresiduals2, nrank2, nsingular_values2, nrcond2 = np.polyfit(np.arange(1024), dx, 2, full=True)
    nf3, nresiduals3, nrank3, nsingular_values3, nrcond3 = np.polyfit(np.arange(1024), dx, 3, full=True)
    nf4, nresiduals4, nrank4, nsingular_values4, nrcond4 = np.polyfit(np.arange(1024), dx, 4, full=True)
    nf20, nresiduals20, nrank20, nsingular_values20, nrcond20 = np.polyfit(np.arange(1024), dx, 20, full=True)

    pl.figure(3)
    pl.clf()
    # print'dx.shape=',dx.shape
    pl.plot(np.arange(1024), dx, 'k.')
    pl.plot(np.arange(1024), n1(np.arange(1024)), 'g-', label='1st order chisq=' + str(nresiduals1))
    pl.plot(np.arange(1024), n2(np.arange(1024)), 'y-', label='2nd order chisq=' + str(nresiduals2))
    pl.plot(np.arange(1024), n3(np.arange(1024)), 'c-', label='3rd order chisq=' + str(nresiduals3))
    pl.plot(np.arange(1024), n20(np.arange(1024)), 'm-', label='20th order chisq=' + str(nresiduals20))
    pl.legend(loc=4)
    pl.show()

    # dx = nf3[0]*new_x**3 + nf3[1]*new_x**2 + nf3[2]*new_x + nf3[3]
    # dx = nf2[0]*new_x**2 ++ nf2[1]*new_x + nf2[2]

    #dx = nf1[0]*new_x + nf1[1]

    #print'nf2=',nf2

    #orig_indices=[]
    #for match in matchesdx:            
    #    orig_indices.append(np.abs(sciorder.dx - match).argmin())

    oh_indices = []
    for omatch in matchesohx:
        oh_indices.append(np.abs(sciorder.dx - omatch).argmin())

        #new_indices = find_peaks(skys, orig_indices)    
    #real_disp, real_offset = nirspec_wavelength_utils.find_orig_lambda_solution(dx)
    #disp_sol = ([self.order, real_disp, real_offset])
    #print'nresiduals2=',nresiduals2
    #print'len(matchesohx)=',len(matchesohx)
    #new_list = ([order, matchesohx, orig_indices])              

    if float(nresiduals2) > 500.:
        found_wavelength = False
        bad_sol = ([order, nf1[0], nf1[1]])

        #print'fw=false'
    #elif float(nresiduals2) > 100. and len(matchesohx) < 7:
    #    found_wavelength=False
    #    bad_sol = ([self.order, nf1[0], nf1[1]])
    else:
        #disp_sol = ([self.order, nf3[0], nf3[1], nf3[2], nf3[3]])
        #disp_sol = ([self.order, nf2[0], nf2[1], nf2[2]])
        disp_sol = ([order, nf1[0], nf1[1]])
        #new_list = ([order, matchesohx, orig_indices])              
        #print'good solution for order='+str(order)

    #print'sobjdx min and max=', dx.min(), dx.max()
    matched = order, matchesdx, matchesohx

    found_wavelength = True


# m2 = fit2dpoly(dxs, order_array, ohs, deg=2)

def fit2dpoly_noXterms(X, Y, Z, deg=3):
    ncols = (2 * deg) + 1
    G = np.zeros((X.size, ncols))
    # ij = itertools.product(range(deg + 1), range(deg + 1))  
    if deg == 2:
        G = np.zeros((X.size, ncols))
        G[:, 0] = X ** 0 * Y ** 0
        G[:, 1] = X ** 0 * Y ** 1
        G[:, 2] = X ** 0 * Y ** 2
        G[:, 3] = X ** 1 * Y ** 0
        # G[:,4] = X**1 * Y**1
        # G[:,5] = X**1 * Y**2
        G[:, 4] = X ** 2 * Y ** 0
        # G[:,7] = X**2 * Y**1   
        # G[:,8] = X**2 * Y**2
    elif deg == 1:
        G = np.zeros((X.size, ncols))
        G[:, 0] = X ** 0 * Y ** 0
        G[:, 1] = X ** 0 * Y ** 1
        G[:, 2] = X ** 1 * Y ** 0

    # if deg==2: print 'G=',G,' X=',X
    m, _, _, _ = np.linalg.lstsq(G, Z)
    return m


def polyval2d_noXterms(X, Y, m):
    deg = (len(m) - 1) / 2
    print 'deg=', deg
    # ij = itertools.product(range(deg+1), range(deg+1))
    Z = np.zeros_like(X)
    if deg == 2:
        ijs = [[0, 0], [0, 1], [0, 2], [1, 0], [2, 0]]
    elif deg == 1:
        ijs = [[0, 0], [0, 1], [1, 0]]
    for a, (i, j) in zip(m, ijs):
        # if i==0 or j==0: #try take out cross-terms
        Z += a * X ** i * Y ** j
    return Z


def second_orderpoly_inpix(a0, a1, a2, a3, a4, a5):
    return lambda X, Y: a0 + a1 * X + a2 * X ** 2 + a3 * Y + a4 * X * Y + a5 * (X ** 2) * Y


def fit2ndorderpoly_inpix(oh):
    params = [20000, 0.3, 0., 20000., 0., 0.]

    errfct = lambda p: ravel(second_orderpoly_inpix(*p)(*indices(oh.shape)) - oh)

    pfit, pcov, infodict, errmsg, success = optimize.leastsq(errfct, params, full_output=1)

    if (len(oh) > len(params)) and pcov is not None:
        s_sq = (errfct(pfit) ** 2).sum() / (len(oh) - len(params))
        print 's_sq=', s_sq
        pcov *= s_sq
    else:
        pcov = inf

    error = []
    for i in range(len(pfit)):
        try:
            error.append(np.absolute(pcov[i][i]) ** 0.5)
        except:
            error.append(0.00)
    pfit_leastsq = pfit
    perr_leastsq = np.array(error)
    print 'pfit=', pfit_leastsq
    print 'perr=', perr_leastsq

    return pfit, success


def full_fit(dx, order_num, oh):
    data = second_orderpoly_inpix()
    p, success = fit2ndorderpoly_inpix(oh)
    # fit = second_orderpoly_inpix(*p)
    # return fit


# def fit_function(p0, datax, datay, dataz, function):
# errfunc = lambda p, x, y, z: function(x,y,p)-y

def fit2dpoly(X, Y, Z, deg=3):
    ncols = (deg + 1) ** 2
    G = np.zeros((X.size, ncols))
    ij = itertools.product(range(deg + 1), range(deg + 1))
    for k, (i, j) in enumerate(ij):
        if deg == 2:
            print k, i, j
        # if i==0 or j==0: #try take out cross-terms
        G[:, k] = X ** i * Y ** j
    # if deg==2: print 'G=',G,' X=',X
    m, _, _, _ = np.linalg.lstsq(G, Z)
    return m


def polyval2d(X, Y, m):
    deg = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(deg + 1), range(deg + 1))
    Z = np.zeros_like(X)
    # ijs=[[0,0],[0,1],[0,2],[1,0],[2,0]]
    for a, (i, j) in zip(m, ij):
        # if i==0 or j==0: #try take out cross-terms
        Z += a * X ** i * Y ** j
    return Z


def wavelength_solution_2d(all_matches):
    if len(all_matches) > 0:

        # all_matches.append((lineobj.matchesdx, lineobj.matchesohx, lineobj.matchesidx, order))
        nord = []
        ohs = []
        dxs = []
        order_array = []
        for point1 in all_matches:
            ohs.append(point1[1])
            dxs.append(point1[2])
            nord.append(point1[3])
            order_array.append((np.zeros(len(point1[2])) + 1. / point1[3]))

        order_array = np.hstack(np.array(order_array))
        ohs = np.hstack(np.array(ohs))
        dxs = np.hstack(np.array(dxs))

        print 'order_array=', order_array, len(order_array)
        print 'ohs=', ohs, len(ohs)
        print 'dxs=', dxs, len(dxs)

        # ohs (real sky line values) is the Z in the fit
        coeffs = fit_plane_ltsq(dxs, order_array, ohs)
        plotplane(coeffs, dxs, order_array, ohs)

        m1 = fit2dpoly(dxs, order_array, ohs, deg=1)
        m2 = fit2dpoly(dxs, order_array, ohs, deg=2)
        print 'm1=', m1
        print 'm2=', m2
        m1_noXterms = fit2dpoly_noXterms(dxs, order_array, ohs, deg=1)
        m2_noXterms = fit2dpoly_noXterms(dxs, order_array, ohs, deg=2)
        m3 = fit2dpoly(dxs, order_array, ohs, deg=3)
        print 'm1_noX=', m1_noXterms
        print 'm2_noX=', m2_noXterms
        m4 = fit2dpoly(dxs, order_array, ohs, deg=4)

        nx, ny = 100, 100
        xx, yy = np.meshgrid(np.linspace(dxs.min() - 5, dxs.max() + 5, nx),
                             np.linspace(order_array.min() - 0.005, order_array.max() + 0.005, ny))

        zz1 = polyval2d(xx, yy, m1)
        zz2 = polyval2d(xx, yy, m2)
        zz3 = polyval2d(xx, yy, m3)

        import matplotlib.pyplot as plt

        fig = plt.figure(1)
        plt.clf()
        plt.hold(True)
        ax = fig.gca(projection='3d')
        ax.plot_wireframe(xx, yy, zz2)
        ax.scatter(dxs, order_array, ohs, color='red')
        ax.set_xlabel('pixel ')
        ax.set_ylabel('1/order')
        ax.set_zlabel('wavlength')
        fig = plt.figure(2)
        plt.clf()
        plt.hold(True)
        ax = fig.gca(projection='3d')
        ax.plot_wireframe(xx, yy, zz3)
        ax.scatter(dxs, order_array, ohs, color='red')
        ax.set_xlabel('pixel ')
        ax.set_ylabel('1/order')
        ax.set_zlabel('wavlength')
        fig = plt.figure(3)
        plt.clf()
        plt.hold(True)
        ax = fig.gca(projection='3d')
        ax.plot_wireframe(xx, yy, zz1)
        ax.scatter(dxs, order_array, ohs, color='red')
        ax.set_xlabel('pixel ')
        ax.set_ylabel('1/order')
        ax.set_zlabel('wavlength')
        plt.show()

        import pylab as pl

        pl.figure(17)
        pl.clf()
        i = 0
        colors = ['r', 'g', 'c', 'k', 'm', 'b', 'y', 'r', 'g', 'c', 'k', 'm', 'b', 'y', 'r', 'g', 'c', 'k', 'm', 'b',
                  'y']
        # all_matches.append((lineobj.matchesdx, lineobj.matchesohx, lineobj.matchesidx, order))

        while i < len(all_matches):
            pl.plot(all_matches[i][2], all_matches[i][1], '*' + str(colors[i]), label='m=' + str(all_matches[i][3]))
            pl.plot(all_matches[i][2], all_matches[i][1], str(colors[i]))
            pl.title('actual sky line values vs their pixel location for each order')
            pl.legend()
            i += 1
        pl.show()

        i = 0
        while i < len(all_matches):
            dx_o = np.array(all_matches[i][2])
            oh_o = np.array(all_matches[i][1])
            order_o = all_matches[i][3]
            order_array = np.zeros(len(dx_o)) + 1. / order_o

            # try seeing how the slope is changing
            import robust

            cc2 = robust.polyfit(dx_o, oh_o, 2, iterMax=25)
            cc3 = robust.polyfit(dx_o, oh_o, 3, iterMax=25)

            # print 'fits to just the oh,orig pix index order='+str(order_o)
            print 'cc2=', cc2
            print 'cc3=', cc3

            # this is the linear solution:
            new_oh_o = coeffs[0] * dx_o + coeffs[1] * order_array + coeffs[2]

            new_oh_o1 = polyval2d(dx_o, order_array, m1)
            new_oh_o1_noXterms = polyval2d_noXterms(dx_o, order_array, m1_noXterms)
            new_oh_o2 = polyval2d(dx_o, order_array, m2)
            new_oh_o2_noXterms = polyval2d_noXterms(dx_o, order_array, m2_noXterms)
            new_oh_o3 = polyval2d(dx_o, order_array, m3)
            new_oh_o4 = polyval2d(dx_o, order_array, m4)
            # new_oh_o10 = polyval2d(dx_o,order_array,m10)
            # new_oh_o20 = polyval2d(dx_o,order_array,m20)

            # print'coeffs=',coeffs
            # print'old oh=',oh_o

            # print'm1 = ',m1
            print 'oldo=', oh_o
            print'new1a=', new_oh_o
            print'new1=', new_oh_o1
            print'new1_noXterms=', new_oh_o1_noXterms

            print'new2 = ', new_oh_o2
            print 'new2_noXterms=', new_oh_o2
            print'new3 = ', new_oh_o3
            print'new4 = ', new_oh_o4

            chi1 = np.sum((new_oh_o1 - oh_o) ** 2)
            chi1_noXterms = np.sum((new_oh_o1_noXterms - oh_o) ** 2)
            chi2 = np.sum((new_oh_o2 - oh_o) ** 2)
            chi2_noXterms = np.sum((new_oh_o2_noXterms - oh_o) ** 2)
            chi3 = np.sum((new_oh_o3 - oh_o) ** 2)
            chi4 = np.sum((new_oh_o4 - oh_o) ** 2)
            # chi10 = np.sum((new_oh_o10-oh_o)**2)
            # chi20 = np.sum((new_oh_o20-oh_o)**2)

            print'1 residuals=', chi1
            print'1 residuals_noXterms=', chi1_noXterms
            print'2 residuals=', chi2
            print'2 residuals noXterms=', chi2_noXterms
            print'3 residuals=', chi3
            print'4 residuals=', chi4
            i += 1
        return m2, m2_noXterms, m1, m1_noXterms
    return np.zeros(1024), np.zeros(1024)


def apply_2d_wavelength_solution(p0, order, all_disp_sol):
    new_zpt = p0[0] * order + p0[1]

    # new_disp = disp_func(order, new_zpt)
    # all_disp_sol = np.append(all_disp_sol,[order, new_disp, new_zpt])
    # print'allds=',all_disp_sol
    popt, order_data, zpt_data, disp_data = fit_curve(all_disp_sol)

    # apply the 2d lambda solution
    # bad solution
    # temp, disp, offset = newpoint(order, three_d_slope, datamean)
    # bad solution

    # real_disp, real_offset = nirspec_wavelength_utils.find_orig_lambda_solution(twod_fit)
    # print 'ads=', all_disp_sol
    # print 'bad linear solution=', temp, disp,offset
    #bad_twod_fit = disp * np.arange(1024) + offset
    # uncomment to show linear (bad fit)
    #dx = bad_twod_fit

    #plot3dline(all_disp_sol, [order,disp,offset])

    #new_disp, new_offset = extend(all_disp_sol, order)

    #twod_fit = new_disp * np.arange(1024) + new_offset
    # uncomment to show extrapolated fit
    #dx = twod_fit

    #self.logger.info('Wavelength dispersion '+str(new_disp)+' and offset '+str(new_offset)+' applied to
    #  order '+str(order))
