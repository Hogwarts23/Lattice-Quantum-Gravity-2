#!/anaconda/bin/python

import numpy as np
from scipy.optimize import least_squares



class Fitter(object):
    """
    Orthogonal distance regression (ODR) fitter object.
    """
    
    def __init__(self,):
        pass
    
    def Model(self, function):
        """
        Set the model to fit to.
        """
        self.model = function
        
    def Beta0(self, init_beta):
        """
        Set the initial parameters.
        """
        self.beta0 = np.array(init_beta)
        self.bs = len(init_beta)
        
    def Data(self, xdata, ydata, sigmax=None, sigmay=None, covx=None, covy=None):
        """
        Input the x and y data, as well as the errors if they exist.
        """
        if (len(xdata) != len(ydata)):
            raise ValueError("x and y data must have same shape")
        self.x = xdata
        self.y = ydata
        self.xs = len(xdata)
        self.ys = len(ydata)
        self.delta = np.zeros(shape=(self.xs,))
        if (sigmax is not None) and (sigmay is not None):
            assert (covx is None) and (covy is None)
            Omega = np.block([[np.diag(1./sigmay**2), np.zeros((self.ys, self.ys))],
                              [np.zeros((self.xs, self.xs)), np.diag(1./sigmax**2)]])
        elif (sigmax is not None) and (sigmay is None):
            assert covx is None
            Omega = np.block([[np.eye(self.xs), np.zeros((self.ys, self.ys))],
                              [np.zeros((self.xs, self.xs)), np.diag(1./sigmax**2)]])
        elif (sigmax is None) and (sigmay is not None):
            assert covy is None
            Omega = np.block([[np.diag(1./sigmay**2), np.zeros((self.ys, self.ys))],
                              [np.zeros((self.xs, self.xs)), np.eye(self.xs)]])
        elif (covx is not None) and (covy is not None):
            assert (sigmax is None) and (sigmay is None)
            Omega = np.block([[np.linalg.inv(covy), np.zeros((self.ys, self.ys))],
                              [np.zeros((self.xs, self.xs)), np.linalg.inv(covx)]])
        elif (covx is None) and (covy is not None):
            assert sigmay is None
            Omega = np.block([[np.linalg.inv(covy), np.zeros((self.ys, self.ys))],
                              [np.zeros((self.xs, self.xs)), np.eye(self.xs)]])
        elif (covx is not None) and (covy is None):
            assert sigmax is None
            Omega = np.block([[np.eye(self.xs), np.zeros((self.ys, self.ys))],
                              [np.zeros((self.xs, self.xs)), np.linalg.inv(covx)]])
        else:
            Omega = np.eye(2*self.xs)
        
        self.L = np.linalg.cholesky(Omega)

        
    def Residuals(self, z):
        """
        Calculates the sum of residuals for ODR.
        """
        para = z[:self.bs]
        arguments = z[self.bs:]
        
        rx = arguments
        ry = (self.model(self.x + rx, *para) - self.y)
        r = np.hstack((ry, rx))
        
        return np.dot(r, self.L)
    
    def Run(self,):
        """
        Runs the fitter.
        """
        self.whole = np.append(self.beta0, self.delta)
        self.out = least_squares(self.Residuals, self.whole, method='lm')
        self.chi2 = np.sum(self.out.fun**2) / float(len(self.x)-len(self.beta0))
        