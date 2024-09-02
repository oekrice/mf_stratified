import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.linalg import eigh_tridiagonal

class compute_initial_condition():
    def __init__(self, grid, lbound_fn, run, xbasis = [], ybasis = [], m2 = 0.0, max_mode = 1e6):
        #use the same grid notation as in the proper diagnostic calculator, even though it can be bodged for now
        self.xs = grid.xs
        self.ys = grid.ys
        self.nx = grid.nx
        self.ny = grid.ny

        self.xc = np.zeros(self.nx + 2)
        self.yc = np.zeros(self.ny + 2)

        self.xc[1:-1] = 0.5*(self.xs[1:] + self.xs[:-1])
        self.yc[1:-1] = 0.5*(self.ys[1:] + self.ys[:-1])
        self.xc[0] = self.xc[1] - (self.xc[2] - self.xc[1])
        self.yc[0] = self.yc[0] - (self.yc[2] - self.yc[2])

        self.xc[-1] = self.xc[-2] + (self.xc[-2] - self.xc[-3])
        self.yc[-1] = self.yc[-2] + (self.yc[-2] - self.yc[-3])

        self.dx = np.sum(self.xs[1:] - self.xs[:-1])/len(self.xs[1:])
        self.dy = np.sum(self.ys[1:] - self.ys[:-1])/len(self.ys[1:])

        self.lbound = lbound_fn(self.xc[1:-1])

        self.max_mode = min(max_mode, self.nx)
        self.lbound_transform = np.zeros(self.nx)

        self.xbasis = np.zeros((self.nx, self.nx+2))
        #find eigenvalues and basis vectors (eigenvalues are m^2 and numbered k)
        self.m2, self.xbasis[:,1:-1] = self.find_eigenstuff()
        self.xbasis[:,0] = self.xbasis[:,1]; self.xbasis[:,-1] = self.xbasis[:,-2]

        self.ybasis = self.find_ys(self.m2)

        for k in range(self.nx):   #basis fourier transform on the two boundaries. This step needs doing twice
            self.lbound_transform[k] = self.coeff(self.lbound, k)

        self.phi = self.make_phi()

        self.bxp = (self.phi[1:,1:-1] - self.phi[:-1,1:-1])/self.dx
        self.byp = (self.phi[1:-1,1:] - self.phi[1:-1,:-1])/self.dy

        self.az = self.find_az()
        #plt.pcolormesh(self.xs, self.ys, self.phi[1:-1,1:-1].T)
        #plt.pcolormesh(self.xs, self.yc, self.byp.T)
        #plt.show()
        self.test_phi()

        np.savetxt('./inits/init%03d.txt' % run, self.az.T, delimiter = ',')

    def find_az(self):
        #Finds the vector potential from the magnetic field
        az = np.zeros((self.nx+1,self.ny+1))
        for i in range(1,self.nx+1):
            az[i,0] = az[i-1,0] + self.dx*self.byp[i-1,0]
            for j in range(1,self.ny+1):
                az[i,j] = az[i,j-1] - self.dy*self.bxp[i,j-1]
        return az


    def find_ys(self, eigenvalues):
        #finds suitable (approximately hyperbolic) functions in the y direction
        #CAN SPEED THIS UP WITH ARRAYS!
        ybasis = np.zeros((self.nx, self.ny+2))  #number of x modes plus dimension (with ghosts) in the y direction
        for k in range(0,self.nx-1):  #run through the modes
            ybasis[k][-1] = 1.0   #set to zero at the top (roughly)
            ybasis[k][-2] = -1.0
            m2 = eigenvalues[k]
            for i in range(self.ny, 0, -1):
                ybasis[k][i-1] = m2*ybasis[k][i]*self.dy**2
                ybasis[k][i-1] += 2*ybasis[k][i] - ybasis[k][i+1]
            dfact = (ybasis[k][1] - ybasis[k][0])/self.dy
            ybasis[k] = ybasis[k]/dfact
        return ybasis

    def make_phi(self):
        phi = np.zeros((self.nx+2, self.ny+2))
        for k in range(self.nx-1):
            phi = phi + self.lbound_transform[k]*self.xbasis[k,:][:,np.newaxis]*self.ybasis[k,:][np.newaxis,:]   #lower bound
        return phi

    def find_eigenstuff(self):
        # Uses scipy tridiagonal solver to find the numerical approximations to the sine functions that have the desired properties.
        #Generates a matrix etc. then solves. Should only need to do this once for a given resolution. Doesn't depend on boundary conditions etc.
        d = 2*np.ones(self.nx)
        e = -1*np.ones(self.nx-1)
        d[0] = 1.0; d[-1] = 1.0
        w, v = eigh_tridiagonal(d, e)
        m2 = w/(self.dx**2)
        return m2, v.T

    def second_derivative(self, array, d):
        return (array[:-2] - 2*array[1:-1] + array[2:])/d**2

    def fcheck(self, bound_trans):
        bcheck = 0.0*ubound
        for k in range(self.nx -1):
            bcheck = bcheck + bound_trans[k]*self.xbasis[k,1:-1]
        return bcheck

    def mode(self, m):
        return np.sin(0.5*m*np.pi*self.xc/self.xs[-1])

    def coeff(self, bound, k):
        return np.sum(bound*self.xbasis[k,1:-1])/np.sum(self.xbasis[k,1:-1]*self.xbasis[k,1:-1])

    def test_phi(self):
        d2x = (self.phi[:-2,1:-1] - 2*self.phi[1:-1,1:-1] + self.phi[2:,1:-1])/self.dx**2
        d2y = (self.phi[1:-1,:-2] - 2*self.phi[1:-1,1:-1] + self.phi[1:-1,2:])/self.dy**2
        print('Max laplacian', np.max(np.abs(d2x + d2y)))
        lbound_test = (self.phi[:,1] - self.phi[:,0])/self.dy

        print('Max lbound error', np.max(np.abs(self.lbound - lbound_test[1:-1])))
        print('Max left side flux', np.max(np.abs(self.phi[1,1:-1] - self.phi[0,1:-1])))
        print('Max right side flux', np.max(np.abs(self.phi[-1,1:-1] - self.phi[-2,1:-1])))

        jz =  (self.byp[1:,1:-1] - self.byp[:-1,1:-1])/self.dx - (self.bxp[1:-1,1:] - self.bxp[1:-1,:-1])/self.dy
        print('Max perpendicular current', np.max(np.abs(jz)))

class make_test_grid():
    def __init__(self, nx, ny):
        self.nx = nx + 1
        self.ny = ny + 1   #number of grid points in the main block
        self.data = [[],[]]
        self.data[0] = np.linspace(-1.,1.,self.nx)
        self.data[1] = np.linspace(0.,1.,self.ny)

