import numpy as np
from scipy.interpolate import griddata, LinearNDInterpolator

class EdsViscosity():
    
    def __init__(self):
        import os 
        fdir, fname = os.path.split(__file__)
        self.data = np.genfromtxt(fdir+'/data/WCA_visco_study.dat', 
                                  names=True)
        self.rho = self.data['density'] 
        self.T = self.data['Temperature'] 
        self.mu = self.data['viscosity']
        self.Trho = np.transpose(np.array([self.T,self.rho]))
        self._interp = LinearNDInterpolator(self.Trho, self.mu)
    
    def interpolate(self, T, rho):
        return self._interp(T,rho)

    def contour(self):
        rhomin = np.min(self.rho)
        rhomax = np.max(self.rho)
        Tmin = np.min(self.T)
        Tmax = np.max(self.T)
        xi = np.linspace(Tmin, Tmax,100)
        yi = np.linspace(rhomin,rhomax,100)
        grid_x, grid_y = np.meshgrid(xi, yi)
        grid_z = griddata(self.Trho, self.mu, 
                          (grid_x, grid_y), 
                          method='linear')
        return xi, yi, grid_z
