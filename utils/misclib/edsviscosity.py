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
        self.Re = self.data['Re']
        self.Trho = np.transpose(np.array([self.T,self.rho]))
        self._interp = LinearNDInterpolator(self.Trho, self.mu)
    
    def interpolate(self, T, rho):
        return self._interp(T,rho)

    def contour(self,res=100,conttype='mu'):
        rhomin = np.min(self.rho)
        rhomax = np.max(self.rho)
        Tmin = np.min(self.T)
        Tmax = np.max(self.T)
        xi = np.linspace(Tmin, Tmax,res)
        yi = np.linspace(rhomin,rhomax,res)
        grid_x, grid_y = np.meshgrid(xi, yi)

        if conttype is 'mu':
            plotvar = self.mu
        elif conttype is 'Re':
            plotvar = self.Re

        grid_z = griddata(self.Trho, plotvar, 
                          (grid_x, grid_y), 
                          method='linear')
        return xi, yi, grid_z


if __name__ == "__main__":

    import matplotlib.pyplot as plt
    import matplotlib

    obj = EdsViscosity()

    #Plot viscosity
    x, y, mu = obj.contour()
    plt.pcolormesh(x, y, mu,vmin=0.1,vmax=100.,norm=matplotlib.colors.LogNorm(), cmap=plt.cm.RdYlBu_r); 
    plt.colorbar(); 
    plt.show()

    #Plot Reynolds
    x, y, Re = obj.contour(conttype='Re')    
    plt.pcolormesh(x, y, Re,vmin=0.0,vmax=2., cmap=plt.cm.RdYlBu_r); 
    plt.colorbar(); 
    plt.show()

