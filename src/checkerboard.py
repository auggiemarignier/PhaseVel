import numpy as np
import healpy as hp
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt
from pylab import cm


from utils import pixels_in_range, pixelise

def build_bivariate_normal_pdf(x_range, y_range, mean=[0.0, 0.0], cov=np.eye(2)):
    x, y = np.meshgrid(np.linspace(-1,1,x_range), np.linspace(-1,1,y_range))
    pos = np.empty(x.shape + (2,))
    pos[:, :, 0] = x
    pos[:, :, 1] = y
    rv = multivariate_normal(mean, cov)
    pdf = rv.pdf(pos)
    return pdf

class CheckerBoard:
    def __init__(self, base_size=15, Nside=32, gaussian=True, mean=None, cov=None):
        self.base_size = base_size
        self.gaussian = gaussian
        self.Nside = Nside
        if self.gaussian:
            assert mean is not None
            assert cov is not None
            self.mean = mean
            self.cov = cov

    def build_base_board(self):
        if self.gaussian:
            n_in_row = 360 // self.base_size
            n_in_col = 180 // self.base_size
            
            pdf = build_bivariate_normal_pdf(self.base_size, self.base_size, self.mean, self.cov)

            row = np.hstack([pdf, -pdf] * (n_in_row // 2))
            base_board = np.vstack([row, -row] * (n_in_col // 2))

            lons = np.linspace(-180,180,base_board.shape[1])
            lats = np.linspace(-90,90,base_board.shape[0])
            lons, lats = np.meshgrid(lons, lats)
            lons, lats = lons.flatten(), lats.flatten()
            base_board = pixelise(base_board.flatten(),self.Nside,lons,lats)
        else:
            lats = np.arange(-90,90+self.base_size,self.base_size)
            lons = np.arange(-180,180+self.base_size,self.base_size)
            base_board = np.ones(hp.nside2npix(self.Nside))
            i = 0
            for lat in lats:
                for lon in lons:  
                    i += 1
                    pixels = pixels_in_range(lat, lon, self.base_size, self.base_size, Nside=self.Nside)
                    if i%2 == 0:
                        base_board[pixels] = -1
        self.board = base_board

    def add_feature(self, feature, lat, lon):
        '''
        Feature is an array of shape (lat_step, lon_step)
        TODO: see if this works with different Nsides
        '''
        lat_step, lon_step = feature.shape
        pixels = pixels_in_range(lat, lon, lat_step, lon_step, Nside=self.Nside)
        self.board[pixels] = feature

if __name__ == '__main__':
    gaussian_board = CheckerBoard(mean=[0.0, 0.0], cov=np.eye(2)*0.15)
    gaussian_board.build_base_board()
    hp.write_map(f'../../data/PhaseVelMaps/gaussian_chkrbrd{gaussian_board.base_size}.fits', gaussian_board.board, overwrite=True)

    hp.mollview(gaussian_board.board,cmap=cm.seismic_r,flip='geo')
    hp.graticule(gaussian_board.base_size)
    plt.savefig(f'../../figs/gaussian_chkrbrd{gaussian_board.base_size}.png')

    chckrbrd = CheckerBoard(gaussian=False)
    chckrbrd.build_base_board()
    hp.write_map(f'../../data/PhaseVelMaps/chkrbrd{chckrbrd.base_size}.fits', chckrbrd.board, overwrite=True)

    hp.mollview(chckrbrd.board,cmap=cm.seismic_r,flip='geo')
    hp.graticule(chckrbrd.base_size)
    plt.savefig(f'../../figs/chkrbrd{chckrbrd.base_size}.png')