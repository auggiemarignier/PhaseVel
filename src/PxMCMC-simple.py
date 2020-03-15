from PxMCMC import Params, PxMCMC
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pys2let
import healpy as hp

# Simple checkerboard example

class Checkerboard:
    def __init__(self,Nside=32,step=30):
        self.lats = np.arange(-90,90+step,step)
        self.longs = np.arange(-180,180+step,step)
        self.Nside = Nside
        self.step = step

        self.fill_board()

    def fill_board(self):
        self.board = np.zeros(hp.nside2npix(self.Nside))
        i=0
        for lat in self.lats:
            for long in self.longs:  
                i += 1
                pixels = []
                for lt in np.arange(lat-self.step/2,lat+self.step/2):
                    if lt < -90:
                        lt = -90
                    if lt > 90:
                        lt =90
                    for ln in np.arange(long-self.step/2,long+self.step/2):
                        if ln < -180:
                            ln += 360
                        if ln > 180:
                            ln -= 360
                        pixels.append(hp.ang2pix(self.Nside,ln,lt,lonlat=True))
                if i%2 == 0:
                    self.board[pixels] = 1

    def plot_board(self):
        hp.mollview(self.board,cmap=cm.jet)
        plt.show()
                


class Simple(PxMCMC):
    '''
    Need to override forward, and gradg to ignore ylmavs
    '''
    def forward(self,X):
        wav_lm,scal_lm = self.expand_mlm(X)
        scal_lm_hp = pys2let.lm2lm_hp(scal_lm,self.params.L+1)
        wav_lm_hp = np.zeros([(self.params.L+1)*(self.params.L+2)//2,self.params.nscales],dtype=np.complex)
        for j in range(self.params.nscales):
            wav_lm_hp[:,j] = pys2let.lm2lm_hp(np.ascontiguousarray(wav_lm[:,j]),self.params.L+1)
        clm_hp = pys2let.synthesis_axisym_lm_wav(wav_lm_hp,scal_lm_hp,self.params.B,self.params.L+1,self.params.J_min)
        clm = np.real(pys2let.lm_hp2lm(clm_hp,self.params.L+1))
        return clm

    def calc_gradg(self,preds):
        diff = (preds-self.data)/self.params.sig_d
        diffs = np.sum(diff,axis=1)
        arrays = [diffs for _ in range(self.basis.shape[1])]
        diffs = np.concatenate(arrays)
        gradg = self.pf*diffs
        return gradg

if __name__ == '__main__':
    board = Checkerboard()
    board.board += np.random.normal(0,0.1,board.board.shape)
    board.plot_board()
    
