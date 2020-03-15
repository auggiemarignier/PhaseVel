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

    def fill_board(self):
        board = np.zeros(hp.nside2npix(self.Nside))
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
                    board[pixels] = 1
        return board

    def plot_board(self,board,min=0,max=1):
        hp.mollview(board,cmap=cm.jet,min=min,max=max)
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
    clean = board.fill_board()
    noisy = clean + np.random.normal(0,0.1,clean.shape)

    vmin = np.min(noisy)
    vmax = np.max(noisy)
    board.plot_board(clean,min=vmin,max=vmax)
    board.plot_board(noisy,min=vmin,max=vmax)

    B = 1.5
    L = 10
    dirs = 1
    spin=0
    J_min = 2
    lmda=3e-5
    delta=1e-5
    mu=1e4
    sig_d=0.6718
    sig_m=1
    nsamples = int(1e6)
    nburn = int(1e3)
    ngap = int(1e2)
    hard = False

    print(' Setting parameters...')
    params = Params(L,B,dirs,spin,J_min,lmda,delta,mu,sig_d,sig_m,nsamples,nburn,ngap,hard)

    algo='MYULA'
    print(f' Setting up {algo}...')
    mcmc = Simple(algo,data,Ylmavs,p_weights,params)
    print(f' Running {algo}...')
