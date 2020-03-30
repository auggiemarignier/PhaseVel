import numpy as np
import healpy as hp


def flatten_mlm(wav_lm,scal_lm):
    '''
    Takes a set of wavelet and scaling coefficients and flattens them into a single vector
    '''
    buff = wav_lm.ravel(order='F')
    mlm = np.concatenate((scal_lm,buff))
    return mlm

def expand_mlm(mlm,nscales):
    '''
    Sepatates scaling and wavelet coefficients from a single vector to separate arrays.
    '''
    v_len = mlm.size//(nscales+1)
    assert v_len > 0
    scal_lm = mlm[:v_len]
    wav_lm = np.zeros((v_len,nscales),dtype=np.complex)
    for i in range(nscales):
        wav_lm[:,i] = mlm[(i+1)*v_len:(i+2)*v_len]
    return wav_lm, scal_lm

def soft(X,T=0.1):
    '''
    Soft thresholding of a vector X with threshold T.  If Xi is less than T, then soft(Xi) = 0, otherwise soft(Xi) = Xi-T. 
    '''
    t = np.zeros_like(X)
    for i,x in enumerate(X):
        if abs(x) == 0:
            continue
        t[i] = x*max(abs(x)-T,0)/abs(x)
    return t

def hard(X,T=0.1): 
    '''
    Hard thresholding of a vector X with fraction threshold T. T is the fraction kept, i.e. the largest 100T% absolute values are kept, the others are thresholded to 0.
    TODO: What happens when all elements of X are equal?
    '''
    X_srt = np.sort(abs(X))
    thresh_ind=int(T*len(X))
    thresh_val = X_srt[-thresh_ind]
    X[abs(X)<thresh_val]=0
    return X

def pixels_in_range(lat, lon, lat_range, lon_range, Nside=32):
    pixels = []
    for lt in np.arange(lat-lat_range/2,lat+lat_range/2):
        if lt < -90:
            lt = -90
        if lt > 90:
            lt =90
        for ln in np.arange(lon-lon_range/2,lon+lon_range/2):
            if ln < -180:
                ln += 360
            if ln > 180:
                ln -= 360
            pixels.append(hp.ang2pix(Nside,ln,lt,lonlat=True))
    return pixels

def pixelise(signal,Nside,longs,lats):
    Npix   = hp.nside2npix(Nside)
    pixnum = hp.ang2pix(Nside, longs, lats, lonlat=True)
    amap   = np.zeros(Npix)
    count  = np.zeros(Npix)
    nsample = len(signal)
    for i in range(nsample):
        pix = pixnum[i]
        amap[pix] += signal[i]
        count[pix]+= 1.
    for i in range(Npix):
        if count[i]>0:
            amap[i] = amap[i]/count[i]
        else:
            amap[i] = hp.UNSEEN
    return amap