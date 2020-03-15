import pytest
from utils import *

def test_flatten_expand():
    nparams = 100
    nscales = 3
    mlm_flat = np.random.randn(nparams) + np.random.randn(nparams)*1j  
    mlm_wav,mlm_scal = expand_mlm(mlm_flat,nscales)
    assert (flatten_mlm(mlm_wav,mlm_scal) == mlm_flat).all()

def test_soft():
    X = np.array([-2-2j,-1,1j,0,2,3+1j,3])
    expected = np.array([-2/np.sqrt(8)*(1+1j)*(np.sqrt(8)-1),0,0,0,1,(3+1j)*(np.sqrt(10)-1)/np.sqrt(10),2])
    T = 1
    assert np.allclose(expected,soft(X,T))

def test_hard():
    '''
    TODO: test scenarios such as all equal values
    '''
    X = np.arange(-10,10,1)
    expected = np.zeros(20)
    expected[0] = -10
    expected[1] = -9
    expected[-1] = 9
    assert np.allclose(expected,hard(X))