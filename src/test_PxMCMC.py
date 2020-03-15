from PxMCMC import *
import unittest
import pytest

npaths = 10
L = 4
n_lm = (L+1)**2
B = 1.5
dirs = 1
spin=0
J_min = 2
J_max = pys2let.pys2let_j_max(B,L,J_min)
nscales = J_max-J_min+1
lmda=1e-1
delta=lmda/4
mu=lmda
sig_d=0.6
sig_m=0.1
nsamples = 10
nburn = 2
ngap = 1
nparams = (nscales+1)*n_lm
hard = False

t_data = np.random.randn(npaths)
p_weights = np.random.randint(low=1,high=10,size=npaths)
t_Ylmavs = np.random.randn(npaths,n_lm)


def test_Params_init():
    Params(L,B,dirs,spin,J_min,lmda,delta,mu,sig_d,sig_m,nsamples,nburn,ngap,hard)
    with pytest.raises(AssertionError):
        Params(L,-B,dirs,spin,J_min,lmda,delta,mu,sig_d,sig_m,nsamples,nburn,ngap,hard)
        Params(L,B,0,spin,J_min,lmda,delta,mu,sig_d,sig_m,nsamples,nburn,ngap,hard)
        Params(L,B,dirs,-2,J_min,lmda,delta,mu,sig_d,sig_m,nsamples,nburn,ngap,hard)
        Params(L,-B,dirs,spin,3.4,lmda,delta,mu,sig_d,sig_m,nsamples,nburn,ngap,hard)
        Params(L,-B,dirs,spin,J_min,lmda,delta,mu,sig_d,sig_m,0,nburn,ngap,hard)
        Params(L,-B,dirs,spin,J_min,lmda,delta,mu,sig_d,sig_m,nsamples,nburn,-1,hard)

def test_PxMCMC_init():
    params = Params(L,B,dirs,spin,J_min,lmda,delta,mu,sig_d,sig_m,nsamples,nburn,ngap,hard)
    mcmc = PxMCMC('MYULA',t_data,t_Ylmavs,p_weights,params)
    mcmc = PxMCMC('PxMALA',t_data,t_Ylmavs,p_weights,params)
    with pytest.raises(AssertionError):
        mcmc = PxMCMC('blah',t_data,t_Ylmavs,p_weights,params)
        mcmc = PxMCMC('MYULA',t_data,t_data,p_weights,params)

def test_flatten_expand():
    params = Params(L,B,dirs,spin,J_min,lmda,delta,mu,sig_d,sig_m,nsamples,nburn,ngap,hard)
    mcmc = PxMCMC('MYULA',t_data,t_Ylmavs,p_weights,params)
    mlm_flat = np.random.randn(nparams) + np.random.randn(nparams)*1j  
    mlm_wav,mlm_scal = mcmc.expand_mlm(mlm_flat)
    assert (mcmc.flatten_mlm(mlm_wav,mlm_scal) == mlm_flat).all()

def test_proxf():
    params = Params(L,B,dirs,spin,J_min,lmda,delta,mu,sig_d,sig_m,nsamples,nburn,ngap,hard)
    mcmc = PxMCMC('MYULA',t_data,t_Ylmavs,p_weights,params)
    mcmc.params.lmda = 1
    mcmc.params.mu = 2
    X = np.array([-2-2j,-1,1j,0,2,3+1j,3])
    expected = np.array([-2/np.sqrt(8)*(1+1j)*(np.sqrt(8)-1),0,0,0,1,(3+1j)*(np.sqrt(10)-1)/np.sqrt(10),2])
    proxf = mcmc.calc_proxf(X)
    assert np.allclose(expected,proxf)

def test_gradg():
    params = Params(L,B,dirs,spin,J_min,lmda,delta,mu,1,sig_m,nsamples,nburn,ngap,hard)
    Ylmavs = np.ones(t_Ylmavs.shape)
    mcmc = PxMCMC('MYULA',t_data,Ylmavs,p_weights,params)
    mcmc.basis = np.ones(mcmc.basis.shape)
    mcmc.calc_prefactors()
    mcmc.p_weights = np.ones(npaths)
    preds = t_data + 1
    gradg = mcmc.calc_gradg(preds)
    expected = len(t_data)*mcmc.pf
    assert (gradg == expected).all()


def test_forward():
    '''
    Got equal the same results with matmul than the naive looping method.
    Take results from linear inversion and see if predictions match data well.
    NOTE: Failing this test could just mean the linear model is shit
    '''
    L=10
    read_data = np.loadtxt('../inputs/hvh.000S061.asc.rwt3_std4')
    data = read_data[:,4]
    p_weights = 1/read_data[:,-1]
    Ylmavs = (np.loadtxt('../outputs/000S061YlmavL10'))[:,7:]
    
    phasevel_lin = hp.read_map('../../data/PhaseVelMaps/hvh.000S061.fits',verbose=False)
    p_lm = hp.map2alm(phasevel_lin,lmax=L+1)
    f_wav_lm_hp, f_scal_lm_hp = pys2let.analysis_axisym_lm_wav(p_lm,B,L+1,J_min)
    f_scal_lm = pys2let.lm_hp2lm(f_scal_lm_hp,L+1)
    f_wav_lm = np.zeros((121,5),dtype=np.complex)
    for i in range(5):
        f_wav_lm[:,i] = pys2let.lm_hp2lm(np.ascontiguousarray(f_wav_lm_hp[:,i]),L+1)

    params = Params(L,B,dirs,spin,J_min,lmda,delta,mu,sig_d,sig_m,nsamples,nburn,ngap,False)
    mcmc = PxMCMC('MYULA',data,Ylmavs,p_weights,params)
    X = mcmc.flatten_mlm(f_wav_lm,f_scal_lm)
    preds = mcmc.forward(X)
    np.savetxt('../outputs/000S061_synth',preds)
    def l2(a,b=None): # l2 difference between two vectors
        if b is not None:
            return np.sqrt(sum((a-b)*(a-b)))
        else:
            return np.sqrt(sum(a*a))
    l2_diff = 100*l2(data,preds)/l2(data)
    assert l2_diff < 10

def test_logpi():
    params = Params(L,B,dirs,spin,J_min,lmda,delta,1,np.sqrt(0.5),sig_m,nsamples,nburn,ngap,False)
    mcmc = PxMCMC('MYULA',t_data,t_Ylmavs,p_weights,params)
    X = np.ones(mcmc.nparams)
    preds = mcmc.data - 1
    mcmc.p_weights = np.ones(npaths)
    expected = - mcmc.nparams - mcmc.npaths
    assert mcmc.logpi(X,preds) == expected

def test_myula_mcmc():
    params = Params(L,B,dirs,spin,J_min,lmda,delta,mu,sig_d,sig_m,nsamples,nburn,ngap,hard)
    mcmc = PxMCMC('MYULA',t_data,t_Ylmavs,p_weights,params)
    mcmc.mcmc()

def test_pxmala_mcmc():
    params = Params(L,B,dirs,spin,J_min,lmda,delta,mu,sig_d,sig_m,nsamples,nburn,ngap,hard)
    mcmc = PxMCMC('PxMALA',t_data,t_Ylmavs,p_weights,params)
    mcmc.mcmc()

def test_outfiles(tmpdir):
    logpost = np.random.randn(5)
    preds = np.random.randn(5,15)
    chain = np.random.randn(5,10)
    dictnry = {'logposterior':logpost,'predictions':preds,'chain':chain}
    out = Outfile(logpost,preds,chain,tmpdir)
    out.write_outfiles()
    for file in dictnry:
        print(dictnry[file])
        print(np.load(f'{tmpdir}/{file}.npy'))
        assert (dictnry[file] == np.load(f'{tmpdir}/{file}.npy')).all()
    out = Outfile(logpost,preds,chain,tmpdir,False)
    out.write_outfiles()
    for file in dictnry:
        assert (dictnry[file] == np.loadtxt(f'{tmpdir}/{file}')).all()

