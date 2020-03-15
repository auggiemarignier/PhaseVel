import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy
from scipy.stats import norm, multivariate_normal, expon
import pys2let
import healpy as hp
import cmath, math
from multiprocessing import Pool
from tqdm import tqdm


class Params:
    def __init__(self,L,B,dirs,spin,J_min,lmda,delta,mu,sig_d,sig_m,nsamples,nburn,ngap,hard):
        '''
        Sets up parameters including wavelet parameters, tuning parameters and chain length, thinning and burn
        '''
        self.L = L  # maximum angular order
        self.B = B  # wavelet parameter
        self.dirs = dirs  # number of directions for directional wavelets
        self.spin = spin  # spin of the field
        self.J_min = J_min  # minimum wavelet scale
        self.J_max = pys2let.pys2let_j_max(self.B,self.L,self.J_min)  # maximum wavelet scale
        self.nscales=self.J_max-self.J_min+1  # number of wavelet scales
        self.lmda = lmda  # proximity parameter. tuned to make proxf abritrarily close to f
        self.delta = delta  # Forward-Euler approximation step-size
        self.mu = mu  # regularization parameter
        self.sig_d = sig_d  # data errors, could be estimated hierarchically
        self.sig_m = sig_m  # model parameter errors
        self.nsamples = nsamples  # number of desired samples
        self.nburn = nburn  # burn-in size
        self.ngap = ngap  # Thinning parameter=number of iterations between samples. reduces correlations between samples
        self.hard = hard

        assert type(self.L) == int and self.L >= 0
        assert self.B > 1
        assert type(self.dirs) == int and self.dirs > 0
        assert type(self.spin) == int
        assert type(self.J_min) == int and self.J_min >= 0
        assert type(self.nscales) == int and self.nscales >= 0
        assert self.lmda > 0
        assert self.delta > 0
        assert self.mu >= 0
        assert self.sig_d > 0
        assert self.sig_m > 0
        assert type(self.nsamples) == int and self.nsamples > 0
        assert type(self.nburn) == int and self.nburn >= 0
        assert type(self.ngap) == int and self.ngap >=0
        assert type(self.hard) == bool
        

class PxMCMC:
    def __init__(self,algo,data,Ylmavs,p_weights,params):
        '''
        Initialises proximal MCMC algorithm.  Reads in data and path averaged spherical harmonics, which, if necessary, is rearranged to match indexing of s2let. STILL NEED TO CHECK INDEXING OF YLMAVS GENERATED.  Sets up the wavelet basis functions.  Calculates prefactors of the gradg function which are constant throughout the chain.
        '''
        assert algo == 'MYULA' or algo == 'PxMALA'
        self.algo = algo
        self.data = data
        self.npaths = len(data)
        self.p_weights = p_weights
        self.params = params

        assert len(self.p_weights) == self.npaths

        B = self.params.B
        L = self.params.L
        dirs = self.params.dirs
        spin = self.params.spin
        J_min = self.params.J_min

        phi_l, psi_lm = pys2let.wavelet_tiling(B,L+1,dirs,spin,J_min)
        psi_lm = psi_lm[:,J_min:]
        phi_lm = np.zeros(((L+1)**2,1), dtype=np.complex)
        for ell in range(L+1):
            phi_lm[ell*ell+ell] = phi_l[ell]     
        self.basis = np.concatenate((phi_lm,psi_lm),axis=1)   
    
        self.n_lm = self.basis.shape[0]
        self.nb = self.basis.shape[1]
        self.nparams = self.n_lm*self.nb

        # SETUP IN CASE ORIGINAL INDEXING IS DIFFERENT TO S2LET
        self.Ylmavs = np.zeros_like(Ylmavs)
        s2let_indexes = []
        for ell in range(L+1):
            for em in range(-ell,ell+1):
                s2let_indexes.append(ell*ell+ell+em)  #THIS DOESN'T MOVE ANYTHING AROUND.  CHECK ORIGINAL INDEXING
        s2let_indexes = np.array(s2let_indexes)
        for i,ylmavs in enumerate(Ylmavs.T):
            ind = s2let_indexes[i]
            self.Ylmavs.T[ind] = ylmavs

        assert type(self.nparams) == int and self.nparams > 0
        assert self.Ylmavs.shape == (self.npaths,self.n_lm)

        self.calc_prefactors()

    def calc_prefactors(self):
        '''
        Calculates prefactors of gradg which are constant throughout the chain, and so only need to be calculated once at the start.
        '''    
        prefactors = np.zeros(np.prod(self.basis.shape))
        for i,base in enumerate(self.basis.T):
            base_l0s = [base[l**2+l] for l in range(self.params.L)]
            for ell in range(self.params.L):
                prefactors[i*len(base)+ell**2:i*len(base)+(ell+1)**2] = np.sqrt(4*np.pi/(2*ell+1))*np.real(base_l0s[ell]) 
        self.pf = prefactors
        
    def flatten_mlm(self,wav_lm,scal_lm):
        '''
        Takes a set of wavelet and scaling coefficients and flattens them into a single vector
        '''
        buff = wav_lm.ravel(order='F')
        mlm = np.concatenate((scal_lm,buff))
        return mlm

    def expand_mlm(self,mlm):
        '''
        Sepatates scaling and wavelet coefficients from a single vector to separate arrays.
        '''
        v_len = self.nparams//(self.params.nscales+1)
        scal_lm = mlm[:v_len]
        wav_lm = np.zeros((v_len,self.params.nscales),dtype=np.complex)
        for i in range(self.params.nscales):
            wav_lm[:,i] = mlm[(i+1)*v_len:(i+2)*v_len]
        return wav_lm, scal_lm

    def forward(self,X):
        '''
        Forward modelling.  Takes a vector X containing the scaling and wavelet coefficients generated by the chain and predicts path averaged phase velocity.  Possible extension of this is to include finite frequency kernels.
        '''
        wav_lm,scal_lm = self.expand_mlm(X)
        scal_lm_hp = pys2let.lm2lm_hp(scal_lm,self.params.L+1)
        wav_lm_hp = np.zeros([(self.params.L+1)*(self.params.L+2)//2,self.params.nscales],dtype=np.complex)
        for j in range(self.params.nscales):
            wav_lm_hp[:,j] = pys2let.lm2lm_hp(np.ascontiguousarray(wav_lm[:,j]),self.params.L+1)
        clm_hp = pys2let.synthesis_axisym_lm_wav(wav_lm_hp,scal_lm_hp,self.params.B,self.params.L+1,self.params.J_min)
        clm = np.real(pys2let.lm_hp2lm(clm_hp,self.params.L+1)) # check complexity
        preds = np.matmul(clm,self.Ylmavs.T)
        return preds

    def calc_gradg(self,preds):
        '''
        Calculates the gradient of the data fidelity term, which should guide the MCMC search.
        '''
        diff = (self.p_weights**2)*(preds-self.data)/self.params.sig_d
        diffYlmavs = np.sum(diff*self.Ylmavs.T,axis=1)
        arrays = [diffYlmavs for _ in range(self.basis.shape[1])]
        diffYlmavs = np.concatenate(arrays)
        gradg = self.pf*diffYlmavs
        return gradg

    def soft(self,X,T):
        '''
        Soft thresholding of a vector X with threshold T.  If Xi is less than T, then soft(Xi) = 0, otherwise soft(Xi) = Xi-T. 
        '''
        t = np.zeros_like(X)
        for i,x in enumerate(X):
            if abs(x) == 0:
                continue
            t[i] = x*max(abs(x)-T,0)/abs(x)
        return t
    
    def calc_proxf(self,X):
        '''
        Calculates the prox of the sparsity regularisation term.
        '''
        return self.soft(X,self.params.lmda*self.params.mu/2) 
            
    def hard(self,X,T=0.1): 
        '''
        Hard thresholding of a vector X with fraction threshold T. T is the fraction kept, i.e. the largest 100T% values are kept, the others are thresholded to 0.
        '''
        X_srt = np.sort(abs(X))
        thresh_ind=int(T*self.nparams)
        thresh_val = X_srt[-thresh_ind]
        X[abs(X)<thresh_val]=0
        return X

    def chain_step(self,X,proxf,gradg):
        '''
        Takes a step in the chain.
        '''
        w = np.random.randn(self.nparams)
        return (1-self.params.delta/self.params.lmda)*X + (self.params.delta/self.params.lmda)*proxf - self.params.delta*gradg + np.sqrt(2*self.params.delta)*w

    def logpi(self,X,preds):
        '''
        Calculates the log(posterior) of a model X.  Takes in the predictions, preds, of X.
        '''
        return -self.params.mu*sum(abs(X))-sum((self.data-preds)**2)/(2*self.params.sig_d**2)

    def calc_logtransition(self,X1,X2,proxf,gradg):
        '''
        Calculates the transition probability of stepping from model X1 to model X2 i.e. q(X2|X1).  TO BE REWRITTEN
        '''
        gradlogpiX1 = -(1/self.params.lmda)*(X1-proxf) - gradg
        return -(1/2*self.params.delta)*sum((X2-X1-(self.params.delta/2)*gradlogpiX1)**2) # not sure about sum of squares here

    def accept_prob(self,X_curr,curr_preds,X_prop,prop_preds,proxf,gradg):
        '''
        Calculates the acceptance probability of the propsed model X_pop, as a ratio of the transtion probabilities times the ratio of the posteriors.  Strictly speaking, the returned value should be min(0,p)=min(1,e^p) but this makes no difference in the MH acceptance step.
        '''
        logtransXcXp = self.calc_logtransition(X_curr,X_prop,proxf,gradg)
        logtransXpXc = self.calc_logtransition(X_prop,X_curr,proxf,gradg)
        logpiXc = self.logpi(X_curr,curr_preds)
        logpiXp = self.logpi(X_prop,prop_preds)
        p = np.real(logtransXpXc+logpiXp-logtransXcXp+logpiXc)
        assert not np.isnan(p)
        return p

    def accept(self,alpha):
        '''
        Metropolis-Hastings acceptance step.  Accept if the acceptance probability alpha is greater than a random number.
        '''
        u = np.log(np.random.rand())
        if u <= alpha:
            return True
        else:
            return False

    def mcmc(self):
        '''
        Runs MCMC.  At present, logposteriors are becoming more and more negative and converging abnormally quickly.
        '''
        logPi = np.zeros(self.params.nsamples+1)
        preds = np.zeros((self.params.nsamples+1,self.npaths))
        chain = np.zeros((self.params.nsamples+1,self.nparams),dtype=np.complex)
        # X_curr = np.random.normal(0,self.params.sig_m,self.nparams) + np.random.normal(0,self.params.sig_m,self.nparams)*1j
        X_curr = self.calc_proxf(np.random.uniform(-0.5,0.5,self.nparams) + np.random.uniform(-0.1,0.1,self.nparams)*1j)
        if self.params.hard:
            X_curr = self.hard(X_curr)
        curr_preds = self.forward(X_curr)
        i = 0
        while i<self.params.nsamples:
            if i >= self.params.nburn:
                if self.params.ngap==0 or (i-self.params.nburn)%self.params.ngap==0:
                    logPi[i] = self.logpi(X_curr,curr_preds)
                    preds[i] = curr_preds
                    chain[i] = X_curr
            gradg = self.calc_gradg(curr_preds)
            proxf = self.calc_proxf(X_curr)
            X_prop = self.chain_step(X_curr,proxf,gradg)
            if self.params.hard:
                X_prop = self.hard(X_prop)
            prop_preds = self.forward(X_prop)

            if self.algo =='PxMALA':
                alpha = self.accept_prob(X_curr,curr_preds,X_prop,prop_preds,proxf,gradg)
                print(alpha)
                if self.accept(alpha):
                    X_curr = X_prop
                    curr_preds = prop_preds
                    # i += 1
            if self.algo == 'MYULA':
                X_curr = X_prop
                curr_preds = prop_preds
            print(f'{i+1}/{self.params.nsamples} - logposterior: {logPi[i]}')
            i += 1
        self.logPi = logPi
        self.preds = preds
        self.chain = chain

class Outfile:
    def __init__(self,logpost,preds,chain,outpath,binary=True):
        '''
        Initialises required output values, path and whether or not to save as a binary.
        '''
        self.dictnry = {'logposterior':logpost,'predictions':preds,'chain':chain}
        self.outpath = outpath
        self.binary = binary

    def write_outfile(self,key):
        '''
        Writes out an output value identified by key.
        '''
        outfile = f'{self.outpath}/{key}'
        if self.binary:
            np.save(outfile,self.dictnry[key])
            print(f'{key} written to {outfile}.npy')
        else:
            np.savetxt(outfile,self.dictnry[key])
            print(f'{key} written to {outfile}')
        
    def write_outfiles(self):
        'Writes out all the output files.'
        for key in self.dictnry:
            self.write_outfile(key)

if __name__=="__main__":
    print(' Reading inputs...')
    # read_data = np.loadtxt('../inputs/hvh.000S061.asc.rwt3_std4')
    # data = read_data[:,4]
    # p_weights = 1/read_data[:,7]
    data = np.loadtxt('../outputs/000S061_synth')
    p_weights = np.ones(data.shape)
    Ylmavs = (np.loadtxt('../outputs/000S061YlmavL10'))[:,7:]
    print(f'    Data shape: {data.shape}, Ylmavs shape: {Ylmavs.shape}')

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
    mcmc = PxMCMC(algo,data,Ylmavs,p_weights,params)
    print(f' Running {algo}...')
    mcmc = mcmc.mcmc()

    writer = Outfile(mcmc.logPi,mcmc.preds,mcmc.chain,'../outputs')
    writer.write_outfiles()

    '''
    Seems to be falling into a local minimum very quickly and can't get out.  Oddly enough posteriors initially go down then back up to a stable value.  Need to find a way to explore the parameter space better.

    PxMALA never accepts a new model, so never progresses.  Check MH step, in particular transition probability, which still doesn't have a test.  Also not clear on how the transition probability ratio implemented in the MATLAB version is obtained from the definition in the paper.

    Increase area of initial random model search.  Maybe pick a huge ensemble of random points and build a chain off the most likely one...

    Implement parallel chains.  These would hopefully help find the global minimum.

    Uncertainty quantification: pixel credible intervals

    Try finite-frequency kernels

    Decide what to do with complex numbers.  Either just take the real part of harmonic coefficients (currently implemented) or assign 0 probability to models that would result in a complex signal.

    '''