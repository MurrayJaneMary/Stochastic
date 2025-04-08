def draw_bb18(ndraw,lmax=10,z3=True):
    import numpy as np
    """
    return BB18 realizations of shape: [number gauss coeff, number realizations]
    set z3=True for BB18.Z3; z3=False for BB18
    """
    # inner functions
    def get_numgc(lmax=10):
        """
        returns number of gauss coeffs for given max degree
        """
        ngc = lmax**2 + 2*lmax
        return ngc

    def get_g10dist(mu,sig,nn):
        """
        returns vector of g10 gauss coeffs given mu, sigma and number of draws
        """
        g10dist = np.random.normal(mu,sig,nn)
        return g10dist

    def coefnum(l,m,gh):
        """
        returns coefficient number following sequence 0:g10,1:g11,2:h11...X:glm,Y:hlm
        """
        if l<1:
            print('l must be >=1')
            return None
        if m>l: 
            print('m must be <= l')
            return None
        if gh not in ['g','h']: 
            print('gh must be either g or h')
        if m==0 and not gh=='g':
            print('gh must be g if m==0')
            return None
        cntr = 0
        for ll in np.arange(1,l+1):
            for mm in np.arange(0,ll+1):
                cntr += 1
                if ll==l and mm==m and gh=='g':
                    return cntr-1
                if mm>0:
                    cntr += 1
                    if ll==l and mm==m and gh=='h':
                        return cntr-1

    def get_corrlm(lmax=10):
        """
        returns correlation matrix as defined in Bono et al
        positive definite symmetric matrix
        elements follow numbering of 0:g10,1:g11,2:h11...
        """
        # defined constants from Bono et al
        g10g30 = (coefnum(1,0,'g'),coefnum(3,0,'g'),0.51)
        g11g31 = (coefnum(1,1,'g'),coefnum(3,1,'g'),0.55)
        h11h31 = (coefnum(1,1,'h'),coefnum(3,1,'h'),0.53)
        g20g40 = (coefnum(2,0,'g'),coefnum(4,0,'g'),0.14)
        g21g41 = (coefnum(2,1,'g'),coefnum(4,1,'g'),0.60)
        h21h41 = (coefnum(2,1,'h'),coefnum(4,1,'h'),0.58)
        g22g42 = (coefnum(2,2,'g'),coefnum(4,2,'g'),0.42)
        h22h42 = (coefnum(2,2,'h'),coefnum(4,2,'h'),0.37)

        ghterms = [g10g30, g11g31, h11h31, g20g40, g21g41, h21h41, g22g42, h22h42]

        # number of gauss coeffs given max degree
        ngc = get_numgc(lmax)

        # def identity matrix
        corrlm = np.eye(ngc)

        # assign cov terms
        for gh in ghterms:
            corrlm[gh[0],gh[1]] = gh[2]
            corrlm[gh[1],gh[0]] = gh[2]

        return corrlm

    def get_alpha(g10,beta):
        """
        returns alpha (following tk03) following polynomial fit to beta, scaled by g10
        """
        p1 = -0.01042
        p2 = 0.04757 
        p3 = -0.6069

        f1 = lambda b: p1*b**2 + p2*b + p3
        alpha = f1(beta)*-np.abs(g10)
        return alpha

    def get_sigmal(g10,beta,lmax=10):
        """
        returns vector of sigma for each gauss coeff following TK03
        needs median g10 (for alpha) and beta
        element follows sequence of 0:g10, 1:g11, 2:h11, etc up to lmax(default 10)
        """
        ca=0.547 # ratio of Rcmb/Re
        a2 = get_alpha(g10,beta)**2

        s_even = lambda ll: np.sqrt(ca**(2*ll)*a2/((ll+1)*(2*ll+1)))

        sigmal = []
        for ll in np.arange(1,lmax+1):
            for mm in np.arange(0,ll+1):
                if np.mod(ll-mm,2)==0: #even
                    cc = s_even(ll)
                    sigmal.append(cc)
                else: # odd
                    cc = beta*s_even(ll)
                    sigmal.append(cc)
                if mm>0: # h term
                    if np.mod(ll-mm,2)==0: #even
                        cc = s_even(ll)
                        sigmal.append(cc)
                    else: # odd
                        cc = beta*s_even(ll)
                        sigmal.append(cc)

        return sigmal

    def get_covlm(beta,g10,lmax=10):
        """
        returns covariance matrix following Bono et al
        where Cov(i,j) = corr(i,j)*sigma(i)*sigma(j)
        needs max num degrees, beta, median g10, sigma following TK03 and correlation matrix
        """
        ngc = get_numgc(lmax)
        sigmal = get_sigmal(beta=beta,g10=g10,lmax=lmax)
        corrlm = get_corrlm(lmax)

        # convert correlation matrix rho_bar to covariance (scaled by TK03 variance)
        covlm = np.nan*np.zeros([ngc,ngc])
        for ii in range(ngc):
            for jj in range(ngc):
                covlm[ii,jj] = corrlm[ii,jj]*sigmal[ii]*sigmal[jj]

        return covlm

    def get_mus(ndraw,med_g10,s10,g20,g30,lmax=10):
        """
        return vector of mean terms for mvrnd
        all zero-mean except for g10, and g20,g30 for BB18.Z3
        """
        ngc = get_numgc(lmax)
        g10v = get_g10dist(mu=med_g10,sig=s10,nn=ndraw)
        mus = np.zeros([ngc,len(g10v)]) # [number of terms, number of realizations]
        mus[0,:] = g10v # g10 mean term
        mus[3,:] = med_g10*g20 # g20 mean term
        mus[8,:] = med_g10*g30 # g30 mean term

        return mus

    def get_draw(mus,covlm):
        """
        helper function for parallel draws, returns single draw given vector
        of means and covariance matrix
        """
        ght = np.random.multivariate_normal(mus,covlm)
        return ght
    
    # inits

    med_g10 = -22.04
    
    if z3:  # BB18.Z3
        beta = 2.82
        s10 = 10.74
        g20 = -0.65
        g30 = 0.29
    else: # BB18
        beta = 2.82
        s10 = 10.80
        g20 = 0
        g30 = 0    
        
    covlm = get_covlm(lmax=10,beta=beta,g10=med_g10)
    numgc = get_numgc(lmax)
    mus = np.zeros(numgc)
    mus[0] = med_g10
    mus[3] = g20
    mus[8] = g30
    
    
    bb18_draw = np.random.multivariate_normal(mus,covlm,ndraw).T
    # replace g10 with defined mu, sigma
    bb18_g10 = np.random.normal(loc=med_g10,scale=s10,size=[1,ndraw])
    bb18_draw[0,:] = bb18_g10
    
    return bb18_draw
