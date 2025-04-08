ERDA 2420: Associated code and demo for BB18 Giant Gaussian Process models, Bono et al.

Included in this ZIP file:
- BB18_demo.ipynb: Jupyter notebook demo showing draw_bb18 function usage
- bb18_draws.csv: CSV file containing 10000 realizations of BB18 up to degree 10, rows are Gauss coeffs sequenced g10, g11, h11, g20 ... g10-10, h10-10, columns are realizations. CSV shape is (120,10000).
- bb18z3_draws.csv: CSV file containing 10000 realizations of BB18.Z3 up to degree 10, rows are Gauss coeffs sequenced g10, g11, h11, g20 ... g10-10, h10-10, columns are realizations. CSV shape is (120,10000).
- draw_bb18.py: Python 3 module containing functions used to generate BB18 GGP models. Requires Numpy.
	Internal functions:
		get_numgc(lmax=10): returns number of gauss coeffs for given max degree 
		get_g10dist(mu,sig,nn): returns vector of g10 gauss coeffs given mu, sigma and number of draws 
		coefnum(l,m,gh): returns coefficient number following sequence 0:g10,1:g11,2:h11...X:glm,Y:hlm 
		get_corrlm(lmax=10): returns correlation matrix as defined in Bono et al, positive definite symmetric matrix elements follow numbering of 0:g10,1:g11,2:h11... 
		get_alpha(g10,beta): returns alpha (following tk03) following polynomial fit to beta, scaled by g10 
		get_sigmal(g10,beta,lmax=10): returns vector of sigma for each gauss coeff following TK03 needs median g10 (for alpha) and beta element follows sequence of 0:g10, 1:g11, 2:h11, etc up to lmax(default 10) 
		get_covlm(beta,g10,lmax=10): returns covariance matrix following Bono et al where Cov(i,j) = corr(i,j)*sigma(i)*sigma(j) needs max num degrees, beta, median g10, sigma following TK03 and correlation matrix 
		get_mus(ndraw,med_g10,s10,g20,g30,lmax=10): return vector of mean terms for mvrnd, all zero-mean except for g10, and g20,g30 for BB18.Z3 
		get_draw(mus,covlm): helper function for parallel draws, returns single draw given vector of means and covariance matrix 
