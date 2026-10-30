# MH-MCMC
This is a demonstration of the Metropolis-Hastings algorithm for Bayesian inference. It uses Supernovae Ia data to infer cosmological parameters, available from [http://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html](http://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html).

## Run the code
A CBLAS implementation (the Makefile links 'OpenBLAS', [https://www.openblas.net](https://www.openblas.net); change `LINK_LIB` to use another, e.g. `-lcblas`) and 'gcc' are a pre-requisite to use this code. To run the code, type `make` and hit Enter. This will compile the code and generate an executable 'mcmc'. We need to provide an initial guess for $\Omega_m$ and $h$ anywhere strictly between $0.20 < \Omega_m < 0.40$ and $0.60 < h < 0.80$. For example, to give an initial guess of $\Omega_m = 0.30$ and $h = 0.67$, the executable can be run as `./mcmc 0.30 0.67`. This generates a file 'walker_chain.txt', containing the random-walker chain for $\Omega_m$ values (first column) and $h$ values (second column).
An example plot ('inference.png') has been provided in the directory 'examples', which shows the inferred parameters with initial guess $\Omega_m = 0.30$ and $h=0.67$.

## How the code works
It reads redshifts and distance modulus ($\mu$) values from the file 'jla_mub.txt' and the corresponding covariance matrix from the file 'jla_mub_covmatrix.dat'. The model for distance modulus is used as $\mu = 25 - 5\log_{10}(h) + 5\log_{10}(D_L(z;\Omega_m, h=1))$ while estimating the log-likelihood, with $D_L(z;\Omega_m, h)$ being the luminosity distance in Mpc. A flat prior over the box above and a flat proposal distribution have been used for the mcmc random-walker.

Note that supernovae alone do not measure $h$: it is degenerate with the absolute magnitude of SNe Ia, which was fixed when the binned JLA distance moduli were computed. The inferred $h \approx 0.70$ therefore reflects that choice, while $\Omega_m$ is a genuine constraint.
