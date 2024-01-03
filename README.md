# MH-MCMC
This is a demonstration of the Metropolis-Hastings algorithm for Bayesian inference. It uses Supernovae Ia data to infer cosmological parameters, available from [http://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html](http://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html).

## Run the code
'GSL' ([https://www.gnu.org/software/gsl](https://www.gnu.org/software/gsl)) and 'gcc' are a pre-requisite to use this code. To run the code, type 'Make' and hit Enter. This will compile the code and generate an executable 'mcmc'. We need to provide an initial guess for $\Omega_m$ and $h$ anywhere between $0.20 \leq \Omega_m \leq 0.40$ and $0.60 \leq h \leq 0.80$. For example, to give an initial guess of $\Omega_m = 0.30$ and $h = 0.67$, the executable can be run as './mcmc 0.30 0.67'. This generates a file 'walker_chain.txt', containing the random-walker chain for $\Omega_m$ values (first column) and $h$ values (second column).
An example plot ('inference.png') have been provided in the directory 'examples', which shows the infered parameters with initial guess $\Omega_m = 0.30$ and $h=0.67$.

## How the code works
It reads redshifts and distance modulus ($\mu$) values from the file 'jla_mub.txt' and the corresponding covariance matrix from the file 'jla_mub_covmatrix.dat'. The model for distance modulus is used as $\mu = 25 - 5\log(h) + 5\log(D_L(h=1,\Omega_m, z))$ while estimating the log-likelihood, with $D_L(h, \Omega_m, z)$ being the luminosity distance. A flat proposal distribution has been used for the mcmc random-walker.
