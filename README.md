# MH-MCMC
This is a demonstration of the Metropolis-Hastings algorithm for Bayesian inference. It uses Supernovae Ia data to infer cosmological parameters, available from [http://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html](http://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html).

## Run the code
To run the code, type 'Make' and hit Enter. This will compile the code and generate an executable 'mcmc'. We need to provide an initial guess for $\Omega_m$ and $h$ anywhere between $0.20 \lesseq \Omega_m \lesseq 0.40$ and $0.60 \lesseq h \lesseq 0.80$. For example, to give an initial guess of $\Omega_m = 0.27$ and $h = 0.67$, the executalbe can be run as './mcmc 0.30 0.67'. This generates a file 'walker_chain.txt', containing the random-walker chain for $\Omega_m$ values (first column) and $h$ values (second column).
