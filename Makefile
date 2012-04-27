default:
	cd pycosmc/derivers/hubble_theta && make
	cd pycosmc/models/camb/camb && make
	cd pycosmc/likelihoods/wmap && make

clean:
	cd pycosmc/derivers/hubble_theta && make clean
	cd pycosmc/models/camb/camb && make clean
	cd pycosmc/likelihoods/wmap && make clean

pycosmc/likelihoods/spt_k11/bandpowers:
	cd pycosmc/likelihoods/spt_k11 && ./get

