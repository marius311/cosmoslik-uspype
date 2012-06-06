default: pycosmc/likelihoods/spt_k11/bandpowers
	cd pycosmc/derivers/hubble_theta && make
	cd pycosmc/models/camb/camb && make
	cd pycosmc/likelihoods/wmap && make
	cd pycosmc/likelihood/xfaster && make

clean:
	cd pycosmc/derivers/hubble_theta && make clean
	cd pycosmc/models/camb/camb && make clean
	cd pycosmc/likelihoods/wmap && make clean
	cd pycosmc/likelihood/xfaster && clean

pycosmc/likelihoods/spt_k11/bandpowers:
	cd pycosmc/likelihoods/spt_k11 && ./get

