hubble_theta.so: hubble_theta_convert.f90 hubble_theta_wrap.f90
	gfortran -fPIC -O3 -c hubble_theta_convert.f90
	f2py hubble_theta_convert.o -c hubble_theta_wrap.f90 -m hubble_theta

clean:
	-rm -rf *.o *.so *.mod
