import os

if not os.path.exists('spt_multif_0809'):
    os.system("wget http://lambda.gsfc.nasa.gov/data/suborbital/SPT/spectra_2011/bandpowers_spt20082009.tar.gz")
    os.system("tar zxvf bandpowers_spt20082009.tar.gz")

