# -*- coding:utf-8 -*-
from distutils.core import setup
#from distutils.extension import Extension
import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages
#try:
#    from Cython.Distutils import build_ext
#except:
#    print "You don't seem to have Cython installed. Please get a"
#    print "copy from www.cython.org and install it"
#    sys.exit(1)




setup(name='BIP', 
        version  = '0.3.12',
        author = 'Flavio Codeco Coelho', 
        author_email = 'fccoelho@gmail.com', 
        url = 'http://code.google.com/p/bayesian-inference/',
        description = 'Bayesian Inference Tools for Python',
        zip_safe = False,
        packages = find_packages(),#['','BIP','BIP.SDE','BIP.Bayes','BIP.SMC','BIP.Bayes.general','BIP.Bayes.conjugate','BIP.Bayes.tests','BIP.Viz'],
        test_suite = 'nose.collector', 
        license = 'GPL',
#        cmdclass = {'build_ext': build_ext},
        #ext_modules=[Extension('BIP/SDE/gillespie', ['BIP/SDE/gillespie.c'])],
        
      )
