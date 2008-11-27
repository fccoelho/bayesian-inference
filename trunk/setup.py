# -*- coding:utf-8 -*-
import ez_setup
ez_setup.use_setuptools()
from setuptools import setup

setup(name='BIP', 
        version  = '0.1a',
        author = 'Flavio Codeco Coelho', 
        author_email = 'fccoelho@gmail.com', 
        url = 'http://code.google.com/p/bayesian-inference/',
        description = 'Bayesian Inference Tools for Python', 
        test_suite = 'nose.collector', 
        license = 'GPL', 
        
      )
