from __future__ import absolute_import
from __future__ import print_function

from setuptools import setup, Extension
from setuptools import find_packages
from Cython.Build import cythonize
from Cython.Distutils import Extension
from Cython.Distutils import build_ext

extensions = [Extension('BIP.SDE.cgillespie', ["BIP/SDE/cgillespie.pyx"]), ]
print(cythonize(extensions))

setup(name='BIP',
      version='0.6.12',
      author='Flavio Codeco Coelho',
      author_email='fccoelho@gmail.com',
      url='http://code.google.com/p/bayesian-inference/',
      description='Bayesian Inference Tools for Python',
      zip_safe=False,
      packages=find_packages(),
      # ['','BIP','BIP.SDE','BIP.Bayes','BIP.SMC','BIP.Bayes.general','BIP.Bayes.conjugate','BIP.Bayes.tests','BIP.Viz'],
      package_data={'': ['*.txt']},
      ext_modules=cythonize(extensions),
      setup_requires=["cython", "numpy", "cythongsl"],
      install_requires=["numpy", "scipy", "openopt", "liveplots", "cython"],
      test_suite='nose.collector',
      license='GPLv3',
      cmdclass={'build_ext': build_ext},
      # ext_modules=[Extension('BIP/SDE/gillespie', ['BIP/SDE/gillespie.c'])],
      )

import os

# This to avoid creating the log file with super-user privileges during instalation.
try:
    os.unlink('/tmp/BIP.log')
except:
    pass
