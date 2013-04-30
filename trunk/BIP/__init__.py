"""
Bayesian Inference Package containing usefull classes and functions
for doing inference in various applications.

This is not a "Bayesian Statistics" package, i.e., it was not conceived  
to provide data analysis methods such as the one you can find 
on a statistical package such as R (for example).

Here, you will find some the basic building blocks those 
sophisticated Bayesian regression methods are built from, such as 
likelihood functions, MCMC samplers, SMC samplers, etc..

This package exists because such basic tools are not readily accessible 
in Task-oriented statistical software. From these tools, as this package matures,
you will be able to easily build a solution for your own inferential inquiries, 
a solution which may not be available on standard statistical packages.
"""

__docformat__ = "restructuredtext en"

import logging
from logging.handlers import RotatingFileHandler


logger = logging.getLogger("BIP")
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = RotatingFileHandler("/tmp/BIP.log", maxBytes=500000, backupCount=2)
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter and add it to the handlers
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch.setFormatter(formatter)
fh.setFormatter(formatter)
# add the handlers to logger
logger.addHandler(ch)
logger.addHandler(fh)
