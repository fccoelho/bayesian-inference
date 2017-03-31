Overview
========


The Bip Package is a collection of useful classes for basic Bayesian inference. Currently, its main goal is to be a tool
for learning and exploration of Bayesian probabilistic calculations.

Currently it also includes subpackages for stochastic simulation tools which are not strictly related to Bayesian
inference, but are currently being developed within BIP. One such package is the BIP.SDE which contains a parallelized
solver for stochastic differential equations, an implementation of the Gillespie direct algorithm.

The Subpackage Bayes also offers a tool for parameter estimation of Deterministic and Stochastic Dynamical Models. This
tool will be fully described briefly in a scientific paper currently submitted for publication.

Instalation
-----------

BIP only works on Linux environments. Here  we will include brief instruction on how to install on Ubuntu Linux.
These should be easy to adapt to other distributions.

Binary dependencies:
^^^^^^^^^^^^^^^^^^^

BIP depends on some packages which must be installed by the OS package manager the package itself can be installed with
the pip command.

.. code-block:: bash

   $ sudo apt install gnuplot5-X11 build-essential libgsl-dev python3-pip
   $ sudo pip3 install -U BIP
