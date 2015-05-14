"""
Basic Likelihood tools such as functions for computing likelihoods, Latin Hypercube sampling (efficient random sampling) and other tools which don't belong on other packages, or apply to multiple packages.
"""
from BIP.Bayes import Melding
from BIP.Bayes import like
from BIP.Bayes import PlotMeld

__all__ = ['Melding', 'like', 'PlotMeld']
