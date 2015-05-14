# To change this template, choose Tools | Templates
# and open the template in the editor.
"""
Conjugate prior classes
"""
__author__ = "fccoelho"
__date__ = "$01/12/2009 12:12:39$"
__docformat__ = "restructuredtext en"

class Beta:
    def __init__(self, support, moments):
        """Conjugate Prior Distribution class.
        Provides efficient sampling from prior, and posterior

        :Parameters:
            -`support`: tuple with (lower,upper) limits of the variable distribution support
            -`moments`: list with moments of the prior distribution
        """
        pass
    def set_likelihood(self):
        pass
    def get_posterior(self):
        """
        :Return:
            posterior in the form of a frozen scipy rng
            
        """
        pass

if __name__ == "__main__":
    pass
