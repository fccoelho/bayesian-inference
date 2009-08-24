# -*- coding:utf-8 -*-
#-----------------------------------------------------------------------------
# Name:        condensation.py
# Project:  Bayesian-Inference
# Purpose:     implementation of the Condensation algorithm, originally conceived by Michael Isard
#                   http://www.cs.ubc.ca/~nando/smc/index.htm
#
# Author:      Flávio Codeço Coelho<fccoelho@gmail.com>
#
# Created:     2008-11-26
# Copyright:   (c) 2008 by the Author
# Licence:     GPL
#-----------------------------------------------------------------------------
__docformat__ = "restructuredtext en"
'''
implementation of the Condensation algorithm, originally conceived by Michael Isard
http://www.cs.ubc.ca/~nando/smc/index.htm
'''

from numpy.random import uniform, normal,  random
from numpy import  array, ndarray,  zeros, nan_to_num
from math import pi, exp, sqrt
import pylab as P

class Model(object):
    """
    Example Model Specification
    """
    def __init__(self):
        self.gdata = GlobalData((0, .2), (-.1, .4, .075), (-.1, .4, .075, .03), .03, 1000, 100)
        self.data = IterationData()
        self.out= zeros((100, 3), dtype = float)
    def setupPriorConditions(self):
        for n in xrange(self.gdata.nsamples):
            self.data.oldpos[n] = self.gdata.PriorModel[0] +self.gdata.PriorModel[1]*normal()
            # The probabilities are not normalised. 
            self.data.cumul_prob_array[n] =  float(n)
            self.data.sample_weights[n] = 1.0
        # The probabilities are not normalised, so store the largest value here 
        # (for simplicity of the binary search algorithm,
        # cumul_prob_array[0] = 0). This can then be used as a
        # multiplicative normalisation constant for the sample_weights
        # array, as well. 
        self.data.largest_cumulative_prob = float(n)
        # This is the initial positions the simulated object. 
        self.data.meas[0] = 0.0;
    def iterate(self, previous, process):
        '''
        The process model for a first-order auto-regressive process is:

        x_{t+1} - mean = (x_t - mean)*scaling + sigma*w_t

        where w_t is unit iid Gaussian noise.
        
        :Parameters:
            -  `previous`: previous data
            -  `process`: processmodel parameter tuple
        '''
        return process[0]+((previous-process[0])*process[1])+process[2]*normal()
        
    def predictSamplePosition(self, new_sample, old_sample):
        '''
        This routine samples from the distribution

        p(x_t | x_{t-1} = oldpos[old_sample])

        and stores the result in new_positions[new_sample]. This is
        straightforward for the simple first-order auto-regressive process
        model used here, but any model could be substituted.
        '''
        self.data.newpos[new_sample] = self.iterate(self.data.oldpos[old_sample], self.gdata.ProcessModel)
        
    def evaluateObservationDensity(self, new_sample):
        '''
        This routine evaluates the observation density

        p(z_t|x_t = newpos[new_sample])

        The observation model in this implementation is a simple mixture of
        Gaussians, where each simulated object is observed as a 1d position
        and measurement noise is represented as Gaussian. For a
        visual-tracking application, this routine would go and evaluate the
        likelihood that the object is present in the image at the position
        encoded by new_positions[new_sample]. 
        '''
        return evaluate_gaussian(self.data.newpos[new_sample]-self.data.meas[1], self.gdata.ObservationModel )
    def obtainObservations(self):
        '''
        In a real implementation, this routine would go and actually make
        measurements and store them in the data.meas structure. This
        simulation consists of an object moving around obeying a
        first-order auto-regressive process, and being observed with its
        true positions coorrupted by Gaussian measurement noise.
        Accordingly, this routine calculates the new simulated true and
        measured position of the object.
        '''
        self.data.meas[0] = self.iterate(self.data.meas[0], self.gdata.SceneModel)
        self.data.meas[1] = self.data.meas[0]+self.gdata.SceneModel[3]*normal()
    def display(self, iteration):
        aggregate = 0.
        aggregate =  sum(self.data.newpos*self.data.sample_weights)/ self.data.largest_cumulative_prob
        
        self.out[iteration, :] = (self.data.meas[1], self.data.meas[0], aggregate)
        #print "%04d: Measured pos. % 3.4lf True pos. % 3.4lf Est. position % 3.4lf\n"%(iteration, self.data.meas[1], self.data.meas[0], aggregate)
        #print "==>Error: ",  self.data.meas[0]-aggregate
            
class GlobalData:
    def __init__(self,  prior, process, scene, observ, nsam, nit):
        '''
        Class to hold global data for the simulation
        
        :Parameters:
            -  `prior`: parameter tuple specifying the model of the prior distribution for the first step.
            -  `process`: the parameters specifying the process model. (mean,scaling, sigma)
            -  `scene`: parameters for the simlation model used to track the process (mean,scaling, sigma, sigma)
            -  `observ`: sigma of the observation model
            -  `nsam`: number of samples
            -  `nit`: number  of terations of the model
        '''
        #The prior distribution of the state is taken to be Gaussian with the parameters stored in this structure.
        self.PriorModel = prior
        #The parameters specifying the process model.
        self.ProcessModel = process
        self.SceneModel = scene
        self.ObservationModel = observ        
        self.nsamples = nsam
        self.niterations = nit
        
class IterationData(object):
    def __init__(self):
        self.newpos = 0.
        self.oldpos = 0.
        self.sample_weights = None
        self.cumul_prob_array = None
        self.largest_cumulative_prob = 1.
        self.meas = [0., 0.] #(true,observed)
        
class Condensation(object):
    def __init__(self,  model):
        self.globaldata = model.gdata
        self.iterdata = model.data
        self.model = model
        self.iterdata.newpos = zeros((self.globaldata.nsamples, 1), dtype = float)
        self.iterdata.oldpos = zeros((self.globaldata.nsamples, 1), dtype = float)
        self.iterdata.sample_weights = zeros((self.globaldata.nsamples, 1), dtype = float)
        self.iterdata.cumul_prob_array =  zeros((self.globaldata.nsamples, 1), dtype = float)
        
        self.model.setupPriorConditions()
    
    def pickBaseSample(self):
        '''
        This is binary search using cumulative probabilities to pick a base
        sample. The use of this routine makes Condensation O(NlogN) where N
        is the number of samples. It is probably better to pick base
        samples deterministically, since then the algorithm is O(N) and
        probably marginally more efficient, but this routine is kept here
        for conceptual simplicity and because it maps better to the
        published literature.
        '''
        choice = random()*self.iterdata.largest_cumulative_prob
        low = 0
        high = self.globaldata.nsamples
        
        while high>(low+1):
            middle = (high+low)/2
            if choice > self.iterdata.cumul_prob_array[middle]:
                low = middle
            else:
                high = middle
                
        return low
        
    def predictNewBases(self):
        '''
        This method computes all of the new (unweighted) sample
        positions. For each sample, first a base is chosen, then the new
        sample position is computed by sampling from the prediction density
        p(x_t|x_t-1 = base). predict_sample_position is obviously
        model-dependent and is found in Model, but it can be
        replaced by any process model required.
        '''
        for n in xrange(self.globaldata.nsamples):
            base = self.pickBaseSample()
            self.model.predictSamplePosition(n, base)
    
    def runFilter(self):
        '''
        '''
        for i in xrange(self.globaldata.niterations):
            self.model.obtainObservations()#Go make necessary measurements
            self.predictNewBases()#Push previous state through process model
            self.calculateBaseWeights() #Apply Bayesian measurement weighting
            self.updateAfterIterating(i)#
        #print self.model.out.shape
        P.plot(nan_to_num(self.model.out[:, 0]), '-^')
        P.plot(nan_to_num(self.model.out[:, 1]), '-^')
        P.plot(nan_to_num(self.model.out[:, 2]), '-^')
        P.legend(['Measured', 'True', 'Estimated'])
        P.show()
    
    def calculateBaseWeights(self):
        '''
        Once all the unweighted sample positions have been computed using
        predict_new_bases, this routine computes the weights by evaluating
        the observation density at each of the positions. Cumulative
        probabilities are also computed at the same time, to permit an
        efficient implementation of pick_base_sample using binary
        search. evaluate_observation_density is obviously model-dependent
        and is found in the Model class, but it can be replaced by any
        observation model required.
        '''
        cumul_total = 0.0
        for n in xrange(self.globaldata.nsamples):
            self.iterdata.sample_weights[n] = self.model.evaluateObservationDensity(n)
            self.iterdata.cumul_prob_array[n] = cumul_total;
            cumul_total += self.iterdata.sample_weights[n];
        self.iterdata.largest_cumulative_prob = cumul_total

    def updateAfterIterating(self, iteration):
        '''
        Go and output the estimate for this iteration (which is a
        model-dependent routine found in Model) and then swap
        over the arrays ready for the next iteration.
        '''
        self.model.display(iteration)
        temp = self.iterdata.newpos
        self.iterdata.newpos = self.iterdata.oldpos
        self.iterdata.oldpos = temp
        
def evaluate_gaussian(val, sigma):
    return 1.0/(sqrt(2.0*pi) * sigma) * exp(-0.5 * (val*val / (sigma*sigma)));


    
    
normal
    
if __name__=="__main__":
    M = Model()
    C = Condensation(M)
    C.runFilter()
    
