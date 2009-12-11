# To change this template, choose Tools | Templates
# and open the template in the editor.
"""
Module implementing MCMC samplers 

"""

__author__="fccoelho"
__date__ ="$09/12/2009 10:44:11$"
__docformat__ = "restructuredtext en"

#TODO: Finish implementation and test
class Metropolis(object):
    """
    Metropolis Hastings sampler class
    """
    def __init__(self, proposal_dist,likfun):
        self.salt_band = 0.05
        pass
    def propose(self):
        """generates proposals"""
        pvar = self.proposal_variance*self.adaptscalefactor
        theta = [self.theta_dists[dist]() for dist in self.q1theta.dtype.names]
        prop = self.model(*theta)
        if t == 1:
            prop=[(v,) for v in prop]
        return theta,prop
    def step(self,n=1):
        """
        Does the actual sampling loop for 
        n steps
        """
        ptheta = recarray((self.K),formats=['f8']*self.ntheta, names = self.post_theta.dtype.names)
        i=0;j=0;rej=0 #total samples,accepted samples, rejected proposals
        last_lik = None
        liklist = []
        while j < n:
            theta,prop = self.propose()
            #calculate likelihoods
            lik=0 
            for k in xrange(self.nphi): #loop on series
                if self.q2phi.dtype.names[k] not in data:
                    continue#Only calculate liks of series for which we have data
                obs = data[self.q2phi.dtype.names[k]]
                for p in xrange(t): #Loop on time
                    lik += likfun(obs[p],prop[p][k],1./variance)
            #store samples
            if last_lik == None: #on first sample
                last_lik = lik
                continue
            
            if not (lik-last_lik) > numpy.log(1):
                if not numpy.log(random())< lik-last_lik:
                    rej +=1 #adjust rejection counter
                    i +=1
                    #print theta
                    continue
            if i >= burnin:#store only after burnin samples
                self.phi[j] = prop[0] if t==1 else [tuple(point) for point in prop]
                ptheta[j] = tuple(theta)
                liklist.append(lik)
                j += 1 #update good sample counter 
            last_lik = lik
            i+=1
            
        self.post_theta = self.imp_sample(self.L,ptheta,liklist)
        ar = (i-rej)/float(i) 
        if ar > 0.9:
            if self.salt_band < 0:
                self.salt_band = 0.05
            elif self.salt_band < 30:
                self.salt_band *= 10
        if ar < 0.4:
            self.salt_band = 0.001 #no more expansion
        print "Total steps(i): ",i,"rej:",rej, "j:",j
        print ">>> Salt-Band: %s"%self.salt_band
        print ">>> Acceptance rate: %s"%ar
        
        return 1
    def tune(self, ar):
        """
        Tune the proposal distribtion variance
        in the case of using a Normal proposal distribution 
        """
        if self.proposal_dist == "prior":
            return
        if ar<0.05:
            self.adaptscalefactor *= 0.5
        elif ar <0.2:
            self.adaptscalefactor *= 0.9
        elif ar >0.9:
            self.adaptscalefactor *= 10
        elif ar >0.75:
            self.adaptscalefactor *= 2
        elif ar >0.5:
            self.adaptscalefactor *= 1.1
    def add_salt(self,dataset,band):
        """
        Adds a few extra uniformly distributed data 
        points beyond the dataset range.
        This is done by adding from a uniform dist.
        
        :Parameters:
            -`dataset`: vector of data
            -`band`: Fraction of range to extend: [0,1[
            
        :Returns:
            Salted dataset.
        """
        dmax = max(dataset)
        dmin = min(dataset)
        drange = dmax-dmin
        hb = drange*band/2.
        d = numpy.concatenate((dataset,stats.uniform(dmin-hb,dmax-dmin+hb).rvs(self.K*.05)))
        return d

if __name__ == "__main__":
    pass
