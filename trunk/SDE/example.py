"""
Example of an SEIR model with two Infectious classes: subclinical(Is) and clinical(Ic)
        Is
       /  \
S -> E     R
       \  /
        Ic
        
States:
S: Susceptible
E: Exposed
Is: Infectious subclinical
Ic: Infectious clinical
R: Recovered

Transition rates:
b,ks,kc,rs,rc = (0.001, 0.1, 0.1, 0.01, .01)
Transitions:
S -> E : b*S*(Is+Ic)
E -> Is : ks*E
E -> Ic : kc*E
Is -> R : rs*Is
Ic -> R : rc*Ic

"""
from gillespie import Model
import time
from numpy import array
vars = ['S','E','Is','Ic','R']
#rates: b,ks,kc,rs,rc

r = (0.001, 0.1, 0.1, 0.01, .01)
ini = (490,0,10,0,0)
prop = (lambda r,ini:r[0]*ini[0]*(ini[2]+ini[3]),
        lambda r,ini:r[1]*ini[1],
        lambda r,ini:r[2]*ini[1],
        lambda r,ini:r[3]*ini[2],
        lambda r,ini:r[4]*ini[3]
        )

tmat = array([[-1,0,0,0,0],
            [1,-1,-1,0,0],
            [0,1,0,-1,0],
            [0,0,1,0,-1],
            [0,0,0,1,1]
            ])
#for e in prop:
#    print e()
M=Model(vnames=vars,rates = r,inits=ini,tmat=tmat,propensity=prop)
t0 = time.time()
M.run(tmax=80,reps=100)
print 'total time: ',time.time()-t0
t,series,steps = M.getStats()
print steps,'steps'
from pylab import plot , show, legend, errorbar
plot(t,series.mean(axis=2),'-o')
legend(vars,loc=0)
show()
