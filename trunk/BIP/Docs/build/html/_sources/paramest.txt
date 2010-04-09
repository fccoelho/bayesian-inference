**************************************
Parameter Estimation in Dynamic Models
**************************************


A growing theme in mathematical modeling is uncertainty analysis. The Melding Module provides a Bayesian framework to analyze uncertainty in mathematical models. It includes tools that allow modellers to integrate Prior information about the model's parameters and variables into the model, in order to explore the full uncertainty associated with a model.

This framework is inspired on the original Bayesian Melding paper by Poole and Raftery [1]_, but extended to handle dynamical systems and different posterior sampling mechanisms, i.e., the user has the choice to use Sampling Importance resampling, Approximate Bayesian computations or MCMC.

Once a model is thus parameterized, we can simulate the model, with full uncertainty representation and also fit the model to available data to reduce that uncertaity. Markov chain Monte Carlo algorithms are at the core of the framework, which requires a large number of simulations of the models in order to explore parameter space.

Single Session Retrospective estimation
=======================================

Frequently, we have a complete time series corresponding to one or more state variables of our dynamic model. In such cases it may be interesting to use this information, to estimate the parameter values which maximize the fit of our model to the data. Below are examples of such inference situations.

Example Usage
-------------

This first example includes a simple ODE (an SIR epidemic model) model which is fitted against simulated data to which noise is added:

.. literalinclude:: ../../examples/Bayes/deterministic.py


The code above starts by defining the models parameters and initial conditions, and a function which takes in the parameters runs the model and returns the output.

After that, we Instantiate our fitting Object::

    F = FitModel(300,model,y0,tf,['beta'],['S','I','R'],
                wl=36,nw=1,verbose=False,burnin=100)

Here we have to pass a few arguments: the first (``K=300``) is the number of samples we will take from the joint prior distribution of the parameters to run the inference. The second one (``model``) is the callable(function) which corresponds to the model you want to fit to data. Then you have the initial condition vector(``inits=y0``), the list of parameter names (``thetanames = ['beta']``), the list of variable names (``phinames=['S','I','R']``), inference window length (``wl=36``), number of juxtaposed windows (``nw=1``), verbosity flag (``verbose=False``) and finally the number of burnin samples (``burnin=1000``), which is only needed for if the inference method chosen is ``MCMC``.

One should always have ``verbose=True`` on a first fitting run of a model or if the simulations seems to be taking longer than expected. When ``verbose`` is true, printed and graphical is generated regarding the behavior of fitting, which can be useful to fine tune its parameters.

The next line of code also carries a lot of relevant information about the inference: the specification of the prior distributions. By now you must have noticed that not all parameters included in the model need to be included in the analysis. any number of them except for one can be set constant, which is what happens with the parameter ``tau`` in this example::

    F.set_priors(tdists=[st.norm],tpars=[(1.1,.2)],tlims=[(0.5,1.5)],
        pdists=[st.uniform]*3,ppars=[(0,.1),(0,.1),(.8,.2)],plims=[(0,1)]*3)

here we set the prior distributions for the theta (the model's parameters) and phi (the model's variables). ``tdists``, ``tpars`` and ``tlims`` are theta's distributions, parameters, and ranges. For example here we use a Normal distribution (``st.norm``) for ``beta``, with mean and standard deviation equal to 1.1 and .2, respectively. we also set the range of ``beta`` to be from 0.5 to 1.5. We do the same for phi.

The remaining lines just generate some simulated data to fit the model with, run the inference and plot the results which should include plots like this:

.. figure:: images/fit_series.png
   :width: 15cm
   
   Series posterior distributions. Colored areas represent 95% credible intervals.

.. figure:: images/fit_par.png
   :width: 15cm
   
    Parameters prior and posterior distributions.
    
    
One important argument in the run call, is the likvar, Which is the initial value for the likelihood variance. Try to increase its value if the acceptance ratio of the markov chain is too llow. Ideal levels for the acceptance ratio should be between 0.3 and 0.5.

The code for the above example can be found in the examples directory of the BIP distribution as :download:`deterministic.py <../../examples/Bayes/deterministic.py>`

Stochastic Model Example
------------------------

This example fits a stochastic model to simulated data. It uses the :ref:`SDE <SDE>` package of BIP:

.. literalinclude:: ../../examples/Bayes/flu_stochastic.py

This example can be found in the examples folder of BIP under the name of :download:`flu_stochastic.py <../../examples/Bayes/flu_stochastic.py>`.

Iterative Estimation and Forecast
=================================

In some other types of application, one's data accrue gradually and it may be interesting to use newly available data to improve previously obtained parameter estimations.

Here we envision two types of scenarios: one assuming constant parameters and another where parameter values can actually vary with time. These two scenarios lead to the two fitting strategies depicted on figure

.. figure:: images/Inference_scenarios2.png
    :width: 15cm
    
    Fitting scenarios: Moving windows and expanding windows.

References
==========
.. [1] Poole, D., & Raftery, A. E. (2000). Inference for Deterministic Simulation Models: The Bayesian Melding Approach. Journal of the American Statistical Association, 95(452), 1244-1255. doi:10.2307/2669764
