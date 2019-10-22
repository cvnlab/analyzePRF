# forwardModel

forwardModel is a MATLAB toolbox for non-linear fitting of models to fMRI
data. It started life as analyzePRF (http://kendrickkay.net/analyzePRF/)
by Kendrick Kay. It has undergone refactoring to operate as an object-
oriented system for fitting models. Much of the clever structure that 
Kendrick had created (e.g., for cross-validation) was stripped out to
create a simpler, more general-purpose code set.

The legacy of Kendrick's code is most clearly seen in the model for pRF
mapping in retinotopic data. Specifically, the approach to creating
"super grid seeds", and the inclusion of a compressive non-linearity
in the modeled neural response, which is taken from:

  Kay KN, Winawer J, Mezer A and Wandell BA (2013) 
    Compressive spatial summation in human visual cortex.
    J. Neurophys. doi: 10.1152/jn.00105.2013

The implementation of pRF mapping implemented here differs from Kendrick's
original codes in a few ways:
  - The HRF is modeled as a double-gamma, which can be modified under the
    control of parameters
  - Multiple, cascading stages of non-linear fitting are supported with the
    ability to define sets of parameters that are fixed or float in a given
    search, with the results of a search passing to initialize the next
    stage.
  - For retinotopic mapping designs that play the same stimulus forward and
    reverse in time, the model will estimate a shift of the HRF time-to-peak
    to best fit the data.
  - Upper and lower bounds are enforced with an fmincon search.

More generally, the code supports the creation of object-oriented model
classes that can be evaluated within a common framework.
