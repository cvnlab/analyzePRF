# analyzePRF

analyzePRF is a MATLAB toolbox for fitting population receptive field (pRF) models
to fMRI data.  It is developed by Kendrick Kay (kendrick@post.harvard.edu).

The toolbox has several dependencies:
- MATLAB Optimization Toolbox
- GLMdenoise (necessary only if you use the GLMdenoise feature; to download
              GLMdenoise, see http://kendrickkay.net/GLMdenoise/)

To use the toolbox, add it to your MATLAB path:
  addpath(genpath('analyzePRF'));

To try the toolbox on an example dataset, change to the analyzePRF directory 
and then type:
  example1;
This script will download the example dataset (if it has not already been
downloaded) and will run the toolbox on the dataset.

For additional information, please visit:
  http://kendrickkay.net/analyzePRF/

Terms of use: This content is licensed under a Creative Commons Attribution 3.0 
Unported License (http://creativecommons.org/licenses/by/3.0/us/). You are free 
to share and adapt the content as you please, under the condition that you cite 
the appropriate manuscript (see below).

If you use analyzePRF in your research, please cite the following paper:
  Kay KN, Winawer J, Mezer A and Wandell BA (2013) 
    Compressive spatial summation in human visual cortex.
    J. Neurophys. doi: 10.1152/jn.00105.2013

History of major code changes:
- 2014/06/17 - Version 1.1.

## CONTENTS

Contents:
- analyzePRF.m - Top-level function that you want to call
- analyzePRFcomputeGLMdenoiseregressors.m - Helper function
- analyzePRFcomputesupergridseeds.m - Helper function
- example*.m - Example scripts
- exampledataset.mat - Example dataset
- README - The file you are reading
- setup.m - A simple script that downloads the example dataset
            and adds analyzePRF to the MATLAB path
- utilities - A directory containing various utility functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Copyright (c) 2014, Kendrick Kay
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

The names of its contributors may not be used to endorse or promote products 
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
