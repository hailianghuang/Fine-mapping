----------------------
Copyright
----------------------
Copyright (c) 2015, Hailiang Huang and Mark J Daly, Analytic and Translational Genetics Unit, Massachusetts General Hospital, Boston MA
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

----------------------
Synopsis
----------------------
A Fine-mapping method using flat prior with steepest descent approximation.  Method details available at: HTTP://biorxiv.org/content/early/2015/10/20/028688

----------------------
Prerequisites
----------------------
The R nnet package is required to run this method

----------------------
Usage
----------------------
A sample dataset was provided at "dat/", including:
1. example.raw: genotype data generated using the --recodeA option in PLINK (HTTP://pngu.mgh.harvard.edu/~purcell/plink/).  Each row is a sample and each column is a SNP
2. example.pheno: phenotype for the samples. 0 is control, 1 is UC and 2 is CD.
3. pcs.txt: covariate file. Each row is a sample and each column is covariate you want to adjust for. 
Fine-mapping can be performed by executing source("code/finemapping.r") in R console. Output files will be saved in the "result/", including:
1. repos.txt: detailed summary statistics. Can usually be ignored.
2. credible.txt: a list of SNPs in the 99% credible set. Most column headers should be self explanatory.  The "probNormed" column has the posterior probability for the SNPs, and the "R" column is the correlation coefficient between the SNP and the lead SNP in the signal. 










