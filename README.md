General Description
___________________

This directory contains the source code of the application diffExprPermutation.

diffExprPermutation Calculates P values for expression values per gene. Once the 
p value is calculated then is corrected by False Discovery Rate (Benjamini-Hochber). 

diffExprPermutation has three modes of execution (SPLIT,PERM,FDR). You can run directly
PERM mode over your input data and you will get the expected results (pValues and Adjusted pValue).
If you want to speed up the process you must use SPLIT and FDR modes. SPLIT creates chunks of genes
with the goal of performing a permutation test in parallel. Per each created chunk you can run the application 
in PERM mode. Once PERM step is performed then FDR must be called to get all genes in a file and an adjusted pvalue
per Benjamini-Hochber method. 

Permutation Test
________________

To calculate the p value it is used a MonteCarlo Sampling Permutation Test. A permutation test is a type of statistical 
significance test in which the distribution of the test statistic under the null hypothesis is obtained 
by calculating all possible values of the test statistic under rearrangements of the labels on the 
observed data points. In other words, the method by which treatments are allocated to subjects in an 
experimental design is mirrored in the analysis of that design. If the labels are exchangeable under 
the null hypothesis, then the resulting tests yield exact significance levels.

The test proceeds as follows. First, the difference in medians between the two samples is calculated: 
this is the observed value of the test statistic, T(obs). Then the observations of groups A and B are pooled.

Next, the difference in sample medians is calculated and recorded for every possible way of dividing these 
pooled values into two groups of size n_{A} and n_{B}. 
The set of these calculated differences is the exact distribution of possible differences under the null hypothesis 
that group label does not matter.

The one-sided p-value of the test is calculated as the proportion of sampled permutations where the difference 
in medians was greater than or equal to T(obs). The two-sided p-value of the test is calculated as the proportion 
of sampled permutations where the absolute difference was greater than or equal to ABS(T(obs)).

If the only purpose of the test is reject or not reject the null hypothesis, we can as an alternative sort 
the recorded differences, and then observe if T(obs) is contained within the middle 95% of them. If it is not, 
we reject the hypothesis of identical probability curves at the 5% significance level.

Monte Carlo testing
___________________

An asymptotically equivalent permutation test can be created when there are too many possible orderings of the 
data to allow complete enumeration in a convenient manner. This is done by generating the reference distribution 
by Monte Carlo sampling, which takes a small (relative to the total number of permutations) random sample of 
the possible replicates. 

After N random permutations, it is possible to obtain a confidence interval for the p-value based on 
the Binomial distribution. 

For example, if after N = 10000 random permutations the p-value is estimated to be {p}=0.05 , then a 99% confidence 
interval for the true p (the one that would result from trying all possible permutations) is  [0.044, 0.056].

The assessment of statistical significance by Monte Carlo simulation may be costly in computer time. The criteria to 
stop the number of permutations for a given gene and reduce the computing time is based on the methods proposed in 
Sequential Monte Carlo p-Values (Julian Besag and Peter Clifford).

P Values are adjusted by Benjamini-Hochberg procedure.

Downloading
-----------

diffExprPermutation can be obtained from:
https://github.com/MarcosFernandez/diffExprPermutation


Compile
_______

You can compile the source code like:

(Release)
1.Compile
g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"main.d" -MT"main.d" -o "main.o" "main.cpp"
2.Link
g++  -o "diffExprPermutation"  main.o

(Debug) Just for testing.
1.Compile
g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"main.d" -MT"main.d" -o "main.o" "main.cpp"
2.Link 
g++  -o "diffExprPermutation"  main.o

This code was compiled using gcc version 4.4.6 over 64 bits Redhat Linux and gcc 4.6.3 over 32 bits 
Ubuntu distribution.


Examples
________

Example A: Run a complete file
diffExprPermutation -m PERM -f geneExpresionFile -o geneAdjuestedPvalues

Exemple B: Create chunks of files to perform a parallel running in case you have a cluster 
environment or a queue manager (SLURM, SGE, Torque...)

B.1.SPLIT your input file in small files where its file contains a set of genes
diffExprPermutation -m SPLIT -f geneExpresionFile -o chunk --n 4 

B.2.Call PERM mode per each chunk
diffExprPermutation -m PERM -f chunk_1 -o permutation/chunk_1
diffExprPermutation -m PERM -f chunk_2 -o permutation/chunk_2
diffExprPermutation -m PERM -f chunk_3 -o permutation/chunk_3
diffExprPermutation -m PERM -f chunk_4 -o permutation/chunk_4

B.3.Merge Results and adjust p-value by Benjamini-Hochberg
diffExprPermutation -m FDR -d permutation/ -o geneAdjuestedPvalues


-------------------------------------------------------------------------------
Copyright (C) 2009-2015  Centro nacional de análisis genómico.

diffExprPermutation is free software; you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

diffExprPermutation is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.
