General Description
___________________

This directory contains the source code of the application diffExprPermutation.

diffExprPermutation Calculates P values for expression values per gene. Once the 
p value is calculated then is corrected by False Discovery Rate (Benjamini-Hochber). 


Permutation Test for difference in percentiles
______________________________________________

When we have a set of N samples (A of group 1 and B of group 2) and we want to test if
some of them are significant we must calculate the absolute difference in a chosen percentile
(median, 1Q, 2Q,...).

We must generate all possible permutations of group membership, and calculate the proportion of 
permutations with a difference >= to the observation.

Such operation implies a huge number of permutations:  n! / (a! b!)

The value of a given percentile depends only on one or two data points where we need to know 
were is which. Many fewer possible combinations of percentile values are needed.

We can cycle through all possible combinations in a reasonable time for most datasets.

Algorithm
_________
1. Find positions of percentile in the two samples, and the number of data points before and after the position.
2. Sort the two samples separately for each gene and calculate absolute percentile difference.
3. For each gene, merge the data for the two samples, keeping the final list sorted.
4. Arrange the data into vectors of length equal to the number of genes, indexed by the position in the sorted list.
5. Cycle through all possible values for the percentile position for sample 1. 
6. Cycle through all possible values for the percentile position for sample 2 conditional on the selected position for sample 1.
7. Calculate probability of this combination of positions.
8. Calculate absolute percentile difference for each gene using these positions.
9. Accumulate probability for all genes where observed difference grater or equal to permutation value. 

Percentile selection if the percentile falls between two ranks
--------------------------------------------------------------
The two flanking values are combined to get the percentile value.

Method Nearest: Closest rank to the chosen percentile. Only uses single values. Does not compute the average of two flanking 
values.
Method Linear: More stable when the number of samples is low. Can be more computational intensive for large sample numbers. (Track all 
combinations of three or four ranks rather than just two. For large sample sizes the difference between linear and nearest modes
are small.
In auto mode, Linear is used if the number of samples is less or equal to one hundred. Otherwise Nearest. 


Downloading
-----------

diffExprPermutation can be obtained from:
https://github.com/MarcosFernandez/diffExprPermutation


Compile
_______

You can compile the source code like:

1.Release
1.1. Compile
g++ -O3 -Wall -c -o main.o main.cpp
gcc --std=c99 -O3 -g -march=native -W -Wall -pedantic   -c -o perm.o perm.c

1.2.Link
g++ -o diffExprPermutation main.o perm.o  -lpthread -lm

2.Debug Just for testing.
2.1.Compile
g++ -O0 -g3 -Wall -c -o main.o main.cpp
gcc --std=c99 -O0 -g -march=native -W -Wall -pedantic   -c -o perm.o perm.c
2.2.Link 
g++  -o "diffExprPermutation"  main.o


3.Compile through make file
make clean
make


This code was compiled using gcc version 4.4.6 over 64 bits Redhat Linux and gcc 4.6.3 over 32 bits 
Ubuntu distribution.


Examples
________

Example A: Run a complete file (condition 1: case condition 2: control statistic median )
diffExprPermutation -f geneExpresionFile -o geneAdjuestedPvalues

Example B: Run a complete file (condition 1: group1 condition 2: group3 statistic perc75 )
diffExprPermutation -f geneExpresionFile -o geneAdjuestedPvalues --c1 group1 --c2 group3 -s per75



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
