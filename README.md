# EBSDrefine
Refine EBSD data by a modified dictionary indexing method. Looks at input orientations (presumably Hough-indexed data) to create a dictionary focused in orientation space and then interpolates to improve the accuracy of the indexed solution. A main application of this method is to resolve pseudosymmetry, but it can also be used to refine the orientations of non-pseudosymmetric materials.

Our paper describing this method and its performance can be found at: https://arxiv.org/abs/2003.04476 (preprint)

*Patterns used in this paper can be found in the folder 'paperdata/'*

This package uses EMsoft, which needs to be installed on your computer. More information about EMsoft can be found at: https://github.com/EMsoft-org/EMsoft

You will also need the following code packages:
1. *pcglobal*: https://github.com/epang22/pcglobal
  * Fit pattern centers using global optimization. 
2. *EMsoft-utilities*: https://github.com/epang22/EMsoft-utilities
  * Store your patterns in a .data file for EMsoft.
  * Figure out some parameters needed for EBSDrefine. *This is an important step, as some combinations of nthreads, numexptsingle/numdictsingle, and chunk (# orientations indexed at once) will cause EMsoft's EMEBSDDI program to give an error (bug in EMsoft)*
  * Combine data from multiple phases, which need to be indexed separately using EBSDrefine.

Program descriptions:
- RunEBSDrefine_checkangles_plotfit.m: Use this program to check what fit values you should restrict when reading in data.
- RunEBSDrefine_checkangles.m: Use this program to check list of compiled euler angles before running EBSDrefine to figure out what values for 'eulerminappear','angletol', 'fitmax', and 'reduceeuler' to use.
- RunEBSDrefine.m: Run EBSDrefine to index your data.

Feel free to email me at epang@mit.edu if you have any questions/difficulties/suggestions.
