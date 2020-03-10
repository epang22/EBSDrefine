# EBSDrefine
Refine EBSD data by a modified dictionary indexing method. Looks at input orientations (presumably Hough-indexed data) to create a dictionary focused in orientation space and then interpolates to improve the accuracy of the indexed solution. A main application of this method is to resolve pseudosymmetry, but it can also be used to refine the orientations of non-pseudosymmetric materials.

Our paper describing this method and its performance can be found at: (preprint)

*Patterns used in this paper can be found in the folder 'paperdata/'*

This package uses EMsoft, which needs to be installed on your computer. More information about EMsoft can be found at: https://github.com/EMsoft-org/EMsoft

Recommended code packages:
1. *pcglobal*. Fit pattern centers using global optimization. https://github.com/epang22/pcglobal
2. *EMsoft-utilities*. Use it to store your patterns in a .data file for EMsoft. Also use it to combine data from multiple phases, which need to be indexed separately using EBSDrefine. https://github.com/epang22/EMsoft-utilities

Program descriptions:
- RunEBSDrefine_checkangles_plotfit.m: Use this program to check what fit values you should restrict when reading in data.
- RunEBSDrefine_checkangles.m: Use this program to check list of compiled euler angles before running EBSDrefine to figure out what values for 'eulerminappear','angletol', 'fitmax', and 'reduceeuler' to use.
- RunEBSDrefine.m: Run EBSDrefine to index your data.

Feel free to email me at epang@mit.edu if you have any questions/difficulties/suggestions.
