# GrapheneFocusing

This repository contains all the code to generate the plots 
contained in the paper "Valley separation via trigonal warping". 
There are two separate pieces of code. The first is used to 
generate basic trajectories. The second is used to create the 
intensity pictures and the multitrajectory analysis with the grounded gates. 

The code is not very efficient, as the implementation of the grounded
gates just runs through a loop, and thus requires looping over many 
randomized trajectories to build a stable result for the intensity. 
