---------------------Compressed Non-Negative CPD Toolbox--------------------------

Work by Jérémy E.Cohen, Rodrigo Cabral Farias and Pierre Comon.
To contact authors : 
- jeremy.cohen(at)gipsa-lab.fr
- rodrigo.cabral-farias(at)gipsa-lab.fr
- pierre.comon(at)gipsa-lab.fr

If you use this work, please cite :

E.Cohen, J., Farias, R. C., & Comon, P. (2015). Fast Decomposition of Large Nonnegative Tensors. IEEE Signal Processing Letters, 22(7), 862-866.


Last update : February 2015
-----------------------------------------------------------------------------------------------------

Why using the Compressed Non-Negative CPD toolbox for MATLAB?

For a given positive tensor T, this toolbox provides an implementation of a fast CPD based on compression and constrained CPD to obtain the non-negative factors of the multilinear model (a complete description is given in the reference at the beginning of this text file).

The constrained CPD itself is computed by two methods, one based on the gradient of the penalized objective function, the other one on a projected ALS. An uncompressed projected ALS is also provided. 

Our implementations of ALS does not use line search. Also, the compression is computed by three randomized SVDs. 

Missing data has to be handled separately. We provide a function that interpolates linearly along the third mode; if this does not bear any physical meaning, we encourage the user to complete the data randomly or by another interpolation paradigm.

This toolbox should run on all versions of MATLAB since R2009b (it uses the command ~ ), but can be easily adapted to older versions.

-----------------------------------------------------------------------------------------------------

Content : 

(For detailed comments on each function and arguments list, use the 'help' MATLAB function or open the .m files)

anls : computes the CPD with a simple projected ALS (no compression).

ccg : computes the CPD with the CCG algorithm. 

proco_als : computes the CPD with the ProCo-ALS algorithm.

eyeT : constructs a unit hyper diagonal tensor of given rank with given dimensions.

HOSVD : computes the high order SVD of the input tensor with given multilinear ranks.

kr : computes the khatri-rao product of two matrices.

miss_data : linearly interpolates missing data along the third mode.

MultProd : multiply a tensor by three matrices, one in each mode.

nlcg_loop : executes one loop of the Polak-Ribière conjugate gradient algorithm.

ParFold : transforms the parameter vector into three factor matrices.

ParUnfold : transforms the factors matrices into a parameter vector.

rsvd : computes the randomized SVD of the input matrix by an iterative algorithm from Tropp, 2009.

Run_decomp : the main file, from where to launch simulations.

sigmoid_diff : applies the derivative of a chosen sigmoid function to some input table with compression constraints.

sidmoid : applies the chosen sigmoid function to some input table.

tensor_gen : generates a random tensor with or without noise with given dimensions.

unvec : transforms a vector into a tensor, inverse function of 'vec'.

vec : transforms a tensor into a vector, row after row, front slice per front slice.

-------------------------------------------------------------------------------------------------------

How to use : 

The file 'Run_decomp' can be run as soon as the toolbox is downloaded. An example of syntax for personal data is given as comments l.23-26. This file contains the few parameters that one should fix before running the CPD algorithms : 

- R : tensor true rank (for simulations only)

- Re : tensor estimated rank

- dim_u : dimensions of the input data

- dim_c : multilinear ranks in the compression. Default values can be [Re+10,Re+10,Re+10] if the rank is not known, or [Re,Re,Re] if the rank is almost surely known. In any case, compressing more than the rank will generate major errors in the compression. Also, the SVDs make it impossible to do compression resulting in Rk>Ri*Rj with any rotations of indices.


CCG has moreover two parameters 'alpha' and 'gamma', respectively the stiffness and the magnitude of the penalization. Optimizing these parameters is important when using CCG for the CPD. 

-------------------------------------------------------------------------------------------------------

For any questions, remarks or improvement suggestions, please contact the authors at the following mail adresses :
- jeremy.cohen(at)gipsa-lab.fr
- rodrigo.cabral-farias(at)gipsa-lab.fr
- pierre.comon(at)gipsa-lab.fr


This work has been funded by the European Research Council under the European Community’s Seventh Framework Programme FP7/20072013 Grant Agreement no. 320594, ``DecoDA’’ project.
