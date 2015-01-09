
3CDE Deconvolution Software:

We have provided the following 4 software:

3CDE: Approximate 3CDE deconvolution for either binary or linear case.
3CDEfrac: Approximate 3CDE deconvolution for fractional class densities.
3CDEilp: Exact 3CDE deconvolution for 
3CDEfracilp: Exact 3CDE deconvolution for

We explain how to run each in detail below.

A) INFORMATION ABOUT THE FILE FORMATS:

1- Frequency matrix file: Interaction matrix is gzipped. Each line starts with the binname such as 'node0', and continues with the number of interactions with the rest of nodes. Unzipped file will look like:

node0\t1\t0\t2\t.... 
node1\t1\t0\t1\t.... 
node2\t2\t3\1\t
....

where the number of interactions between node0 and node2 is 2. This matrix file is symmetric. Sample file is freq.gz

2- Prior domains file: Prior available domains file. These domains can be obtained from existing domain finders over ensemble data~(See Armatus). First bin is 1 not 0! The file is comma separated:

10,13
16,20
...

Sample file is domains.txt


B) Prerequisites:

3CDE have the following prerequisites:

a- IBM Cplex to be installed. It should be run by command 'cplex'.
b- Version of MATLAB newer than 2010a needs to be installed. It should be run by command 'matlab'.

Already provided:
a- SDPT3: Semidefinite Programming solver http://www.math.nus.edu.sg/~mattohkc/sdpt3.html. It is already provided in 'SDPT3-4.0' folder which is in the same directory as software.

If these software are installed in different locations, you can provide matlab and SDPT3 locations by the following options:

-sdpt SDPTPATH
-mat MATLABPATH


C) 3CDE: Deconvolves the ensemble interaction matrix in given number of classes in terms of domains. It uses the approximate method described in our paper. Class densities can be either uniform, or integer.

Sample 3CDErun:

python 3CDE.py -f freq.gz -d domains.txt -c 4 -o deconout -s binary

where parameters are:

-f: Ensemble interaction matrix 'freq.gz'
-d: Prior possible domains file if exists. For instance, 'domains.txt'
-c: Number of classes(default: 5)
-o: output prefix(default: deconout)
-s: scale type binary or integer(default: binary)
-p: Estimation type, ML or MAP problems(default: map)
-l: lambda for prior domain effects(default: 1.0)
-kern: Kernel type, exp or gauss. Leave blank for None(Default: None) 
-coef: Kernel coefficient such as 0.2. Leave blank for None(Default: None) 
-sdpt: SDPT Path(Default: SDPT3-4.0)
-mat: MATLAB Path(Default: matlab)

It generates 2 output files:

deconout_objscore.txt: This file stores the meta information of the method such as:

objlist	29111.0	27639.0
fqsum	29111.0
objval	27639.0
ratio	0.949434921507
time	21.607671976

where objlist is iterated objectives, objval is the last objective, and time is the overall running time

deconout_decon.txt: Deconvolution output. The format is as follows:

density of class 1\t1,3\t5,6\t.....\n
density of class 2\t2,5\t8,11\t.....\n

Each line is decomposition of a class. It starts with the estimated density, and continues with domain start,end positions seperated by tabs. If scale type is uniform, each class density becomes 1.


D) 3CDEfrac:

This method is similar to 3CDE except class densities are fractional instead of integer or uniform.

Sample 3CDEfracruns:

python 3CDEfrac.py -f freq.gz -c 5 -o deconfracout
python 3CDEfrac.py -f freq.gz -d domains.txt -c 5 -o deconfracout

where parameters are:

-f: Ensemble interaction matrix 'freq.gz'
-d: Prior possible domains file if exists. For instance, 'domains.txt'
-c: Number of classes(default: 5)
-o: output prefix(default: deconfracout)
-p: Estimation type, ML or MAP problems(default: map)
-l: lambda for prior domain effects(default: 1.0)
-kern: Kernel type, exp or gauss. Leave blank for None(Default: None) 
-coef: Kernel coefficient such as 0.2. Leave blank for None(Default: None) 
-sdpt: SDPT Path(Default: SDPT3-4.0)
-mat: MATLAB Path(Default: matlab)

It generates 2 output files similar to 3CDE. Please see the explanation of 3CDE.


E) 3CDE_ILP:

This method is an exact version of the approximate method 3CDE for integer and uniform class densities. Its parameters are similar to 3CDE.

python 3CDE_ILP.py -f freq.gz -c 5 -o deconilpout

It generates 2 output files similar to 3CDE. Please see the explanation of the approximate method 3CDE.


F) 3CDEfrac_ILP:

This method is exact method for fractional densities. Its parameters are similar to 3CDEfrac.

python 3CDEfrac_ILP.py -f freq.gz -c 5 -o deconilpfracout

It generates 2 output files similar to 3CDEfrac. Please see the explanation of the approximate methods 3CDEfrac and 3CDE.

