# pcglobal
Algorithm for fitting pattern centers of EBSD patterns. Fits orientation and pattern center by optimization in 6D space using the SNOBFIT global optimization algorithm. Uses the EMsoft package to simulate EBSD patterns.

More information about EMsoft can be found at: https://github.com/EMsoft-org/EMsoft
More information about the SNOBFIT algorithm can be found at: https://www.mat.univie.ac.at/~neum/software/snobfit/

Publication coming soon!


Installation instructions:
1. Download the minq5 package from: https://www.mat.univie.ac.at/~neum/software/minq/minq5.tar.gz. Extract in a folder named 'minq5' located in the 'pcglobal' directory.
2. Install EMsoft. Code and installation instructions can be found at: https://github.com/EMsoft-org/EMsoft


Use instructions:
1. Before performing a fitting run, you need to precompute the master pattern using EMsoft. See here for instructions: https://github.com/EMsoft-org/EMsoft/wiki/EBSD-Example
2. Open 'RunPCfit.m'.
3. Enter the parameter values in the section labeled 'INPUT PARAMETERS'.
    a. homepath: This is EMsoft's EMdatapathname for your installation. If you are not sure where this is, you can find it at: ~/.config/EMsoft/EMsoftConfig.json
    b. Other parameters are explained in the comments in the script
4. Run the script.
*There is also a version 'RunPCfit_loop.m' that allows you to loop through multiple image files.


I have included some examples to learn how to use the package. This data is located in the folder 'testdata'. You need to take a few additional steps to be able to perform a fitting run on this data:
1. Move the Ni.xtal file to your EMsoft xtal folder. If you are not sure where this is, you can find it at: ~/.config/EMsoft/EMsoftConfig.json
2. Download Ni.tar.gz from: https://kilthub.cmu.edu/articles/Data_Sets_and_Analysis_Results_Accompanying_a_Dictionary_Indexing_Tutorial_Paper/7792505. Copy the file 'Ni-master-20kV.h5' to any location within EMsoft's EMdatapathname.


Feel free to email me at epang22@gmail.com if you have any questions/difficulties/suggestions.
