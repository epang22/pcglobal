# pcglobal
Algorithm for fitting pattern centers of EBSD patterns. Fits orientation and pattern center by optimization in 6D space using the SNOBFIT global optimization algorithm. This helps prevent getting trapped in local minima, which is an issue for standard downhill search methods. Uses EMsoft to simulate EBSD patterns. 

Our paper describing this method and its performance can be found at: https://arxiv.org/abs/1908.10692 (preprint) and https://www.sciencedirect.com/science/article/pii/S030439911930292X (final version in Ultramicroscopy)

*Patterns used in this paper can be found in the folder 'paperdata/'*

More information about EMsoft can be found at: https://github.com/EMsoft-org/EMsoft

More information about the SNOBFIT algorithm can be found at: https://www.mat.univie.ac.at/~neum/software/snobfit/

Installation instructions:
1. Install EMsoft. Code and installation instructions can be found at: https://github.com/EMsoft-org/EMsoft
2. Download/clone this repository in any location.
3. Download the minq5 package from: https://www.mat.univie.ac.at/~neum/software/minq/minq5.tar.gz. 

If you are using EMsoft v4.0, extract in a folder named 'minq5/' located in the 'pcglobal/EMsoft_v4_0/' directory. Only use the version of the code in this directory.

If you are using EMsoft v4.3, extract in a folder named 'minq5/' located in the 'pcglobal/EMsoft_v4_3/' directory. Only use the version of the code in this directory.

*Different versions of EMsoft have separate code here because of the differences in pattern orientation. Using the incorrect version can lead to low dot product values (and may not give an error) because the incorrect pattern orientation. We have only tested EMsoft v4.0 and v4.3. Please let us know if one or neither of these versions seems to work for you.*



Use instructions:
1. Before performing a fitting run, you need to precompute the master pattern using EMsoft. See here for instructions: https://github.com/EMsoft-org/EMsoft/wiki/EBSD-Example
2. Open 'RunPCfit.m'. There is also a version 'RunPCfit_loop.m' that allows you to loop through multiple image files.
3. Enter the parameter values in the section labeled 'INPUT PARAMETERS'. Comments in the script explain the meaning of each parameter. Default SNOBFIT values should be fine in most cases, but we have included a program FindParameters.m that can help you determine the optimal values for some of these parameters for your system.
4. Run the script.


I have included some examples to learn how to use the package. This data is located in the folder /testdata/. You need to take a few additional steps to be able to perform a fitting run on this data:
1. Download 'Ni.tar.gz' from: https://kilthub.cmu.edu/articles/Data_Sets_and_Analysis_Results_Accompanying_a_Dictionary_Indexing_Tutorial_Paper/7792505. Copy the file 'Ni-master-20kV.h5' to any location within EMsoft's EMdatapathname.
2. Download 'Xtals.tar.gz' from: https://kilthub.cmu.edu/articles/Data_Sets_and_Analysis_Results_Accompanying_a_Dictionary_Indexing_Tutorial_Paper/7792505
Copy the 'Ni.xtal' file to your EMsoft xtal folder. 
3. Copy the folder 'testdata/' to any location within EMsoft's EMdatapathname.
*You can find the location of your xtal and EMdatapathname folders at: '~/.config/EMsoft/EMsoftConfig.json'*


Feel free to email me at epang@mit.edu if you have any questions/difficulties/suggestions.
