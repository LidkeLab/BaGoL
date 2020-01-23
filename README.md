# BaGoL

"Sub-Nanometer Precision using Bayesian Grouping of Localizations"
Mohamadreza Fazel, Michael J. Wester, Sebastian Restrepo Cruz, Sebastian Strauss, 
Florian Schueder, Jennifer M. Gillette, Diane S. Lidke, Bernd Rieger, Ralf Jungmann, Keith A. Lidke

# Algorithm Overview:

Single molecule localization microscopy super-resolution methods such as DNA-PAINT and (d)STORM generate multiple observed 
localizations over the time course of data acquisition from each dye or binding site that are nor a priori assigned to 
those specific dyes or binding sites. BaGoL implements a Bayesian method of grouping and combining localizations from 
multiple blinking/binding events that can improve localization precision to better than one naometer. BaGoL allows 
inclusion of prior knowledge such as distribution of the number of localizations per emitter and the localization precisions.

The algorithm is comprised of several steps depicted in the figure. First, the list of localizations are split into smaller
subsets. Second, the outliers are recognized as localizations with less than a certain number of neighbors within a certain 
distance. Third, localizations within each subset are further split into preclusters using hierarchical clustering 
algorithm. Fourth, each precluster is processed using RJMCMC. Fifth, the chain from all the preclusters are combined to 
produce posterior and MAPN images.

We tested several other common algorithm for the porpuse of grouping and combining of the localizations and BaGoL did better
than all of them. This method can be used for about a factor of two precision improvement on a typical dSTORM data set and
facilitate further quantitative analysis. When using DNA-PAINT, the method can achieve better than one nanometer precision.
We concieve numerous biological applications of the algorithm, such as inspection of protein-protein interactions, etc. 

<p align="center"><img src="Data Flowt.png" width="80%" height="80%"></p>

# Software Package Description:

The software package contains code and example scripts for the Baysian Grouping of Localizations (BaGoL) analysis method described by: 

Software Package:
The algorithm codes and a pre-compliled mex execuatable needed for frame connection.   
The @BaGoL sub-folder contains a MATLAB class definition for the algorithm

REQUIREMENTS:
Windows 64 bit OS
MATLAB 64 bit
MATLAB Statistics and Machine Learning Toolbox
The algorithm was tested using MATLAB R2016a and will likely run on any later version. 

INSTALLATION:
Download the Software Package.
In MATLAB change directories to the BaGoL folder. 

Data:
dSTORM data of microtubules and EGF receptors and simulated 8-mer data. These data are used by the demos. 

Demos: 
To run the demos open the scripts in MATLAB and run them or type the name of the scripts in a command window.

Microtubules.m: 
Analysis of localizations from a 4800x4800 nm^2 region of dSTORM microtubule data. 
This show a basic BaGoL data flow using a broad gamma prior for localizations per emitter. 
The results will be saved in the BaGoL\Results_Microtubule folder. Scale bars are 1000 nm.
Run time is ~30 min. 

The script should produce: 

PreBaGoL_SRImage.png: 	Traditional Super-resolution Image
PreBaGoL_SRImage_Filtered.png:	Traditional Super-resolution Image after filtering by NN and intensity
Posterior_SRImage.png:	Super-resolution image from Posterior (weighted average) output
MAPN_SRImage.png: 	Super-resolution image from MAPN (most likely) output
Lambda_Hist.png:	Localizations per emitter for prior (curve) and found (histogram)
BaGoL_X-SE.png: 	Localization precision after grouping
BaGoL_Y-SE.png: 	Localization precision after grouping

EGFR.m: 
Analysis of localizations from a 4660x4660 nm^2 region of dSTORM data of A647-EGF bound to EGFR.
This shows a BaGoL data flow where localizations per emitter are estimated from data in a two-step process.  
The results will be saved in BaGoL\Results_EGFR folder. Scalebars are 1000 nm.
Run time is ~20 min. 

The script should produce: 

PreBaGoL_SRImage.png: 	Traditional Super-resolution Image
PreBaGoL_SRImage_Filtered.png:	Traditional Super-resolution Image after filtering by NN and intensity
Posterior_SRImage.png:	Super-resolution image from Posterior (weighted average) output
MAPN_SRImage.png: 	Super-resolution image from MAPN (most likely) output
Lambda_Hist_PriorEst.png:	Localizations per emitter for prior (curve) and found (histogram) for prior estimation
Lambda_Hist.png:	Localizations per emitter for prior (curve) and found (histogram) used for final analysis
BaGoL_X-SE.png: 	Localization precision after grouping
BaGoL_Y-SE.png: 	Localization precision after grouping
NND_Hist.png: 		Nearist Neighbor Distribution Histogram
NND_Hist+Random.png:	Nearist Neighbor Distribution Histogram with theorical curve for spatial Poisson random data


Eight_Mer.m: 
Animation of the jumps in the RJMCMC step for a simulated 8mer data set.
Demonstrates the core BaGoL algorithm on a single cluster of localizations.
The results will be saved in the BaGoL\Results_Eight_Mer folder. Scale bars are 5 nm.

The script should produce: 

An animation of the classification and emitter movement steps
PreBaGoL_SRImage.png: 	Traditional Super-resolution Image
Posterior_SRImage.png:	Super-resolution image from Posterior (weighted average) output
MAPN_SRImage.png: 	Super-resolution image from MAPN (most likely) output

# Parameters and parameters adjusment:
BaGoL has a few parameters that need to be carefully adjusted. A good description of the parameters are included
in the scripts documentation but they are also presented in the following. The unit for all the lengths are in nm.

Lambda:
Lambda can be either a scalar or a vector with two elements. Given a scalar value, BaGoL will implement a Poisson 
prior with mean value of Lambda for average number of localizations per emitter. Given a vector with two elements,
BaGoL wil use a gamma prior for number of localizations per emitter. The product of the vector elements is equal 
to the average of the number of localizations per emitter. These two parameters gives the user the flexibility of 
adjusting the shape of the gamma distribution when the distribution shape is not well characterized.

ROIsize: 
The given coordinates are split into subregions with the size assigned to ROIsize for speed porpuses. The size of 
regions are inversely correlated with the density of localizations. (nm)

Overlap:
The size of overlapping region between adjacant subregions. Subregions are overlapped with their neighbors to 
avoid edge artifacts. The default value usually works for this parameter. (nm)

Cutoff:
The localizations within each subregion are further divided into smaller set using hirerarchical algorithm as 
a pre-clustering algorithm. Cutoff is the size of the pre-clusters produced. (nm)

NNN and NNR:
To remove the ouliers within the list of localizations, BaGoL implement a nearest neighbor filter. The localizations
with less than NNN neighbors within NNR distance are considered outliers. These can be adjusted by inspecting the 
nearest neighbor distribution of the data.

Drift:
Drift may be presented in the data due to different reasons. BaGoL is able to handle movement of individual emitters
where Drift is the maximum movement of an emitter per frame. Use zero when there is no drift or residual drift. (nm/frame)

N_Burin: 
Number of jumps within burnin portion of the chain for each pre-cluster.
N_Trials: 
Number of the jumps within the post-burnin chain for each pre-cluster. The post-burnin part of the chain are returned  
for further analysis.  

P_Jumps: 
Probabilities of proposing different jumps [Move Allocate Add Remove], add up to one.

PImageFlag:
1 produces the posterior image. Default is 0.

PixelSize: 
The pixel size of the output images. (nm)

PImageSize:
Size of the produced posterior image, which is the same as the range of the input data set. (nm)

ChainFlag:
1 saves the output chain. default is 0. It is recommended not to save the chain because it can take a very large chunk 
of the memory.






