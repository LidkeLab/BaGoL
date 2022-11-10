# BaGoL

"High-Precision Estimation of Emitter Positions using Bayesian Grouping of Localizations"

Mohamadreza Fazel, Michael J. Wester, David J. Schodt, Sebastian Restrepo Cruz, Sebastian Strauss, Florian Schueder, 
Thomas Schlichthaerle, Jennifer M. Gillette, Diane S. Lidke, Bernd Rieger, Ralf Jungmann, Keith A. Lidke

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
produce posterior and MAPN images. The figure "Data_Flow" describes these steps. 

We tested several other common algorithm for the porpuse of grouping and combining of the localizations and BaGoL did better
than all of them. This method can be used for about a factor of two precision improvement on a typical dSTORM data set and
facilitate further quantitative analysis. When using DNA-PAINT, the method can achieve better than one nanometer precision.
We concieve numerous biological applications of the algorithm, such as inspection of protein-protein interactions, etc. 

# Software Package Description:

The software package contains code and example scripts for the Baysian Grouping of Localizations (BaGoL) analysis method.

Software Package:
The algorithm codes and a pre-compliled mex execuatable needed for frame connection.   
The @BaGoL sub-folder contains a MATLAB class definition for the algorithm

REQUIREMENTS:
Windows 64 bit OS
MATLAB 64 bit
MATLAB Statistics and Machine Learning Toolbox
The algorithm was tested using MATLAB R2018a and will likely run on any later version. 

INSTALLATION:
Download the Software Package.
In MATLAB change directories to the BaGoL folder. 

Data:
dSTORM data of EGF receptors, DNA Origami and simulated 8-mer data. These data are used by the demos. 
More data are available at the Nature Communication website published along the paper. 

#Demos: 
To run the demos open the scripts in MATLAB and run them or type the name of the scripts in a command window.

BaGoL_MPI_Origami.m:
Analysis of localizations from an experimental MPI DNA-Origami  
structure data. This show a basic hierarchical BaGoL data flow. 
The results will be saved in the Results_MPI folder.  Scale bars 
are 20 nm. Run time is ~5 min.  

The script should produce: 
SR_Im.png:                 Traditional super-resolution image. 
Post-Im.png:               Posterior image or histogram image of the chain
                           (weighted average over all models).
MAPN-Im.png:               MAPN image which is the image of localizations from the
                           most likely model. 
Overlay_SR_Map.png:        Overlay of grayscale SR-image and color MAPN image.
Overlay_SR_Post.png:       Overlay of grayscale SR-image and color posterior image. 
Overlay_SR_Map_circle.png: Overlay of the SR & MAPN coordinates where every coordinate
                           is represented by a circle located at the given location 
		           and a radius of double of the given precision.
Xi.png:                    Distribution of localizations per emitter.
NND.png:                   Histogram of nearest neighbor distances from
                           MAPN-coordinates. 
BaGoL_X-SE.png:            Histogram of X-localization precisions after grouping. 
BaGoL_Y-SE.png:            Histogram of Y-Localization precisions after grouping.
LocsScatter-MAPN.fig:      Plot of time color-coded localizations and
                           MAPN-coordinates.
MAPN.mat:                  Structure containing the MAPN-coordinates of emitters.

BaGoL_EGFR_dSTORM.m:
Analysis of localizations from a 4660x4660 nm^2 region of dSTORM EGFR data.  
This example demonstrates use of hierarchial Bayes to infer number of localizations 
per emitter distribution from part of the data and then use it to process the entire 
data set. It takes ~15 mins in total to run this script. The results will 
be saved in Results_EGFR folder. The outputs are similar to what was 
described previously with the inclusion of

MAPN_Hist+Random.png    Histogram of found NNDs compared to the curve for a
                        random distribution.

Eight_Mer.m: 
Animation of the RJMCMC chian for a simulated 8mer data set. 
It demonstrates the core BaGoL algorithm on a single cluster of localizations.
The results will be saved in the Results_Eight_Mer folder.  Scale bars
are 5 nm.  Run time is ~3 min.

The script should produce:
An animation of the chain,
PreBaGoL_SRImage.png:   Traditional super-resolution image.
Posterior_SRImage.png:  Super-resolution image from Posterior (weighted average
                        over all models).
MAPN_SRImage.png:       Super-resolution image from MAPN (most likely model).

# Input Parameters and Parameter Adjustments:
BaGoL has a few parameters that need to be carefully adjusted. A good description of the parameters are included
in the scripts documentation but they are also presented in the following. The unit for all the lengths are in nm.

SMD:                                                                                                                   
Structure containing input data with the following fields:                                                                       
   X:    Vector of X localizations (nm),                                                                                           
   Y:    Vector of Y localizations (nm),                                                                                        
   X_SE: Vector of X precisions (nm),                                                                                                    
   Y_SE: Vector of Y precisions (nm).                               

ROIsize:                                                                                                                    
The given coordinates are split into subregions with the size assigned to ROIsize for speed porpuses. The size of 
regions are inversely correlated with the density of localizations. (nm)
We recommend to pick the ROIsize so that there is not more than ~1000 localizations per ROI.

Overlap:                                                                                                                        
The size of overlapping region between adjacant subregions. Subregions are overlapped with their neighbors to 
avoid edge artifacts. The default value usually works for this parameter. (nm)
Overlap is often picked to be 10-20 nm depending on the ROIsize and the size of clusters.

Cutoff:                                                                                                                    
The localizations within each subregion are further divided into smaller set using hirerarchical algorithm as 
a pre-clustering algorithm. Cutoff is the size of the pre-clusters produced. (nm)
If your data is not too dense we suggest use whatever value larger than your ROIsize. If your dataset is dense then you
need to set it to a value smaller than your ROIsize so that the localizations within each ROI can be break further into
smaller pieces. However, "Cutoff" must not be too small so that this pre-clustering step starts breaking up clusters.

Drift:                                                                                                                                       
Drift may be presented in the data due to different reasons. BaGoL is able to handle movement of individual emitters
where Drift is the maximum movement of an emitter per frame. Use zero when there is no drift or residual drift. (nm/frame).

SE_Adjust:                                                                                                                                 
Localization precisions are often under-estimated in the loclization step. As such, we need to inflate the precisions by a
small value of SE_Adjust. Default is zero. (nm)

N_Burin:                                                                                                                            
Number of jumps within burnin portion of the chain for each pre-cluster.

N_Trials:                                                                                                                           
Number of the jumps within the post-burnin chain for each pre-cluster. The post-burnin part of the chain are returned  
for further analysis.  

PixelSize: 
The pixel size of the output images. (nm)

PImageSize:                                                                                                                           
Size of the produced posterior image, which is the same as the range of the input data set. (nm)

Xi:                                                                                                                                            
The algorithm can either learn this parameter from the data itself or take it as an input.
The inpout Xi can be either a scalar or a vector with two elements. Given a scalar value, BaGoL will implement a Poisson 
prior with mean value of Xi for average number of localizations per emitter. Given a vector with two elements,
BaGoL wil use a gamma prior for number of localizations per emitter. The product of the vector elements is equal 
to the average of the number of localizations per emitter. These two parameters gives the user the flexibility of 
adjusting the shape of the gamma distribution when the distribution shape is not well characterized. When learning Xi, the
given values will be used to initialize the corresponding chain. Again if Xi is a scalar it is used to initialize a Poisson prior
otherwise a gamma prior.

HierarchicalFlag:                                                                                                                             
0 do not learn Xi. 1 learn Xi. Default 0.

PImageFlag:                                                                                                                       
1 produces the posterior image. Default is 0.

ChainFlag:                                                                                                                 
1 saves the output chain. default is 0. It is recommended not to save the chain because it can take a very large chunk 
of the memory.

# Outputs:

MAPN:                                                                                             
Structure containing some results:                                                                               
   X:      Vector of found emitter X positions (nm),                                                                                                                                                                
   Y:      Vector of found emitter Y positions (nm),                                                                                          
   X_SE:   Vector of precisions for found X emitter positions (nm),                                                                          
   Y_SE:   Vector of precisions for found Y emitter positions (nm),   
   AlphaX: Vector of found X-drift velocities for each emitter (nm/frame),                                                                
   AlphaY: Vector of found Y-drift velocities for each emitter (nm/frame),                                                     
   Nmean:  Vector of mean number of localizations allocated to each found emitter.
   
PImage:
Posterior image

Chain:                                                                                                 
Cell array containing the BaGoL chain for each pre-cluster

XiChain:                                                                         
Chain of Xi samples

Note: The software package contain multiple functions to visualize and siplay the results including: 
makeIm(), dispIm(), genSRMAPNOverlay(), plotMAPN(), plotNND_PDF(), saveBaGoL().
The easiest way to generate reults is using the function "saveBaGoL()", which generates the following plots and images:
histogram of NND, histograms of X and Y precisions, plot of Xi chain, SR image using the input localizations, MAPN image
using the found emitter positions within the MAPN structure, Posterior image, Overlay image.
   
If you have questions, please feel free to shoot us an email:                                                                     
Mohamadreza Fazel: fazel.mohamadreza@gmail.com,                                                                                                                                                             
Michael Wester: wester@math.unm.edu,                                                                                                  
Keith Lidke: klidke@unm.edu.








