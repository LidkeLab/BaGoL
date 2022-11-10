OVERVIEW:
 
This software package contains code and example scripts for the Baysian
Grouping of Localizations (BaGoL) analysis method described by

"Sub-Nanometer Precision using Bayesian Grouping of Localizations"
Mohamadreza Fazel, Michael J. Wester, David Shodt, Sebastian Restrepo Cruz, Sebastian
Strauss, Florian Schueder, Thomas Schlichthaerle Jennifer M. Gillette, Diane S. Lidke, 
Bernd Rieger, Ralf Jungmann, Keith A. Lidke

The Zip-file contains example scripts, data and a pre-compiled Windows mex
executable needed for frame connection.
The @BaGoL sub-folder contains a MATLAB class definition for the algorithm.

REQUIREMENTS: 

Windows 64 bit OS, 
MATLAB 64 bit, 
MATLAB Statistics and Machine Learning Toolbox.
The algorithm was tested using MATLAB R2017b and will likely run on any later
version.


INSTALLATION:
 
Unzip the files. 
In MATLAB, change directories to the BaGoL folder.


EXAMPLES: 

To run the demos, change directory to the 'Software Package' folder, open 
the scripts in MATLAB and run them or type the name of the scripts in a 
command window.  The 'Expected Results' folder contains the output produced
by the authors running the examples below. Due to the random nature of Monte 
Carlo technique, the BaGoL results won't be identical to those in the 
Expected_Results folder.
 
 
BaGoL_MPI_Origami.m
------------

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
-------

Analysis of localizations from a 4660x4660 nm^2 region of dSTORM EGFR data.  
This example demonstrates use of hierarchial Bayes to infer number of localizations 
per emitter distribution from part of the data and then use it to process the entire 
data set. It takes ~15 mins in total to run this script. The results will 
be saved in Results_EGFR folder. The outputs are similar to what was 
described previously with the inclusion of

MAPN_Hist+Random.png    Histogram of found NNDs compared to the curve for a
                        random distribution.

Eight_Mer.m: 
------------

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

