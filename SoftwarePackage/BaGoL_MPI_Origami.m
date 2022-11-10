%% Bayesian Grouping of Localizations (BaGoL) Example for MPI DNA-Origami
%  BaGoL is run for MPI DNA-origami data. We use a hierarchical Bayes 
%  procedure to learn the average number of localizations per emitter, \xi,
%  simultaneously along with other parameters.
%
% Requirements and Setup:
%   1. Windows 64 bit OS 
%   2. MATLAB 2017b or higher versions
%   3. Statistics and Machine Learning Toolbox
%   4. BaGoL class
%   5. Set MATLAB directory to the Software Package directory.
%
% Results include:
%   Saved Results:
%     SR_Im.png:                 Traditional super-resolution image. 
%     Post-Im.png:               Posterior image or histogram image of the chain
%                                (weighted average over all models).
%     MAPN-Im.png:               MAPN image which is the image of localizations from the
%                                most likely model. 
%     Overlay_SR_Map.png:        Overlay of grayscale SR-image and color MAPN image.
%     Overlay_SR_Post.png:       Overlay of grayscale SR-image and color posterior image. 
%     Overlay_SR_Map_circle.png: Overlay of the SR & MAPN coordinates where 
%                                every coordinate is represented by a circle  
% 		                         located at the given location and a radius 
%                                of double of the given precision.
%     Xi.png:                    Number of localizations per emitter dist.
%     NND.png:                   Histogram of nearest neighbor distances from
%                                MAPN-coordinates. 
%     BaGoL_X-SE.png:            Histogram of X-localization precisions after grouping. 
%     BaGoL_Y-SE.png:            Histogram of Y-Localization precisions after grouping.
%     LocsScatter-MAPN.fig:      Plot of time color-coded localizations and
%                                MAPN-coordinates.
%     MAPN.mat:                  Structure containing the MAPN-coordinates of emitters.
%
%   Output available on work space:
%     MAPN: Clusters information are stored in this property:
%     MAPN.X: X-Centers (nm)
%     MAPN.Y: Y-Centers (nm)
%     MAPN.X_SE: X-Centers precisions (nm)
%     MAPN.Y_SE: Y-Centers precisions (nm)
%     MAPN.AlphaX: X-Drifts of clusters (nm/frame)
%     MAPN.AlphaY: Y-Drifts of clusters (nm/frame)
%     MAPN.AlphaX_SE: X-Drift precisions (nm/frame)
%     MAPN.AlphaY_SE: Y-Drift precisions (nm/frame)
%     MAPN.Nmean: Mean number of binding events per docking strand

warning('OFF', 'stats:kmeans:FailedToConvergeRep');

%% Important Parameters

Xi = [1 60]; %shape and scale parameters of gamma prior on number of localizations per emitter
SE_Adjust = 0.4; %adjustment of localization accuracies (nm)
ClusterDrift = 0; %emitters drift velocities (nm/frame)
ImSize = 100; %Size of output images (nm)
PixSize = 0.5; %pixel size of output images (nm) 
XStart = 30085; %minimum of X-coordinates (nm)
YStart = 58105; %minimum of Y-coordinates (nm)
ROISZ = 100; %data is split into regions with this size (nm)
Overlap = 0; %size of overlapping regions
NBurnin = 20000; %Number of burn-in samples
NTrils = 10000; %Number of samples after burn-in

SaveDir = fullfile('Results_MPI'); %Save director
if ~isdir(SaveDir)
   mkdir(SaveDir) 
end

%% Load data

DataDir = fullfile('Data'); %Data directory
FileName = 'SMD_DNA-Origami_MPI.mat'; %Data file
load(fullfile(DataDir,FileName))

%% Filter based on number of localizations within 3 median of precisions

Prec_Median = median([SMD.X_SE;SMD.Y_SE]);
[~,D]=knnsearch([SMD.X,SMD.Y],[SMD.X,SMD.Y],'K',length(SMD.X));
D(:,1)=[];

ID = D < 3*Prec_Median;
N = sum(ID,2);
Ind = N>24;

SMDFilt.X = SMD.X(Ind);
SMDFilt.Y = SMD.Y(Ind);
SMDFilt.X_SE = SMD.X_SE(Ind);
SMDFilt.Y_SE = SMD.Y_SE(Ind);
SMDFilt.FrameNum = SMD.FrameNum(Ind);

%% Setting the class properties and running BaGoL

MPI = BaGoL;
MPI.SMD = SMDFilt;
MPI.ROIsize = ROISZ;
MPI.Overlap = Overlap;
MPI.Xi = Xi;
MPI.Beta_Xi = 75;
MPI.N_Burnin = NBurnin;
MPI.N_Trials = NTrils;
MPI.Drift = ClusterDrift;
MPI.SE_Adjust = SE_Adjust;
MPI.PImageFlag = 1; %Produce Posterior image
MPI.HierarchFlag = 1; %learn Xi
MPI.NSamples = 20;
MPI.PImageSize = ImSize; 
MPI.PixelSize = PixSize; 
MPI.XStart = XStart;
MPI.YStart = YStart;

%Analyzing the data
tic;
MPI.analyze_all();
T = toc;
fprintf('It took %g seconds to process the MPI DNA-Origami.\n',T)

%Generating output plots and images.
ScaleBar = 20; %length of scale bars (nm)
OverlayFlag = 1;
RadiusScale = 2;
MPI.saveBaGoL(ScaleBar,SaveDir,OverlayFlag);
MPI.plotMAPN(SaveDir)
PixSizeC = 'rescale';
BaGoL.genSRMAPNOverlay(MPI.SMD, MPI.MAPN, ImSize, ImSize, PixSizeC, SaveDir, ...
                       XStart, YStart, RadiusScale, ScaleBar);
     
%%

SaveDir = fullfile('Results_MPI_NoFilt'); %Save director
if ~isdir(SaveDir)
   mkdir(SaveDir) 
end

MPI.SMD = SMD;
tic;
MPI.analyze_all();
T = toc;
fprintf('It took %g seconds to process the MPI DNA-Origami.\n',T)

%Generating output plots and images.
ScaleBar = 20; %length of scale bars (nm)
OverlayFlag = 1;
RadiusScale = 2;
MPI.saveBaGoL(ScaleBar,SaveDir,OverlayFlag);
MPI.plotMAPN(SaveDir)
PixSizeC = 'rescale';
BaGoL.genSRMAPNOverlay(MPI.SMD, MPI.MAPN, ImSize, ImSize, PixSizeC, SaveDir, ...
                       XStart, YStart, RadiusScale, ScaleBar);

                   