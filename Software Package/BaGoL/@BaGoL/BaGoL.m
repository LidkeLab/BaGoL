classdef BaGoL < handle
    %Bayesian Grouping of Localizations (BaGoL) takes advantage of the
    %localization accuracies and the prior information on the number of
    %blinking/binding events to identify groups of localizations that are
    %associated to specific fluorophores/docking strands. It then employs
    %those groups of localizations to improve the localizations acuuracies.
    %BaGoL splits the coordinates into smaller sub-regions for speed 
    %porpuses. Next, it uses nearest neighbor distance to remove the 
    %outliers. BaGoL then takes advantage of heirarchical clustering, which  
    %is a fast matlab function, to identify preclusters, then the
    %preclusters are sent to the core RJMCMC algorithm to find the groups.
    %  
    % INPUTS:
    %   SMD:        Structure containing inputs such as localization
    %               coordinates, accuacies, time information. The unit for 
    %               the coordinates and accuracies are nm.
    %               SMD.X, SMD.Y, SMD.Z: localization coordinates.
    %               SMD.X_SE, SMD.Y_SE, SMD.Z_SE: localization accuracies.
    %               SMD.T: time (frmae) information.
    %   Lambda:     This can be either a scalar or an array with two 
    %               elements. Scalar would be the average for a Poisson 
    %               prior on the number of localizations per emitters. 
    %               Array of two number will be used as parameters of a 
    %               gamma prior on number of localizations per emitters.
    %   Cutoff:     BaGoL uses a simple and fast clustering algorithm,
    %               which is a matlab function, for some pre-processing. 
    %               This parameter determines the size of the clusters from 
    %               this algorithm. (Default: 1 nm)
    %   NNN:        Number of Nearest Neighbors. BaGoL uses nearest  
    %               neighbor distance to remove the outliers (background).  
    %               This parameter is number of localizations within a  
    %               given radius. The localizations with less number of  
    %               neighbors within this distance will be considered as 
    %               outlier. (Default=5).
    %   NNR:        Nearest Neighbor Radius. There should be at least NNN  
    %               number of localizations within this radius to be     
    %               considered as signal. (nm),(Default = Inf)
    %   ROIsize:    The given coordinates are split to small regions for
    %               speed porpuses. ROIsize is the size of those regions.
    %               (Default: 200 nm)
    %   Overlap:    The sub-regions are overlapped with their adjacent sub-
    %               regions to avoid artifacts at the edges. This parameter
    %               is the size of overlaping regions. (Default: 50 nm)  
    %   Drift:      This is the average drifting of clusters per frame. In 
    %               other words, this is the entire movement of an emitter 
    %               divided by number of the frames. (Default=0) (nm/frame)
    %   N_Burnin:   Number of the jumps in the burnin part of the chain.
    %   N_Trials:   Number of the jumps in the chain after the burnin part. 
    %   P_Jumps:    Probabilities for proposing different jumps.
    %               [Reassign, Redistribute, Birth, Death], the sum
    %               of these probabilities adds up to one. (Default: [0.25, 
    %               0.25, 0.25, 0.25])
    %   PixelSize:  The pixel size for posterior image. (Default: 1 nm)
    %   PImageFlag: The flag for the Posterior image. Since the posterior
    %               image can be very large specially when we have a 
    %               sub-nanometer pixel size, its size might be very large
    %               and the machine will run out of memory. (Default: 0).
    %   PImageSize: Size of the Posterior image. (nm) This usually matches
    %               the size of the area that will be processed. (Scalar)
    %   ChainFlag:  The chain is not saved by default. User may change this
    %               parameter to save the chain.
    % OUTPUTS:
    %   MPAN:       Structure containing the group locations and their 
    %               associated accuracies extracted from the most repeated
    %               model in the chain.
    %               MAPN.X, MAPN.Y, MAPN.Z: groups localizations
    %               MAPN.X_SE, MAPN.Y_SE, MAPN.Z_SE: accuracies
    %               MAPN.Nmean: mean number of localizations in the groups
    %   PImage:     Posterior image containing all the possible models
    %   Chain:      MATLAB-cell containing the accepted chains after burnin 
    %
    % REQUIRES:
    %   MATLAB 2016 or higher versions.
    %   Statistics and Machine learning toolbox.
    %
    % Created by:
    %   Mohamadreza Fazel (Lidke lab 2019)
    %
    properties
        SMD %Single Molecule Data structure (SMD.X, SMD.Y, ..., SMD.X_SE)
        ClusterSMD %Array of SMD - one per cluster from heirarchical 
                   %clustering.   
        MAPN %Structure containing group locations and accuracies.
        PImage; %Posterior image
        PImageFlag = 0; %Flag to make posterior image
        PImageSize = []; %Size of the posterior image (nm)
        PixelSize = 1; %pixel size for posterior images (nm).
        P_Jumps=[0.25 0.25 0.25 0.25]; %Probabilities of proposing different jumps
        N_Burnin=500; %Number of jumps in burnin chain.
        N_Trials=500; %Number of jumps in post-burnin chain.
        ROIsize=200; %Size of the ROIs (nm)
        Overlap=50; %Width of the overlaping region for ROIs (nm).
        NNR = Inf; %Nearest Neighbor Radius for filtering
        NNN = [];  %Number of Nearest Neighbors 
        Drift = 0; %Drift used to inflate the localization precisions and 
                   %to claculate the prior on slope of the drift in
                   %th likelihood. (nm)
        Cutoff=1; %cutoff size of the clustering in preclustering, used 
                  %for heirarchical clustering (nm).
        Lambda; %This can be either a scalar or a vector. Scalar means the
               %parameter for Poisson dist. and vector means array of 
               %distribution starting from zero and increment by one.
        ChainFlag = 0; %Save chain or not.
        Chain = {}; %The cell array of the accepted chain.
    end
    
    methods
        
        function obj=analyze_all(obj)
            %analyze_all() is the main function and all the other functions
            %are called inside this function. (%genROIs, precluster, run 
            %all clusters, gen all statisitics)

            Ind = obj.SMD.X_SE == 0 | obj.SMD.Y_SE == 0;
            obj.SMD.X(Ind) = [];
            obj.SMD.Y(Ind) = [];
            obj.SMD.X_SE(Ind) = [];
            obj.SMD.Y_SE(Ind) = [];
            obj.SMD.FrameNum(Ind) = [];
            if isempty(obj.SMD.FrameNum)
                obj.SMD.FrameNum = zeros(size(obj.SMD.X));
            else
                obj.SMD.FrameNum = single(obj.SMD.FrameNum);
            end
            ROIs = obj.genROIs();
            if obj.PImageFlag == 1
                if isempty(obj.PImageSize) || ~isscalar(obj.PImageSize)
                    error('PImageSize must be given as a scalar');
                end
                SZ = ceil(obj.PImageSize/obj.PixelSize);
                PostIm = zeros(SZ,'single');
            end
            obj = obj.precluster(ROIs);
            if obj.ChainFlag == 1
                obj.Chain = cell(length(obj.ClusterSMD),1);
            end
            obj.MAPN.X = [];
            obj.MAPN.Y = [];
            obj.MAPN.Z = [];
            obj.MAPN.X_SE = [];
            obj.MAPN.Y_SE = [];
            obj.MAPN.Z_SE = [];
            obj.MAPN.AlphaX = [];
            obj.MAPN.AlphaY = [];
            obj.MAPN.AlphaX_SE = [];
            obj.MAPN.AlphaY_SE = [];
            obj.MAPN.Nmean = [];
            obj.MAPN.N = 0;
            ClustNumHeirar = length(obj.ClusterSMD);
            if isempty(obj.SMD.FrameNum)
               MaxAlpha = 0; 
            elseif max(obj.SMD.FrameNum~=0)
               MaxAlpha = obj.Drift / max(obj.SMD.FrameNum); 
            else
               MaxAlpha = 0; 
            end
            for nn = 1:ClustNumHeirar
                if nn/50 == ceil(nn/50)
                    fprintf('Cluster: %g out of %g\n',nn,ClustNumHeirar);
                end
                [TChain]=BaGoL.BaGoL_RJMCMC(obj.ClusterSMD(nn),obj.Lambda,MaxAlpha,obj.P_Jumps,obj.N_Trials,obj.N_Burnin,0);
                if obj.PImageFlag == 1
                    PostIm = obj.genPosterior(PostIm,SZ,TChain,ROIs,nn);
                end
                obj = obj.genMAPN(TChain,ROIs,nn);
                if obj.ChainFlag == 1
                    obj.Chain{nn} = TChain; 
                end
            end
            if obj.PImageFlag == 1
                obj.PImage = PostIm;
            end

            if obj.ChainFlag == 1
                for mm = length(obj.ClusterSMD):-1:1
                     if isempty(obj.Chain{mm})
                         obj.Chain(mm) = []; 
                         obj.ClusterSMD(mm) = [];
                     end
                end
            end
        end
        Ind = removeOverlap(obj,ROIs,MeanX,MeanY,ii)
        obj = genMAPN(obj,TChain,ROIs,nn);
        ROIs = genROIs(obj)
        [PostIm] = genPosterior(obj,PostIm,SZ,Chain,ROIs,ii)
        plotNND(obj,SaveDir)
    end
    
    methods (Static) 
       [SMD,SMC]=genCluster(StuctType,Scale,Ndist,PSF,MeanPhoton,DriftVec,PlotFlag); 
       [SMD,SMD_combined]=frameConnect(SMDin,LOS,MaxDistance,MaxFrameGap,FitType);
       saveMAPN(Directory,FileType,MAPN)
       errPlot(SMD);
       SMD=loadH5(DataDir,FileName)
       [Chain]=BaGoL_RJMCMC(SMD,Lambda,MaxAlpha,PMove,NChain,NBurnin,DEBUG)
       [SRIm,MapIm]=makeIm(SMD,MAPN,SZ,PixSize,XStart,YStart)
       [Alpha,XShift,YShift,Aligned,Chain] = align_template(Temp,Input,Start,Cutoff,NChain,PlotFlag)
       [ImageOut,Image] = scalebar(Image,PixelSize,Length,Location)
       dispIm()
    end
    
end
