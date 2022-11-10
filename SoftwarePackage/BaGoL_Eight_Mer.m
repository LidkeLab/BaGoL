%% This is a demo of the RJMCMC step in BaGoL using a simulated 8-mer 
%  The localization per emitter distribution is known and has a 
%  Poisson distribution parameterized by lambda.  
%  
%  The core RJMCMC method is used to analyze the data as a single cluster.
%
%  This script shows a movie of RJMCMC chain exploring possible numbers of
%  emitters and thier locataions. The circles show input localizations and 
%  their colors represent the allocation of localizations to emitters. The  
%  color changes indicate assignment uncertainties. The thick red circles  
%  show the found emitters. The chain is used to make final BaGoL Posterior 
%  and MAPN images. These images are saved under the 'Results_8mer' folder. 
%  
%  We have used a short chain here (1000 jumps in total) to demonstrate the 
%  chain and different jump types used in BaGoL. However, to guarantee the 
%  chain convergence, we suggest a longer chain, say NChain=2000 and 
%  NBurnin = 1000. To speed up the code, you can also set PlotFlag = 0 to  
%  stop displaying the movie.
%

%% RJMCMC Parameter Setup

ImSize = 60; %(nm)
Lambda = 20; %Average number of localizations per emitter (known from simulated data)
MaxAlpha = 0; %Emitter Drift per frame
PMove = [0.25 0.25 0.25 0.25]; %Probability of moves
NChain = 500; %Number of post-burn-in jumps
NBurnin = 500; %Number of burn-in jumps
%PlotFlag can be a) 0: no animation, b) 1: animation with small circles 
%c) 2: animation with circles with radii equal to 2 times of localization
%precisions.
PlotFlag = 1; %Display animation of jumps

%Make Save Directory
SaveDir = fullfile('Results_8mer'); 
if ~isdir(SaveDir)
    mkdir(SaveDir); 
end

%% Load Data and run the core RJMCMC method

DataDir = fullfile('Data');
FileName = 'SMD_Sim_8mer.mat';
load(fullfile(DataDir,FileName))

Chain = BaGoL.BaGoL_RJMCMC(SMD,Lambda,MaxAlpha,PMove,NChain,NBurnin,PlotFlag);

%% Extracting locations from the RJMCMC chain and generating MAPN


%Extracting the chain of emitter positions
X=[];
Y=[];
N=[];
for ii = 1:length(Chain)
    X=cat(1,Chain(ii).X,X);
    Y=cat(1,Chain(ii).Y,Y);
    N=cat(1,length(Chain(ii).X),N);
end

%Generating MAPN
MAPN.X = [];
MAPN.Y = [];
MAPN.X_SE = [];
MAPN.Y_SE = [];
ModeN = mode(N); %Most repeated model
XC = X(N==ModeN);
YC = Y(N==ModeN);
ID = kmeans([XC,YC],ModeN,'Replicates',5,'MaxIter',200);
for ii = 1:ModeN
    MAPN.X = cat(1,mean(XC(ID==ii)),MAPN.X);
    MAPN.Y = cat(1,mean(YC(ID==ii)),MAPN.Y);
    MAPN.X_SE = cat(1,std(XC(ID==ii)),MAPN.X_SE);
    MAPN.Y_SE = cat(1,std(YC(ID==ii)),MAPN.Y_SE);
end

%% Producing images

%Producing Posterior image from the chain
PixSize = 0.5; %PixelSize (nm)
XInd = floor(X/PixSize);
YInd = floor(Y/PixSize);
IndRemove = XInd>ImSize/PixSize | XInd < 1 | YInd>ImSize/PixSize | YInd < 1;
XInd(IndRemove) = [];
YInd(IndRemove) = [];
Posterior = zeros(ImSize/PixSize,'single');
Ind = sub2ind([ImSize/PixSize,ImSize/PixSize],YInd,XInd);
[Counts,Pixel]=hist(Ind,unique(Ind));
IndRemove = Counts == 0;
Counts(IndRemove) = [];
Pixel(IndRemove) = [];
Posterior(Pixel) = Posterior(Pixel) + Counts';

%Posterior Image
Posterior=BaGoL.scaleIm(Posterior,90);
Posterior=BaGoL.scalebar(Posterior,PixSize,5);
imwrite(Posterior,hot(256),fullfile(SaveDir,'Posterior-SRImage.png'));

%MAPN and SR images
XStart = 0;
YStart = 0;
[SRIm]=BaGoL.makeIm(SMD,ImSize,PixSize,XStart,YStart);
[MapIm]=BaGoL.makeIm(MAPN,ImSize,PixSize,XStart,YStart);

MapIm=BaGoL.scaleIm(MapIm,95);
MapIm=BaGoL.scalebar(MapIm,PixSize,5);
imwrite(MapIm,hot(256),fullfile(SaveDir,'MAPN-SRImage.png'));

SRIm=BaGoL.scaleIm(SRIm,99);
SRIm=BaGoL.scalebar(SRIm,PixSize,5);
imwrite(SRIm,hot(256),fullfile(SaveDir,'PreBaGoL-SRImage.png'));



