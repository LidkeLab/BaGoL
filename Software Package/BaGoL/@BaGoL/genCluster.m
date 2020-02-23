function [SMD,SMC]=genCluster(StuctType,Scale,Ndist,PSF,MeanPhoton,DriftVec,PlotFlag) 
%genCluster creates localization data for a patterned cluster.
%   Generates cluster data in the form used by clusterBC()
%   SMD=genCluster(StuctType,Scale,Lambda,DriftVec) 
% INPUTS:
%   StuctType:    '?grid','?mer','user'
%   Scale:         Spatial scale (nm) 
%       '?Grid':   Smallest NND
%       '?mer':    Radius
%       'user':    Ignored (uses Origami distances)
%   Ndist:         This can be either a scalar or a vector. Scalar means 
%                  the parameter for Poisson dist., and vector means array
%                  of distribution starting from zero and increment by one.
%   DriftVec:      Max linear drift [Dx1] where D is dimensionality of data
% OUTPUT:
%   SMD:    SMD structure
%       SMD.X
%       SMD.X_SE
%       SMD.Y
%       SMD.Y_SE
%       SMD.Z (optional)
%       SMD.Z_SE (optional)
%       SMD.FrameNum: time stamp of localization (used for drift)
%   SMC:    SMC structure
%       SMC.X     Cluster X-center
%       SMC.Y     Cluster Y-center
%       SMC.Z     Cluster Z-center
%
% Created by:
%   Mohamadreza Fazel, (Lidke lab 2019)   
%
if isscalar(Ndist)
    Ndist = poisspdf((0:1000),Ndist); 
end
if ~exist('DriftVec','var')
   DriftVec = zeros(1,1000); 
end
if ~exist('PlotFlag','var')
   PlotFlag = 0;  
end

STR = StuctType(end-2:end);
MaxFrame = length(DriftVec);
Nmean = find(Ndist==max(Ndist));
Tau = MaxFrame/Nmean(1);
SMD.X = [];
SMD.Y = [];
SMD.Z = [];
SMD.X_SE = [];
SMD.Y_SE = [];
SMD.Z_SE = [];
SMD.FrameNum = [];
SMC.X = [];
SMC.Y = [];
if strcmp(STR,'rid')
    N = str2double(StuctType(1))*str2double(StuctType(3));
    XC = Scale:Scale:str2double(StuctType(1))*Scale;
    YC = Scale:Scale:str2double(StuctType(3))*Scale;
    SMC.X = repmat(XC,1,length(YC));
    for mm = 1:length(YC)
        SMC.Y = cat(2,SMC.Y,YC(mm)*ones(1,length(XC)));
    end
    SMC.Z = [];
    for nn = 1:N
        NThisEmitter = genRand(Ndist); 
        FrameNum = randperm(MaxFrame,NThisEmitter);
        %FrameNum = round(cumsum(exprnd(Tau,[NThisEmitter,1])));
        %FrameNum(FrameNum>MaxFrame | FrameNum == 0) = [];
        for mm = 1:length(FrameNum)
            Photons = exprnd(MeanPhoton); 
            Prec = PSF/sqrt(Photons);
            SMD.X = cat(1,SMD.X,SMC.X(nn)+Prec*randn()+DriftVec(FrameNum(mm),1));
            SMD.Y = cat(1,SMD.Y,SMC.Y(nn)+Prec*randn()+DriftVec(FrameNum(mm),2));
            SMD.X_SE = cat(1,SMD.X_SE,Prec);
            SMD.Y_SE = cat(1,SMD.Y_SE,Prec);
            SMD.FrameNum = cat(1,SMD.FrameNum,FrameNum);
        end
    end
    Ind = SMD.X_SE < 40 & SMD.Y_SE < 40;
    SMD.X = SMD.X(Ind);
    SMD.Y = SMD.Y(Ind);
    SMD.X_SE = SMD.X_SE(Ind);
    SMD.Y_SE = SMD.Y_SE(Ind);
    SMD.FrameNum = SMD.FrameNum(Ind);
elseif strcmp(STR,'mer') || strcmp(STR,'Mer')
    N = str2double(StuctType(1:end-3));
    Theta = linspace(0,2*pi,N+1)';
    Theta(end) = [];
    Theta = Theta + (Theta(2)-Theta(1))/2;
    R = Scale*ones(N,1);
    [XC,YC] = pol2cart(Theta,R);
    SMC.X = XC + 2*Scale;
    SMC.Y = YC + 2*Scale;
    SMC.Z = [];
    for nn = 1:N
        NThisEmitter = genRand(Ndist);
        FrameNum = randperm(MaxFrame,NThisEmitter);
        %FrameNum = round(cumsum(exprnd(Tau,[NThisEmitter,1])));
        %FrameNum(FrameNum > MaxFrame | FrameNum == 0) = [];
        for mm = 1:length(FrameNum)
            Photons = exprnd(MeanPhoton); 
            Prec = PSF/sqrt(Photons);
            try
                SMD.X = cat(1,SMD.X,SMC.X(nn)+Prec*randn()+DriftVec(FrameNum(mm),1));
                SMD.Y = cat(1,SMD.Y,SMC.Y(nn)+Prec*randn()+DriftVec(FrameNum(mm),2));
            catch
               a = 0; 
            end
            SMD.X_SE = cat(1,SMD.X_SE,Prec);
            SMD.Y_SE = cat(1,SMD.Y_SE,Prec);
            SMD.FrameNum = cat(1,SMD.FrameNum,FrameNum(mm)); 
        end
    end
    Ind = SMD.X_SE < 40 & SMD.Y_SE < 40;
    SMD.X = SMD.X(Ind);
    SMD.Y = SMD.Y(Ind);
    SMD.X_SE = SMD.X_SE(Ind);
    SMD.Y_SE = SMD.Y_SE(Ind);
    SMD.FrameNum = SMD.FrameNum(Ind);
elseif strcmp(STR,'ser') 
    SMC.X = Scale.X;
    SMC.Y = Scale.Y;
    SMC.Z = [];
    for nn = 1:length(SMC.X)
        NThisEmitter = genRand(Ndist);
        FrameNum = randperm(MaxFrame,NThisEmitter);
        %FrameNum = round(cumsum(exprnd(Tau,[NThisEmitter,1])));
        %FrameNum(FrameNum>MaxFrame | FrameNum == 0) = [];
        for mm = 1:length(FrameNum)
            Photons = exprnd(MeanPhoton); 
            Prec = PSF/sqrt(Photons);
            SMD.X = cat(1,SMD.X,SMC.X(nn)+Prec*randn())+DriftVec(FrameNum(mm));
            SMD.Y = cat(1,SMD.Y,SMC.Y(nn)+Prec*randn())+DriftVec(FrameNum(mm));
            SMD.X_SE = cat(1,SMD.X_SE,Prec);
            SMD.Y_SE = cat(1,SMD.Y_SE,Prec);
            SMD.FrameNum = cat(1,SMD.FrameNum,FrameNum(mm)); 
        end
    end
    Ind = SMD.X_SE < 40 & SMD.Y_SE < 40;
    SMD.X = SMD.X(Ind);
    SMD.Y = SMD.Y(Ind);
    SMD.X_SE = SMD.X_SE(Ind);
    SMD.Y_SE = SMD.Y_SE(Ind);
    SMD.FrameNum = SMD.FrameNum(Ind);
end

if PlotFlag
    Ang = 0:0.05:2*pi+0.05;
    for nn = 1:length(SMD.X)
        Prec = 2*sqrt((SMD.X_SE(nn)^2+SMD.Y_SE(nn)^2)/2);
        XCir = Prec*cos(Ang)+SMD.X(nn);
        YCir = Prec*sin(Ang)+SMD.Y(nn);
        if nn == 1
            figure; 
        end
        plot(XCir,YCir,'color',[0 0.4470 0.7410])
        axis equal
        if nn == 1
            hold; 
        end
    end
    plot(SMC.X,SMC.Y,'ok','linewidth',1.5)
end

end

function Out=genRand(Ndist)
CDFarray = cumsum(Ndist);
RandNum = rand();
IndArray = find(CDFarray-RandNum >= 0);
Out = IndArray(1);
end

