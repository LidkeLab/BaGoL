function [SRIm,MapIm]=makeIm(SMD,MAPN,SZ,PixSize,XStart,YStart)
%makeIm() produces SR-Image and Map-Image.
%
% INPUTS:
%   SMD:     Structure containing X, Y, X_SE, Y_SE of SR-localizations. (nm)
%   MAPN:    Structure containing X, Y, X_SE, Y_SE of MAPN-results. (nm)
%   SZ:      Size of the processed area (image). (nm)
%   PixSize: Pixel size of the image. (nm)
%
% OUTPUTS: 
%   SRIm:    Super-resolution image.
%   MapIm:   MAPN image.
%
% Created by:
%   Mohamadreza Fazel (Lidke Lab, 2019).
%
%% SR-Image

if nargin < 5
    MinX = (SZ - (max(SMD.X)-min(SMD.X)))/2 - min(SMD.X); 
    MinY = (SZ - (max(SMD.Y)-min(SMD.Y)))/2 - min(SMD.Y);
else
    MinX = XStart;
    MinY = YStart;
end

BoxSize = floor(10*median([SMD.X_SE;SMD.Y_SE])/PixSize); %Size of box in pixels
if floor(BoxSize/2)~=BoxSize/2
   BoxSize = BoxSize + 1; 
end
XBox = floor((SMD.X+MinX)/PixSize)-floor(BoxSize/2);
YBox = floor((SMD.Y+MinY)/PixSize)-floor(BoxSize/2);
X = (SMD.X+MinX) - PixSize*floor((SMD.X+MinX)/PixSize) + PixSize*floor(BoxSize/2)+1;
Y = (SMD.Y+MinY) - PixSize*floor((SMD.Y+MinY)/PixSize) + PixSize*floor(BoxSize/2)+1;

[Xg,Yg,~]=meshgrid((PixSize/2:PixSize:BoxSize*PixSize-PixSize/2),...
    (PixSize/2:PixSize:BoxSize*PixSize-PixSize/2),(1:length(SMD.X)));
MuX = ones(size(Xg));
SigX = ones(size(Xg));
MuY = MuX;
SigY = SigX;
for ii = 1:length(SMD.X)
     SigX(:,:,ii) = SMD.X_SE(ii);
     SigY(:,:,ii) = SMD.Y_SE(ii);
     MuX(:,:,ii) = X(ii);
     MuY(:,:,ii) = Y(ii);
end
Im = normpdf(Xg,MuX,SigX).*normpdf(Yg,MuY,SigY);

SRImT = zeros(ceil(SZ/PixSize)+100,'single');
for ii = 1:length(SMD.X) 
    try
        SRImT(YBox(ii)+51:YBox(ii)+BoxSize+50,...
            XBox(ii)+51:XBox(ii)+BoxSize+50) = ...
            SRImT(YBox(ii)+51:YBox(ii)+BoxSize+50,...
            XBox(ii)+51:XBox(ii)+BoxSize+50) + Im(:,:,ii);
    catch
    end
end
SRIm = SRImT(51:end-50,51:end-50);
%% Map-Image

Ind = MAPN.X_SE > 25 | MAPN.Y_SE > 25 | MAPN.X_SE == 0 | MAPN.Y_SE == 0;
MAPN.X(Ind) = [];
MAPN.Y(Ind) = [];
MAPN.X_SE(Ind) = [];
MAPN.Y_SE(Ind) = [];

BoxSize = floor(10*median([MAPN.X_SE;MAPN.Y_SE])/PixSize); %Size of box in pixels
XBox = floor((MAPN.X+MinX)/PixSize)-floor(BoxSize/2);
YBox = floor((MAPN.Y+MinY)/PixSize)-floor(BoxSize/2);
X = (MAPN.X+MinX) - PixSize*floor((MAPN.X+MinX)/PixSize) + PixSize*floor(BoxSize/2)+1;
Y = (MAPN.Y+MinY) - PixSize*floor((MAPN.Y+MinY)/PixSize) + PixSize*floor(BoxSize/2)+1;

[Xg,Yg,~]=meshgrid((PixSize/2:PixSize:BoxSize*PixSize-PixSize/2),...
    (PixSize/2:PixSize:BoxSize*PixSize-PixSize/2),(1:length(SMD.X)));
MuX = ones(size(Xg));
SigX = ones(size(Xg));
MuY = MuX;
SigY = SigX;
for ii = 1:length(MAPN.X)
     SigX(:,:,ii) = MAPN.X_SE(ii);
     SigY(:,:,ii) = MAPN.Y_SE(ii);
     MuX(:,:,ii) = X(ii);
     MuY(:,:,ii) = Y(ii);
end
Im = normpdf(Xg,MuX,SigX).*normpdf(Yg,MuY,SigY);
MapImT = zeros(ceil(SZ/PixSize)+100,'single');
for ii = 1:length(MAPN.X)
    try
        MapImT(YBox(ii)+51:YBox(ii)+BoxSize+50,...
            XBox(ii)+51:XBox(ii)+BoxSize+50) = ...
            MapImT(YBox(ii)+51:YBox(ii)+BoxSize+50,...
            XBox(ii)+51:XBox(ii)+BoxSize+50) + Im(:,:,ii);
    catch
    end
end
MapIm = MapImT(51:end-50,51:end-50);
end