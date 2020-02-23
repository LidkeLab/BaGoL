function [PostIm] = genPosterior(obj,PostIm,SZ,Chain,ROIs,ii)
%making the posterior image of subregions using the returned 
%chain for each cluster
XChain = [];
YChain = [];

for nn = 1:length(Chain)
    if ~isempty(Chain(nn).X)
        XChain = cat(1,XChain,Chain(nn).X);
        YChain = cat(1,YChain,Chain(nn).Y);
    end
end
Ind = obj.removeOverlap(ROIs,XChain,YChain,ii);

MinX = (obj.PImageSize - (max(obj.SMD.X)-min(obj.SMD.X)))/2 - min(obj.SMD.X); 
MinY = (obj.PImageSize - (max(obj.SMD.Y)-min(obj.SMD.Y)))/2 - min(obj.SMD.Y);
YInd = floor((YChain(Ind)+MinY)/obj.PixelSize);
XInd = floor((XChain(Ind)+MinX)/obj.PixelSize);

IndRemove = XInd>SZ | XInd < 1 | YInd>SZ | YInd < 1;
XInd(IndRemove) = [];
YInd(IndRemove) = [];

LinInd = sub2ind([SZ,SZ],YInd,XInd);
[Counts,Pixel]=hist(LinInd,unique(LinInd));
IndRemove = Counts == 0;
Counts(IndRemove) = [];
Pixel(IndRemove) = [];
PostIm(Pixel) = PostIm(Pixel) + Counts';

end