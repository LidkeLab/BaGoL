function Ind = removeOverlap(obj,ROIs,MeanX,MeanY,ii)
%remove the clusters with the centers in the overlapping region
ix = obj.ClusterSMD(ii).ix;
iy = obj.ClusterSMD(ii).iy;
XPix = max(obj.SMD.X);
YPix = max(obj.SMD.Y);
[Ysub,Xsub] = size(ROIs); 
MinX = min(obj.SMD.X);
MinY = min(obj.SMD.Y);
if MinX < obj.ROIsize
    MinX = 0;
end
if MinY < obj.ROIsize
    MinY = 0; 
end
sy = obj.ROIsize*(iy-1);
ey = obj.ROIsize*iy;
sx = obj.ROIsize*(ix-1);
ex = obj.ROIsize*ix;
if iy == Ysub
    ey = YPix;
end
if ix == Xsub
    ex = XPix;
end                    
Ind = find(MeanX - MinX >= sx & MeanX - MinX < ex ...
      & MeanY - MinY >= sy & MeanY - MinY < ey);
end