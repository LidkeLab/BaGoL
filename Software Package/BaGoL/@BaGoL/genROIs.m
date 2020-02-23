function ROIs = genROIs(obj)
%data segmentation
XPix = max(obj.SMD.X);
YPix = max(obj.SMD.Y);
MinX = min(obj.SMD.X);
MinY = min(obj.SMD.Y);
if MinX < obj.ROIsize && MinX > 0
    MinX = 0;
end
if MinY < obj.ROIsize && MinY > 0
    MinY = 0; 
end
Xsub = ceil((XPix-MinX)/obj.ROIsize); %number of the sub-regions along the X-axis
Ysub = ceil((YPix-MinY)/obj.ROIsize); %number of the sub-regions along the Y-axis
ROIs = [];
for ix = 1:Xsub
    for iy = 1:Ysub
        ROIex = obj.ROIsize + 2*obj.Overlap;
        ROIey = obj.ROIsize + 2*obj.Overlap;
        sy = obj.ROIsize*(iy-1)-obj.Overlap;
        ey = obj.ROIsize*iy+obj.Overlap;
        sx = obj.ROIsize*(ix-1)-obj.Overlap;
        ex = obj.ROIsize*ix+obj.Overlap;
        if iy == 1 
            ROIey = obj.ROIsize + obj.Overlap;
            sy = 0;
            ey = obj.ROIsize + obj.Overlap;
        end
        if ix == 1
            ROIex = obj.ROIsize + obj.Overlap;
            sx = 0;
            ex = obj.ROIsize + obj.Overlap;
        end
        if iy == Ysub
            %ey = YPix;
            ROIey = YPix-obj.ROIsize*(Ysub-1)+obj.Overlap;
        end
        if ix == Xsub
            %ex = XPix;
            ROIex = XPix-obj.ROIsize*(Xsub-1)+obj.Overlap;
        end
        if Xsub == 1
            ROIex = XPix;
        end
        if Ysub == 1
            ROIey = YPix;
        end
        Ind = obj.SMD.X-MinX >= sx & obj.SMD.X-MinX <= ex ...
             & obj.SMD.Y-MinY >= sy & obj.SMD.Y-MinY <= ey;
        ROIs(iy,ix).X = obj.SMD.X(Ind);
        ROIs(iy,ix).Y = obj.SMD.Y(Ind);
        if isfield(obj.SMD,'Z') && isempty(obj.SMD.Z)
            ROIs(iy,ix).Z = [];
        else
            ROIs(iy,ix).Z = obj.SMD.Z(Ind);
        end
        ROIs(iy,ix).X_SE = obj.SMD.X_SE(Ind);
        ROIs(iy,ix).Y_SE = obj.SMD.Y_SE(Ind);
        if isempty(obj.SMD.Z)
            ROIs(iy,ix).Z_SE = [];
        else
            ROIs(iy,ix).Z_SE = obj.SMD.Z_SE(Ind);
        end
        if ~isempty(obj.SMD.FrameNum)
            ROIs(iy,ix).FrameNum = obj.SMD.FrameNum(Ind);
        else
            ROIs(iy,ix).FrameNum = []; 
        end
        ROIs(iy,ix).Xsize = ROIex;
        ROIs(iy,ix).Ysize = ROIey;
    end
end
end