function obj=precluster(obj,ROIs)
%removing outliers using nearest neighbors distances and 
%preclustering using heirarchical clustering
obj.ClusterSMD = [];
[Ys,Xs] = size(ROIs);
for iy = 1:Ys
    for ix = 1:Xs
        [y,~]=size(ROIs(iy,ix).X);
        if y == 1
            XYT = [ROIs(iy,ix).X',ROIs(iy,ix).Y'];
        else
            XYT = [ROIs(iy,ix).X,ROIs(iy,ix).Y];
        end  
        
        if isempty(XYT)
           continue; 
        end
        
        if ~isempty(obj.NNN)
            NNNn = obj.NNN+1;
            [~,Dis] = knnsearch(XYT,XYT,'k',NNNn);
            Ind = Dis(:,end)<obj.NNR; 
            Coords = XYT(Ind,1:2);
        else
            Ind = ~isnan(XYT(:,1));
            Coords = single(XYT(:,1:2));
        end
        
        %clustering - give back array of SMD - one per cluster
        % There are three functions for heirarchical clustering that we have to 
        % use them in turn. The first one is pdist which calculates the 
        % distances between every two pairs of the data points. This is the 
        % input to the next function which is linkage. The out put of this 
        % function the distance between each two data points which would be the
        % input of the next function which is linkage.
        if size(Coords,1) <= 1
            continue; 
        end
        try
            Dist = pdist(Coords,'euclid');
            TC = linkage(Dist,'single');
            IndClust = cluster(TC,'Cutoff',obj.Cutoff,'Criterion','distance');
        catch
            IndClust = ones(size(Dist(:,1)),'single');
        end
                    
        % cluster is the third function that we need to use. The input to this
        % function is the output of the previous one. This function has several
        % methods that we need to decide which one is suitable for our purpose. 
        % For example the method distance cluster the data based on the distance
        % between the data points.
        NClust = max(IndClust); 
        MeanX = zeros(NClust,1);
        MeanY = zeros(NClust,1);
        for nn = 1:NClust
            MeanX(nn) = mean(Coords(IndClust==nn,1));
            MeanY(nn) = mean(Coords(IndClust==nn,2));
        end
        Len = length(obj.ClusterSMD);
        TX_SE = ROIs(iy,ix).X_SE(Ind);
        TY_SE = ROIs(iy,ix).Y_SE(Ind);
        if ~isempty(ROIs(iy,ix).FrameNum)
            TFrame = ROIs(iy,ix).FrameNum(Ind);
        end
        for nn = 1:max(IndClust)
            Indtrue = IndClust==nn;
            obj.ClusterSMD(Len+nn).X = Coords(Indtrue,1);
            obj.ClusterSMD(Len+nn).Y = Coords(Indtrue,2);
            obj.ClusterSMD(Len+nn).X_SE = TX_SE(Indtrue);
            obj.ClusterSMD(Len+nn).Y_SE = TY_SE(Indtrue);
            if isempty(ROIs(iy,ix).Z)
                obj.ClusterSMD(Len+nn).Z = [];
                obj.ClusterSMD(Len+nn).Z_SE = [];
            else
                obj.ClusterSMD(Len+nn).Z = ROIs(iy,ix).Z(Indtrue);
                obj.ClusterSMD(Len+nn).Z_SE = ROIs(iy,ix).Z_SE(Indtrue);
            end
                if ~isempty(ROIs(iy,ix).FrameNum)
                    obj.ClusterSMD(Len+nn).FrameNum = TFrame(Indtrue);
                else
                    obj.ClusterSMD(Len+nn).FrameNum = [];
                end
                obj.ClusterSMD(Len+nn).ix = ix;
                obj.ClusterSMD(Len+nn).iy = iy;
        end
    end
end

end