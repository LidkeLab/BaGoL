function obj = genMAPN(obj,Chain,ROIs,ii)
    %extracting the most repeated model from the chain, using the
    %moset repeated model to calculate MAPN results.
    TChain = Chain;
    N = [];
    XChain = [];
    YChain = [];
    AlphaX_Chain = [];
    AlphaY_Chain = [];
    NChain = [];
    for nn = 1:length(TChain)
        if TChain(nn).N == 0
            continue;
        end
        N = cat(1,N,length(TChain(nn).X));
    end
    Frequency = hist(N,(0:max(N)));
    Frequency(1) = [];
    MostFrequent = find(Frequency==max(Frequency));
    MostFrequent = MostFrequent(1);
    obj.MAPN.N = obj.MAPN.N+MostFrequent;
    for nn = 1:length(TChain)
        if TChain(nn).N == MostFrequent
            NSum = zeros(MostFrequent,1);
            XChain = cat(1,XChain,TChain(nn).X);
            YChain = cat(1,YChain,TChain(nn).Y);
            AlphaX_Chain = cat(1,AlphaX_Chain,TChain(nn).AlphaX);
            AlphaY_Chain = cat(1,AlphaY_Chain,TChain(nn).AlphaY);
            for pp = 1: MostFrequent(1)
                NSum(pp) = sum(single(TChain(nn).ID)==pp); 
            end
            NChain = cat(1,NChain,NSum);
        end
    end
    Points(:,1)=XChain;
    Points(:,2)=YChain;
    if ~isempty(obj.SMD.Z)
        Points(:,3)=ZChain;                     
    end
    ID = kmeans(Points,MostFrequent,'Replicates',5,'MaxIter',150);
    Map.X = [];
    Map.Y = [];
    Map.X_SE = [];
    Map.Y_SE = [];
    Map.AlphaX = [];
    Map.AlphaY = [];
    Map.AlphaX_SE = [];
    Map.AlphaY_SE = [];
    Map.Nmean = [];
    for nn = 1:MostFrequent
        TX = mean(XChain(ID==nn)); 
        TY = mean(YChain(ID==nn));
        TX_SE = std(XChain(ID==nn));
        TY_SE = std(YChain(ID==nn));
        TAlphaX = mean(AlphaX_Chain(ID==nn));
        TAlphaY = mean(AlphaY_Chain(ID==nn));
        TAlphaX_SE = std(AlphaX_Chain(ID==nn));
        TAlphaY_SE = std(AlphaY_Chain(ID==nn));
        TNmean = mean(NChain(ID==nn));
        Map.X = cat(1,Map.X,TX);
        Map.Y = cat(1,Map.Y,TY);
        Map.X_SE = cat(1,Map.X_SE,TX_SE);
        Map.Y_SE = cat(1,Map.Y_SE,TY_SE);
        Map.AlphaX = cat(1,Map.AlphaX,TAlphaX);
        Map.AlphaY = cat(1,Map.AlphaY,TAlphaY);
        Map.AlphaX_SE = cat(1,Map.AlphaX_SE,TAlphaX_SE);
        Map.AlphaY_SE = cat(1,Map.AlphaY_SE,TAlphaY_SE);
        Map.Nmean = cat(1,Map.Nmean,TNmean);
    end
    Ind = obj.removeOverlap(ROIs,Map.X,Map.Y,ii);
    obj.MAPN.X = cat(1,obj.MAPN.X,Map.X(Ind));
    obj.MAPN.Y = cat(1,obj.MAPN.Y,Map.Y(Ind));
    obj.MAPN.X_SE = cat(1,obj.MAPN.X_SE,Map.X_SE(Ind));
    obj.MAPN.Y_SE = cat(1,obj.MAPN.Y_SE,Map.Y_SE(Ind));
    obj.MAPN.AlphaX = cat(1,obj.MAPN.AlphaX,Map.AlphaX(Ind));
    obj.MAPN.AlphaY = cat(1,obj.MAPN.AlphaY,Map.AlphaY(Ind));
    obj.MAPN.AlphaX_SE = cat(1,obj.MAPN.AlphaX_SE,Map.AlphaX_SE(Ind));
    obj.MAPN.AlphaY_SE = cat(1,obj.MAPN.AlphaY_SE,Map.AlphaY_SE(Ind));
    if ~isempty(obj.SMD.Z)
        obj.MAPN.Z = cat(1,obj.MAPN.Z,Map.Z(Ind));
        obj.MAPN.Z_SE = cat(1,obj.MAPN.Z_SE,Map.Z_SE(Ind));
    end
    obj.MAPN.Nmean = cat(1,obj.MAPN.Nmean,Map.Nmean(Ind));
end
