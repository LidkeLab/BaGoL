function plotCluster(obj,FlagOrDir,FlagPlot)
%Make plots of clusters and data-points.
% The first plot is color-coded data points and the cluster centers.
% The second plot is plot of data points before and after filtering and the
% cluster centers.
FlagSave = 0;
if nargin < 2
   FlagPlot = 3; 
else
    if ischar(FlagOrDir)
        SaveDir = FlagOrDir;
        FlagSave = 1;
    else
        FlagPlot = FlagOrDir; 
    end
end
if ~exist('FlagPlot')
    FlagPlot = 3; 
end

if FlagPlot == 1 || FlagPlot == 3
    figure;hold;
    P1 = scatter(obj.ClusterSMD(1).X,obj.ClusterSMD(1).Y,[],obj.ClusterSMD(1).FrameNum);
    for mm = 2:length(obj.ClusterSMD)
        scatter(obj.ClusterSMD(mm).X,obj.ClusterSMD(mm).Y,[],obj.ClusterSMD(mm).FrameNum)
    end
    P2 = plot(obj.MAPN.X,obj.MAPN.Y,'*k');
    colorbar();xlabel('X(nm)','FontSize',12);ylabel('Y(nm)','FontSize',12)
    axis equal;
    legend([P1,P2],{'Localizations','Groups centers'})    
    if FlagSave
       saveas(gcf,fullfile(SaveDir,'ScatterPlot.fig')) 
    end
end
if FlagPlot == 2 || FlagPlot ==3
    figure;hold;
    P1 = plot(obj.SMD.X,obj.SMD.Y,'.k');
    P2 = plot(obj.ClusterSMD(1).X,obj.ClusterSMD(1).Y,'.','color',[0,0.447,0.741]);
    for mm = 2:length(obj.ClusterSMD)
        plot(obj.ClusterSMD(mm).X,obj.ClusterSMD(mm).Y,'.','color',[0,0.447,0.741]) 
    end
    P3 = plot(obj.MAPN.X,obj.MAPN.Y,'*','color',[0.3,0.15,0.15],'linewidth',1.25);
    axis equal;xlabel('X(nm)','FontSize',12);ylabel('Y(nm)','FontSize',12)
    legend([P1,P2,P3],{'Before filtering','After filtering','Groups centers'})
    if FlagSave
       saveas(gcf,fullfile(SaveDir,'Signal-BG-Cluster.fig')) 
    end
end
end