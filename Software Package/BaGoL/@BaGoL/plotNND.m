function plotNND(obj,SaveDir)
[~,Dis]=knnsearch([obj.MAPN.X,obj.MAPN.Y],[obj.MAPN.X,obj.MAPN.Y],'k',2);
Dis = Dis(:,2);
P = prctile(Dis,99);
figure;hist(Dis(Dis<P),30)
xlabel('NND(nm)','FontSize',15);ylabel('Frequency','FontSize',15)
if nargin > 1
   print(gcf,fullfile(SaveDir,'NND'),'-dpng'); 
   Dis = double(Dis);
   if size(Dis,2) > 1
      Dis = Dis';
    end
    save(fullfile(SaveDir,'NND.txt'),'Dis','-ascii')
end
end