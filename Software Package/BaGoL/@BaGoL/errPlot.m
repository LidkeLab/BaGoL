function errPlot(SMD) 
%plot of the localizations accuracies, where circles represents the errors
%(circle radius = 2*error).

Ang = 0:0.1:2*pi+0.25;
XCir = zeros(length(Ang),length(SMD.X),'single');
YCir = XCir;
for nn = 1:length(SMD.X)
    Prec = 2*sqrt((SMD.X_SE(nn)^2+SMD.Y_SE(nn)^2)/2);
    XCir(:,nn) = Prec*cos(Ang')+SMD.X(nn);
    YCir(:,nn) = Prec*sin(Ang')+SMD.Y(nn);
end
figure;
plot(XCir,YCir,'m')
axis equal
end