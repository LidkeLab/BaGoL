function [Alpha,XShift,YShift,Aligned,Chain] = align_template(Temp,Input,Start,Cutoff,NChain,PlotFlag)
%align_template() aligns the given structure with the given template.
% This function uses nearest neighbor and MCMC to find the rotation and
% transformation to align the given structure with the given template.
%
% INPUTS:
%   Temp:    Structure containing the template structure coordinates.
%            Temp.X and Temp.Y
%   Input:   Structure containing the structure to be aligned with template
%            Input.X and Input.Y
%   Start:   Structure containing the starting values for the rotation
%            angle and translations in X and Y. 
%            Start.Theta, Start.ShiftX, Start.ShiftY
%   Cutoff:  If there is any missing point in the given structure then it
%            will be penalized by this value in the cost function. 
%   NChain:  Number of jumps in MCMC algorithm. (Default = 5000)
%   PlotFlag: If 1 shows an animation of the accepted jumps. (Default = 0)
%
% OUTPUT:
%   Alpha:   The rotation angle in radian.
%   XShift:  The translation along X-axis.
%   YShift:  The translation along Y-axis.
%   Aligned: Structure containing the aligned coordinates. 
%            Aligned.X & Aligned.Y
%   Chain:   Structure containing  
%
% Citation:
%   Mohamadreza Fazel, Lidke Lab 2019.
if nargin < 4
   error('There must be at leat 4 inputs.') 
end

if nargin < 5
   NChain = 5000;
   PlotFlag = 0;
end

if nargin < 6
   PlotFlag = 0; 
end

Aligned.X = [];
Aligned.Y = [];

TempX = Temp.X - mean(Temp.X);
TempY = Temp.Y - mean(Temp.Y);

XLim = [min(TempX)-10 max(TempX)+10];
YLim = [min(TempY)-10 max(TempY)+10];

X = Input.X';
Y = Input.Y';
X = X - mean(X);
Y = Y - mean(Y);

Points = [X;Y];
[~,Dis] = knnsearch(Points',[TempX',TempY']);
Dis = abs(Dis);
Dis(Dis>Cutoff) = Cutoff;
Cost_Current = sum(Dis);

Start_Theta = Start.Theta;
Start_ShiftX = Start.ShiftX;
Start_ShiftY = Start.ShiftY;

Theta = zeros(NChain,1);
DelX = zeros(NChain,1);
DelY = zeros(NChain,1);

Theta(1) = Start_Theta;
DelX(1) = Start_ShiftX;
DelY(1) = Start_ShiftY;
if PlotFlag
    figure;hold;
    V1 = VideoWriter('Alignment','MPEG-4');
    V1.FrameRate = 10; 
    open(V1) 
end
for nn = 2:NChain
    
    if nn < NChain/2
        Theta_P = Theta(nn-1)+randn();
        DelX_P = DelX(nn-1) + 0.5*randn();
        DelY_P = DelY(nn-1) + 0.5*randn();
    else
        Theta_P = Theta(nn-1)+0.1*randn();
        DelX_P = DelX(nn-1) + 0.02*randn();
        DelY_P = DelY(nn-1) + 0.02*randn();
    end
    R = [cos(Theta_P),sin(Theta_P);-sin(Theta_P),cos(Theta_P)];
    
    Points = [X;Y];
    Points = R*Points;
    Points(1,:) = Points(1,:)+DelX_P;
    Points(2,:) = Points(2,:)+DelY_P;
    
    [~,Dis]=knnsearch(Points',[TempX',TempY']);
    Dis = abs(Dis);
    Dis(Dis>Cutoff) = 20;
    Cost_Proposed = sum(Dis);
    
    if Cost_Current - Cost_Proposed > -rand()
        if PlotFlag
            plot(TempX,TempY,'*')
            hold;
            plot(Points(1,:),Points(2,:),'*')
            title(sprintf('Step:%g',nn))
            xlim(XLim);ylim(YLim)
            %pause(0.25)
            SP = sprintf('Data: %g',nn);
            legend(SP)
            Frame = getframe(gcf);
            writeVideo(V1,Frame);
            hold off
        end
        Theta(nn) = Theta_P;
        DelX(nn) = DelX_P;
        DelY(nn) = DelY_P;
        Cost_Current = Cost_Proposed;
    else
        Theta(nn) = Theta(nn-1);
        DelX(nn) = DelX(nn-1);
        DelY(nn) = DelY(nn-1);
    end
    
end
if PlotFlag
    close(V1)
end
Chain.Theta = Theta;
Chain.XShift = DelX;
Chain.YShift = DelY;

Alpha = mean(Theta(round(NChain/2),:));
XShift = mean(DelX(round(NChain/2),:));
YShift = mean(DelY(round(NChain/2),:));

R = [cos(Alpha),sin(Alpha);-sin(Alpha),cos(Alpha)];
Points = [X;Y];
Points = R*Points;
Points(1,:) = Points(1,:)+XShift;
Points(2,:) = Points(2,:)+YShift;

Aligned.X = cat(1,Aligned.X,Points(1,:)'+mean(Temp.X));
Aligned.Y = cat(1,Aligned.Y,Points(2,:)'+mean(Temp.Y));

end