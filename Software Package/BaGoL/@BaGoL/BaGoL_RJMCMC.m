function [Chain]=BaGoL_RJMCMC(SMD,Lambda,MaxAlpha,PMove,NChain,NBurnin,DEBUG)
%BaGoL_RJMCMC The core RJMCMC algoritm

%DEBUG=0;

if nargin<3
    MaxAlpha=0;
end

if nargin<4
    PMove = [.3 .3 .2 .2]; %PMove = [Theta Z Birth Death]
end

if nargin<5
    NChain = 1e3; %Total
end

if nargin<6
    NBurnin = 1e3;
end

%Storage of Chain
N=length(SMD.X);
KChain=[];
ZChain=[];
Xk=cell(N,1); 
Yk=cell(N,1);
AXk=cell(N,1);
AYk=cell(N,1);
for ii=1:N
    Xk{ii}=[];Yk{ii}=[];
    AYk{ii}=[];AYk{ii}=[];
end

%Intial K Guess
K=ceil(N/prod(Lambda));

Wx = max(SMD.X)-min(SMD.X);
Wy = max(SMD.Y)-min(SMD.Y);

%These are the Mu and Alpha prior bounds
MinX=min(SMD.X);
MaxX=max(SMD.X);
MinY=min(SMD.Y);
MaxY=max(SMD.Y);

%Initial Locations
Mu_X =normrnd(mean(SMD.X),std(SMD.X),[1 K]);
Mu_Y = normrnd(mean(SMD.Y),std(SMD.Y),[1 K]);

Mu_X =SMD.X(randi(N,[1 K]))';
Mu_Y =SMD.Y(randi(N,[1 K]))';

%Initial Alphas
Alpha_X = zeros([1 K]);
Alpha_Y = zeros([1 K]);

if N < 100
   LengN = 100; 
else
   LengN = N; 
end
%Calculating the Prior
if length(Lambda)>1
   Gamma_K=Lambda(1);
   Gamma_Theta=Lambda(2);
   Pk=gampdf(N,(1:LengN)*Gamma_K,Gamma_Theta);
else
   Pk=poisspdf(N,(1:LengN)*Lambda); 
end
Pk = Pk/sum(Pk);
% Pk=poisspdf(N,(1:N)*Lambda);
% Pk=Pk/sum(Pk);

% Intial Allocation
Z=Gibbs_Z(SMD,K,Mu_X,Mu_Y,Alpha_X,Alpha_Y);

% V = VideoWriter('BaGoL_TUD','MPEG-4');
% open(V)
% Run Chain
for nn=1:NChain+NBurnin
    %Get move type:
    JumpType=length(PMove)+1-sum(rand<cumsum(PMove));
    K = length(Mu_X);
    for ii = K:-1:1
        if sum(Z==ii)==0
            Mu_X(ii)=[];
            Mu_Y(ii)=[];
            Alpha_X(ii)=[];
            Alpha_Y(ii)=[];
            K = length(Mu_X);%K-1;
            Z(Z>ii) = Z(Z>ii) - 1;
        end 
    end
    switch JumpType
        case 1  %Move Mu, Alpha 
            Mu_XTest=Mu_X;
            Mu_YTest=Mu_Y;
            Alpha_XTest=Alpha_X;
            Alpha_YTest=Alpha_Y;
            
            
            for ID=1:K
                if sum(Z==ID) % If not empty 
                    
                    %Get new Mu and Alpha using Gibbs
                    if MaxAlpha>0
                        [Mu_XTest(ID),Alpha_XTest(ID)]=Gibbs_MuAlpha(ID,Z,SMD.X,SMD.FrameNum,SMD.X_SE);
                        [Mu_YTest(ID),Alpha_YTest(ID)]=Gibbs_MuAlpha(ID,Z,SMD.Y,SMD.FrameNum,SMD.Y_SE);
                    else
                        [Mu_XTest(ID)]=Gibbs_Mu(ID,Z,SMD.X,SMD.X_SE);
                        [Mu_YTest(ID)]=Gibbs_Mu(ID,Z,SMD.Y,SMD.Y_SE);
                    end
                    
                    if isnan(Mu_XTest(ID)) || isnan(Mu_YTest(ID)) || isinf(Mu_XTest(ID)) || isinf(Mu_YTest(ID))
                        Mu_XTest(ID) = [];
                        Mu_YTest(ID) = [];
                        Ind = Z==ID;
                        Z(Z>ID) = Z(Z>ID) - 1;
                        if K>1
                            Z(Ind) = randi([1 K-1],[sum(Ind),1]);
                        else
                            Z(Ind) = 1;
                        end
                        Mu_X=Mu_XTest;
                        Mu_Y=Mu_YTest;
                        Alpha_X=Alpha_XTest;
                        Alpha_Y=Alpha_YTest;
                        K = K - 1;
                        continue;
                    end
                    
                    %Check if random number is in bounds, if not then reject
                    if (Mu_XTest(ID)>MaxX)||(Mu_XTest(ID)<MinX)||(Mu_YTest(ID)>MaxY)||(Mu_YTest(ID)<MinY)||...
                            (Alpha_XTest(ID)>MaxAlpha)||(Alpha_XTest(ID)<-MaxAlpha)||(Alpha_YTest(ID)>MaxAlpha)||(Alpha_YTest(ID)<-MaxAlpha)
                    else
                        Mu_X=Mu_XTest;
                        Mu_Y=Mu_YTest;
                        Alpha_X=Alpha_XTest;
                        Alpha_Y=Alpha_YTest;
                    end
                else %No assignments, draw random number from prior
                    Mu_XTest(ID)=Wx*rand()+MinX;
                    Mu_YTest(ID)=Wy*rand()+MinY;
                end
                
            end
            
            %fprintf('K: %d\n',K);
            
            if nn>NBurnin %Then record in chain
                
                K = length(Mu_X);
                Chain(nn-NBurnin).N = length(Mu_X);
                Chain(nn-NBurnin).X = Mu_X';
                Chain(nn-NBurnin).Y = Mu_Y';
                Chain(nn-NBurnin).AlphaX = Alpha_X';
                Chain(nn-NBurnin).AlphaY = Alpha_Y';
                Chain(nn-NBurnin).ID = Z;
                
            end

        case 2  %Reallocation of Z
            
            [ZTest]=Gibbs_Z(SMD, K,Mu_X,Mu_Y,Alpha_X,Alpha_Y);            
            KTest = K; 
            Mu_XTest = Mu_X;
            Mu_YTest = Mu_Y;
            Alpha_XTest = Alpha_X;
            Alpha_YTest = Alpha_Y;
            for ii = K:-1:1
                if sum(Z==ii)==0
                    Mu_XTest(ii)=[];
                    Mu_YTest(ii)=[];
                    Alpha_XTest(ii)=[];
                    Alpha_YTest(ii)=[];
                    KTest = length(Mu_XTest);
                    ZTest(ZTest>ii) = ZTest(ZTest>ii) - 1;
                end 
             end
            
            NZ_Current = hist(Z,1:K);
            NZ_Test = hist(ZTest,1:KTest);
            
            Top = mnpdf(NZ_Test,1/KTest*ones(1,KTest));
            MR=Top/mnpdf(NZ_Current,1/K*ones(1,K));
            A = MR;
            Accept = isnan(A) & Top == 0;
            if rand<A || Accept
                Z=ZTest;
                Mu_X = Mu_XTest;
                Mu_Y = Mu_YTest;
                Alpha_X = Alpha_XTest;
                Alpha_Y = Alpha_YTest;
                K = KTest;
            end

%             if rand<A
%                 Z=ZTest;
%             end
            %fprintf('K: %d\n',K);
            if nn>NBurnin %Then record in chain
                
                Chain(nn-NBurnin).N = length(Mu_X);
                Chain(nn-NBurnin).X = Mu_X';
                Chain(nn-NBurnin).Y = Mu_Y';
                Chain(nn-NBurnin).AlphaX = Alpha_X';
                Chain(nn-NBurnin).AlphaY = Alpha_Y';
                Chain(nn-NBurnin).ID = Z;
                
            end
            
        case 3  %Add
            
            %Add an emitter
            Mu_XTest = cat(2,Mu_X,MinX+Wx*rand);
            Mu_YTest = cat(2,Mu_Y,MinY+Wy*rand);
            KTest = length(Mu_XTest);
            
            if MaxAlpha>0
                Alpha_XTest = cat(2,Alpha_X,-MaxAlpha+2*MaxAlpha*rand);
                Alpha_YTest = cat(2,Alpha_Y,-MaxAlpha+2*MaxAlpha*rand);
            else
                Alpha_XTest = cat(2,Alpha_X,0);
                Alpha_YTest = cat(2,Alpha_Y,0);
            end
            
            %Gibbs allocation
            [ZTest]=Gibbs_Z(SMD,K+1,Mu_XTest,Mu_YTest,Alpha_XTest,Alpha_YTest);
            
            
            for ii = K+1:-1:1
               if sum(ZTest==ii)==0
                   Mu_XTest(ii)=[];
                   Mu_YTest(ii)=[];
                   Alpha_XTest(ii)=[];
                   Alpha_YTest(ii)=[];
                   KTest = length(Mu_XTest);
                   ZTest(ZTest>ii) = ZTest(ZTest>ii) - 1;
               end 
            end
            
            %fprintf('K: %d',K);
            NZ_Current = hist(Z,1:K);
            NZ_Test = hist(ZTest,1:KTest);
            MR = mnpdf(NZ_Test,1/KTest*ones(1,KTest))/mnpdf(NZ_Current,1/K*ones(1,K));
            
            %Likelihood
            LogL_Test = LogPData(SMD,ZTest,Mu_XTest,Mu_YTest,Alpha_XTest,Alpha_YTest);
            LogL_Current = LogPData(SMD,Z,Mu_X,Mu_Y,Alpha_X,Alpha_Y);
            LR = exp(LogL_Test-LogL_Current);
            
            %Prior Raio
            PR = Pk(KTest)/Pk(K);
            
            %Posterior Ratio
            A = PR*LR*MR/KTest;
            Accept = isinf(LogL_Current) & LogL_Current < 0 & ~isinf(LogL_Test);
            
            if rand<A || Accept
                Z=ZTest;
                K=KTest;%K+1;
                Mu_X=Mu_XTest;
                Mu_Y=Mu_YTest;
                Alpha_X=Alpha_XTest;
                Alpha_Y=Alpha_YTest;
            end
            
            if nn>NBurnin %Then record in chain
                
                Chain(nn-NBurnin).N = length(Mu_X);
                Chain(nn-NBurnin).X = Mu_X';
                Chain(nn-NBurnin).Y = Mu_Y';
                Chain(nn-NBurnin).AlphaX = Alpha_X';
                Chain(nn-NBurnin).AlphaY = Alpha_Y';
                Chain(nn-NBurnin).ID = Z;
                
            end
            
        case 4  %Remove
            
            if K==1 %Then update chain and return
                if nn>NBurnin %Then record in chain
                    
                    Chain(nn-NBurnin).N = K;
                    Chain(nn-NBurnin).X = Mu_X';
                    Chain(nn-NBurnin).Y = Mu_Y';
                    Chain(nn-NBurnin).AlphaX = Alpha_X';
                    Chain(nn-NBurnin).AlphaY = Alpha_Y';
                    Chain(nn-NBurnin).ID = Z;
                    
                end
                continue;
            end
            
            %pick emitter to remove:
            ID =randi(K);
            
            Mu_XTest = Mu_X;
            Mu_YTest = Mu_Y;
            Alpha_XTest = Alpha_X;
            Alpha_YTest = Alpha_Y;
            
            %Remove from list
            Mu_XTest(ID) = [];
            Mu_YTest(ID) = [];
            Alpha_XTest(ID) = [];
            Alpha_YTest(ID) = [];
            KTest = length(Mu_XTest);
            
            %Gibbs allocation
            [ZTest]=Gibbs_Z(SMD,K-1,Mu_XTest,Mu_YTest,Alpha_XTest,Alpha_YTest);
            
            %fprintf('K: %d\n',K);
            for ii = K-1:-1:1
               if sum(ZTest==ii)==0
                   Mu_XTest(ii)=[];
                   Mu_YTest(ii)=[];
                   Alpha_XTest(ii)=[];
                   Alpha_YTest(ii)=[];
                   KTest = K-1;
                   ZTest(ZTest>ii) = ZTest(ZTest>ii) - 1;
               end 
            end
            
            NZ_Current = hist(Z,1:K);
            NZ_Test = hist(ZTest,1:KTest);
            MR = mnpdf(NZ_Test,1/KTest*ones(1,KTest))/mnpdf(NZ_Current,1/K*ones(1,K));
            
            %Likelihood
            LogL_Test = LogPData(SMD,ZTest,Mu_XTest,Mu_YTest,Alpha_XTest,Alpha_YTest);
            LogL_Current = LogPData(SMD,Z,Mu_X,Mu_Y,Alpha_X,Alpha_Y);
            LR = exp(LogL_Test-LogL_Current);
            
            %Prior Raio
            PR = Pk(KTest)/Pk(K);
            
            %Posterior Ratio
            A = PR*LR*MR*K;
            
            if rand<A
                Z=ZTest;
                K=KTest;%K-1;
                Mu_X=Mu_XTest;
                Mu_Y=Mu_YTest;
                Alpha_X=Alpha_XTest;
                Alpha_Y=Alpha_YTest;
            end
            
            if nn>NBurnin %Then record in chain
                
                Chain(nn-NBurnin).N = length(Mu_X);
                Chain(nn-NBurnin).X = Mu_X';
                Chain(nn-NBurnin).Y = Mu_Y';
                Chain(nn-NBurnin).AlphaX = Alpha_X';
                Chain(nn-NBurnin).AlphaY = Alpha_Y';
                Chain(nn-NBurnin).ID = Z;
                
            end
            
    end

    if DEBUG %for testing
        figure(1111)
        %plot(SMD.X,SMD.Y,'ko','linewidth',2)
        scatter(SMD.X,SMD.Y,[],Z)
        hold on
        plot(Mu_X,Mu_Y,'ro','linewidth',4)
        %text(60,50,sprintf('Jump: %g',nn),'FontSize',12)
        legend(sprintf('Jump: %g',nn))
%         xlim([0 80]);
%         ylim([0 55])
        xlabel('X(nm)')
        ylabel('Y(nm)')
        %set(gca,'FontSize',12)
        hold off
        pause(.001)
        %fprintf('Jump: %d\n',JumpType);
%        fprintf('A: %g\n', A);
%         Frame = getframe(gcf);
%         writeVideo(V,Frame);
    end
end
% close(V)
Results.ZChain=ZChain;
Results.KChain=KChain;
Results.Xk=Xk;
Results.Yk=Yk;
Results.AXk=AXk;
Results.AYk=AYk;


end

function LogL=LogPData(SMD,Z,Mu_X,Mu_Y,Alpha_X,Alpha_Y)
    %This is the likelihood P(Data|Model)
    
    X=SMD.X;
    Y=SMD.Y;
    T=SMD.FrameNum;
    SigmaX=SMD.X_SE;
    SigmaY=SMD.Y_SE;
    
    if size(Z,2) == 1
        Z = Z'; 
    end
    
    LogL=sum(log(normpdf(X,Mu_X(Z)'+Alpha_X(Z)'.*T,SigmaX)));
    LogL=LogL+sum(log(normpdf(Y,Mu_Y(Z)'+Alpha_Y(Z)'.*T,SigmaY)));
    
end

function [ZTest,LnP]=Gibbs_Z(SMD,K,Mu_X,Mu_Y,Alpha_X,Alpha_Y)
    %This function calculates updated allocations (Z)
    
    T=SMD.FrameNum;
    N=length(T);
    PX=zeros(N,K);
    PY=zeros(N,K);
    for kk=1:K
        PX(:,kk)=normpdf(SMD.X-(Mu_X(kk)+Alpha_X(kk)*T),0,SMD.X_SE);
        PY(:,kk)=normpdf(SMD.Y-(Mu_Y(kk)+Alpha_Y(kk)*T),0,SMD.Y_SE);
    end
    
    P=PX.*PY+eps;
    PNorm=P./repmat(sum(P,2),[1 K]);

    if sum(sum(isnan(P)))
       [ZTest] = knnsearch([Mu_X',Mu_Y'],[SMD.X,SMD.Y]); 
    else 
        ZTest=K+1-sum(repmat(rand(N,1),[1,K])<cumsum(PNorm,2),2);
    end
    LI=sub2ind(size(P),(1:size(P,1))',ZTest);
    LnP=sum(log(PNorm(LI))); 
          
end

function [Mu,Alpha]=Gibbs_MuAlpha(ID,Z,X,T,Sigma)
    %This function calculates updated Mu and Alpha (1D)
    
    %Get the localizations from the IDth emitter
    Xs=X(Z==ID);
    Sigs = Sigma(Z==ID);
    Ts=T(Z==ID);
    
    A = sum(Sigs.^-2);
    B = sum(Ts./Sigs.^2);
    D = sum((Ts.^2)./(Sigs.^2));
    
    %MLE estimates of Mu and Alpha
    
    [Alpha,Center] = calAlpha(Xs,Sigs,Ts);
    MA=[Center;Alpha];

    %Covariance matrix Sigma
    COV = pinv([A, B;B,D]);
    
    %This draws [Mu,Alpha] from a multivariate normal
    MuAlpha=mvnrnd(MA,COV);
    Mu=MuAlpha(1);
    Alpha=MuAlpha(2);
    
    if Mu == Center
       Mu = Center + sqrt(A)*randn(); 
    end
    
end

function [Mu]=Gibbs_Mu(ID,Z,X,Sigma)
    %This function calculates updated Mu (1D)
    
    %Get the localizations from the IDth emitter
    Xs=X(Z==ID);
    Sigs = Sigma(Z==ID);
  
    A = sum(Xs./(Sigs.^2));
    B = sum(Sigs.^-2);
    XMLE = A/B;
    X_SE = 1/sqrt(B);
    
    Mu=normrnd(XMLE,X_SE);
end

function [Alpha,Center] = calAlpha(Xs,Sigs,Frames)
Frames = single(Frames);
A = sum(Xs./Sigs.^2);
B = sum(Frames./Sigs.^2);
C = sum(Sigs.^-2);
AlphaTop = sum((C*Xs-A).*Frames./Sigs.^2);
AlphaBottom = sum((C*Frames-B).*Frames./Sigs.^2);
Alpha = AlphaTop/AlphaBottom;
Center = (A-Alpha*B)/C;
end

