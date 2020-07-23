% SEIRD model

clear all                   % clear the workspace (memory)
close all                   % close all previous plots

% Input parameters
load FB2404                % loads the adjacency matrix A

N=size(A,1);
beta=0.08;                  % assume 8% spread rate
alpha=1/5;                  % rate E->I assuming 5 days incubation
delta=1/14;                 % rate I->R assuming 14 days recovery
omega=0.04;                 % rate I->D assuming 4% mortality rate
gamma=1/5;                  % rate I->H: people that need hospitalization

% Initialize the state of the population
% Initially exposed
Seeds=10;                   % number of exposed seeds
E=zeros(size(A,1),1);       % initial infected seed
for seed=1:Seeds
    RandSeed=ceil(N*rand);
    E(RandSeed)=1;
end

S=ones(N,1)-E;     % initial healthy people
I=zeros(N,1); D=zeros(N,1); R=zeros(N,1); H=zeros(N,1);
TotalI(1)=sum(I); TotalE(1)=sum(E); TotalR(1)=sum(R); TotalH(1)=sum(H);

mailedTestingWorks = 0;

% graph for lockdown
Arand=A;
for nodei=1:N
    for nodej=nodei+1:N
        edgecoin=rand<0.5;      % remove a percentage of edges
        if edgecoin==1
            Arand(nodei,nodej)=0; Arand(nodej,nodei)=0;
        end
    end
end
t=1;


% Run iterations
while (sum(I)+sum(E)>0)
    
    % Compute total numbers per state
    TotalS(t)=sum(S);
    TotalE(t)=sum(E);
    TotalI(t)=sum(I);
    TotalR(t)=sum(R);
    TotalD(t)=sum(D);
    TotalH(t)=sum(H);
    
    % Always in Lockdown
    lock(t) = 1;
    
    % Transition from S to E
    NI=Arand*I;  % vector with number of infected neighbors
    
    mask(t)=0;
    maskscount=1000;
    peopleWithMaskDis=0;
         
    if maskscount>0 %if maskscount = 0, then the dispenser will no longer run
        maskDis(t)=rand()<0.65;  %random number of people who get mask
        if(maskDis(t)==1)
            mask(t)=1;
            peopleWithMaskDis= peopleWithMaskDis+1; %the people who get a mask from the dispenser increases by 1
            maskscount=maskscount - 1;   %the number of masks availiable decreases when someone takes a mask
        end
    end
       
    if mask(t)==1
        NewE=rand(N,1)<1-(1-(beta-0.14)).^(NI); % mask used
        NewE=and(NewE,maskDis(t));
    else
        NewE=rand(N,1)<1-(1-beta).^(NI);       % no masks used
    end
    
    % Testing mailing implementation
    for nodei=1:N
        if NI(nodei) > 0
            edgecoin=rand<0.98;
            if edgecoin==1
                mailedTestingWorks = mailedTestingWorks + 1;
                for nodej=nodei+1:N
                    Arand(nodei, nodej) = 0; Arand(nodej, nodei) = 0;
                end
            end
        end
    end
    
    % Transition from E to I
    CoinE2I=rand(N,1)<alpha;
    NewI=and(CoinE2I,boolean(E));
    if rand(size(A,1),1) < 0.001             %max rate virus spreads when touching surfaces
        touchInfected(t) = 1;
        NewI = and(NewI, touchInfected(t));
    end
    
    % Transition from I to H
    RandomNumber=rand(N,1);
    CoinI2H=RandomNumber>1-gamma;
    NewH=and(CoinI2H,boolean(I));
    
    if ~CoinI2H
    % Transition from I to R or D
        RandomNumber=rand(N,1);
        CoinI2R=RandomNumber<delta;
        CoinI2D=RandomNumber>1-omega;
        NewR1=and(CoinI2R,boolean(I));
        NewD1=and(CoinI2D,boolean(I));
    else 
        NewR1=and(0,boolean(I));
        NewD1=and(0,boolean(I));
    end     
          
    % Transition from H to R or D
    RandomNumber=rand(N,1);
    CoinH2R=RandomNumber<delta;
    CoinH2D=RandomNumber>1-omega;
    NewR2=and(CoinH2R,boolean(H));
    NewD2=and(CoinH2D,boolean(H));
    
    % Update indicator vectors
    S=S-NewE;
    E=E+NewE-NewI;
    I=I+NewI-NewH-NewR1-NewD1;
    H=H+NewH-NewR2-NewD2;
    R=R+NewR1+NewR2;
    D=D+NewD1+NewD2;
    
    t=t+1

end

TotalDeaths=TotalD(t-1)         % total number of deaths
DaysLocked=sum(lock)            % days in locked down
DaysCOVID=t                     % days before eradication

% Plot results
figure;
plot(TotalS/N,'b'); hold on;
plot(TotalI/N,'g');
figure;
plot(TotalE,'color',[0.9100    0.4100    0.1700]); hold on
plot(TotalI,'r'); 
plot(TotalH,'k'); 
plot(TotalD,'m');
plot(TotalR,'g');
