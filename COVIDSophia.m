% Implementation for Day to day in-person house checks for temperature: if suspected to be infected by COVID19, 
% move to hospital immediately (Sophia)
% SEIRD hospital model

clear all                   % clear the workspace (memory)
close all                   % close all previous plots

% Input parameters
load FB2404                 % loads the adjacency matrix A
N=size(A,1);
%A=ones(2404)-eye(2404);
beta=20/2404;               % spreading rate
alpha=1/5;                  % rate E->I
gamma=0.9;                  % rate I->H
delta=1/10;                 % rate H->R
omega=(0.10/0.95)*delta;    % rate H->D (assume 10% mortality rate)

% Initialize the state of the population
Seeds=floor(0.005*N);                    % number of exposed seeds
E=zeros(N,1);       % initial infected seed
for seed=1:Seeds
    RandSeed=ceil(N*rand);
    E(RandSeed)=1;
end
S=ones(N,1)-E;     % initial healthy people
I=zeros(N,1); R=zeros(N,1); D=zeros(N,1); H=zeros(N,1);

% Run iterations
t=0;
while (sum(E)+sum(I)>0)

    % Update variables
    t=t+1;
    TotalS(t)=sum(S);
    TotalE(t)=sum(E);
    TotalI(t)=sum(I);
    TotalH(t)=sum(H);
    TotalR(t)=sum(R);
    TotalD(t)=sum(D);

    % Transition from S to E
    NI=A*I;  % vector with number of infected neighbors
    NewE=rand(N,1)<1-(1-beta).^(NI);
    NewE=and(NewE,boolean(S));   % boolean vector indicating susceptible people that transfer into Exposed
        
    % Transition from E to I
    CoinE2I=rand(N,1)<alpha;
    NewI=and(CoinE2I,boolean(E));
    
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
end

% Count final results
TotalDeaths=TotalD(t-1)         % total number of deaths
DaysCOVID=t                     % days before eradication

% Plot results
figure;
plot(TotalS,'b'); hold on
plot(TotalR,'g')
figure;
plot(TotalE,'color',[0.9100    0.4100    0.1700]); hold on
plot(TotalI,'r')
plot(TotalH,'k')
plot(TotalD,'m')

save hospital.mat