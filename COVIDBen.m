% Implementation for Mail all residents (residents with connection to infected) 
% COVID19 testing equipment for daily self checks & reports (Ben)

% SEIRD model

clear all                   % clear the workspace (memory)
close all                   % close all previous plots

% Input parameters
load FB2404                % loads the adjacency matrix A

N=size(A,1);
beta=0.08;                  % assume 8% spread rate
alpha=1/5;                  % rate E->I assuming 5 days incubation
delta=1/14;                 % rate I->R assuming 14 days recovery
omega=0.04;    % rate I->D assuming 4% mortality rate

% Initialize the state of the population
% Initially exposed
Seeds=10;                   % number of exposed seeds
E=zeros(size(A,1),1);       % initial infected seed
for seed=1:Seeds
    RandSeed=ceil(N*rand);
    E(RandSeed)=1;
end

S=ones(N,1)-E;     % initial healthy people
I=zeros(size(A,1),1); D=zeros(size(A,1),1); R=zeros(size(A,1),1);
TotalI(1)=sum(I); TotalE(1)=sum(E); TotalR(1)=sum(R);

% People with masks
Masks=floor(0.5*N);
M=zeros(N,1);       % initial infected seed
for mask=1:Masks
    RandMask=ceil(N*rand);
    M(RandMask)=1;
end

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

%%

% Run iterations
while (sum(I)+sum(E)>0)
    
    % Compute total numbers per state
    TotalS(t)=sum(S);
    TotalE(t)=sum(E);
    TotalI(t)=sum(I);
    TotalR(t)=sum(R);
    TotalD(t)=sum(D);
    
    % Always in Lockdown
    lock(t) = 1;
    
    %%  Transition
    
    % Transition from S to E
    NI=Arand*I;  % vector with number of infected neighbors
    NIM=Arand*(I.*M);  % vector with number of infected neigh wearing masks
    NIN=NI-NIM; % vector with number of infected neigh not wearing masks
    
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
    
    % Mask using strategy where a fraction of people wear masks
    if TotalI(t)/N>0.005
        mask(t)=1;
    elseif TotalI(t)/N<0.0001;
        mask(t)=0;
    else
        mask(t)=mask(t-1);
    end
    
    if mask(t)==1
        NewE=rand(N,1)<1-(1-(1-0.3*M)*beta).^NIN.*(1-(0.05-0.035*M)*beta).^NIM;  % mask used
    else
        NewE=rand(N,1)<1-(1-beta).^(NI);                              % no masks used
    end
    
    
    
    NewE=and(NewE,boolean(S));   % boolean vector indicating susceptible people that transfer into Exposed
        
    % Transition from E to I
    CoinE2I=rand(N,1)<alpha;
    NewI=and(CoinE2I,boolean(E));
        
    % Transition from I to R or D
    RandomNumber=rand(N,1);
    CoinI2R=RandomNumber<delta;
    CoinI2D=RandomNumber>1-omega;
    NewR=and(CoinI2R,boolean(I));
    NewD=and(CoinI2D,boolean(I));    
    
    % Update indicator vectors
    S=S-NewE;
    E=E+NewE-NewI;
    I=I+NewI-NewR-NewD;
    R=R+NewR;
    D=D+NewD;
    %%
    t=t+1
end
TotalDeaths=TotalD(t-1)         % total number of deaths
DaysLocked=sum(lock)            % days in locked down
DaysCOVID=t                     % days before eradication

% Plot results
figure;
plot(TotalS/N,'b'); hold on
plot(TotalE/N,'color',[0.9100    0.4100    0.1700])
plot(TotalI/N,'r')
plot(TotalR/N,'g')
plot(TotalD/N,'m')
figure; plot(TotalI/N,'r')
figure; plot(mask,'b'); hold on; plot(lock,'k'); axis([1 t -0.1 1.1])
