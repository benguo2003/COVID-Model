% SEIRD model

clear all                   % clear the workspace (memory)
close all                   % close all previous plots

% Input parameters
load FB2404                % loads the adjacency matrix A

N=size(A,1);
beta=0.056;                 % assume 5.6% spread rate
alpha=1/5;                  % rate E->I assuming 5 days incubation
delta=1/14;                 % rate I->R assuming 14 days recovery
omega=0.04;                 % rate I->D assuming 4% mortality rate
gamma=1/5;                  % rate I->H: people that need hospitalization
epsilon=0.35;               % rate E->A: people that are asymptomatic

% Initialize the state of the population
% Initially exposed
Seeds=100;                   % number of exposed seeds
E=zeros(size(A,1),1);       % initial infected seed
for seed=1:Seeds
    RandSeed=ceil(N*rand);
    E(RandSeed)=1;
end

S=ones(N,1)-E;

I=zeros(N,1); D=zeros(N,1); R=zeros(N,1); H=zeros(N,1); Asym=zeros(N,1);
TotalE(1)=sum(E); TotalR(1)=sum(R); TotalH(1)=sum(H); TotalA(1)=sum(Asym);

withMasks = zeros(N,1);
maskscount=1000;
while maskscount>0    %if maskscount = 0, then the dispenser will no longer run
    x = randi(N);
    if withMasks(x) ~= 1
        withMasks(x) = 1;
        maskscount = maskscount - 1;
    end
end

temp = rand(N,1) < 0.001;
I = and(temp, withMasks);
TotalI(1)=sum(I);

S=S-I;     % initial healthy people

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
while (sum(I)+sum(E)+sum(Asym)+sum(H)>0)
    
    % Compute total numbers per state
    TotalS(t)=sum(S);
    TotalE(t)=sum(E);
    TotalI(t)=sum(I);
    TotalR(t)=sum(R);
    TotalD(t)=sum(D);
    TotalH(t)=sum(H);
    TotalA(t)=sum(Asym);
    
    % Always in Lockdown
    lock(t) = 1;
    
    % Transition from S to E
    
    NI=Arand*or(I,Asym);             % vector with number of infected neighbors
    
    for i=1:N
        if withMasks(i) == 1         % yes masks
            NewE(i) = rand() < 1-(1-(beta-(0.8*beta))).^(NI(i));
        else                         % no masks
            NewE(i) = rand() < 1-(1-beta).^(NI(i));
        end
    end
    
    if t == 1
        NewE = NewE';
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
    CoinE2I=rand(N,1)<(1-epsilon)*alpha;
    NewI=and(CoinE2I,boolean(E));
        
    % Transition from E to Asym
    CoinE2Asym=rand(N,1)<epsilon*alpha;
    NewA=and(CoinE2Asym,boolean(E));
    NewA=and(NewA,not(NewI));
    
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
    
    % Transition from Asym to R
    RandomNumber=rand(N,1);
    CoinAsym2R=RandomNumber<delta;
    NewR3=and(CoinAsym2R,boolean(Asym));
    
    % Update indicator vectors
    S=S-NewE;
    E=E+NewE-NewI-NewA;
    I=I+NewI-NewH-NewR1-NewD1;
    H=H+NewH-NewR2-NewD2;
    R=R+NewR1+NewR2;
    D=D+NewD1+NewD2;
    Asym=Asym+NewA-NewR3;
    
    t=t+1

end

TotalDeaths=TotalD(t-1)         % total number of deaths
DaysLocked=sum(lock)            % days in locked down
DaysCOVID=t                     % days before eradication

% Plot results
figure;
plot(TotalS/N,'b'); hold on;
plot(TotalI/N,'g');
plot(TotalE/N,'color',[0.9100    0.4100    0.1700])
plot(TotalI/N,'r'); 
plot(TotalH/N,'k'); 
plot(TotalD/N,'m');
plot(TotalR/N,'g');
plot(TotalA/N, 'c');