% Implementation of Mask dispensing machines on streets, hotels, schools… (Audrey)

load FB2404copie  

A=zeros(11);
   
I=zeros(size(A,1),1); 
E=zeros(size(A,1),1);
S=ones(size(A,1),1)-E;
TotalI(1)=sum(I); TotalE(1)=sum(E);
TotalS(1)=sum(S);


N=2404;
t=1;


%  Mask using strategy (reduction beta)
%     if TotalI(t)/N>0.005
%         mask(t)=1;
%         betaeff(t)=0.4*beta;        % mask used
%     elseif TotalI(t)/N<0.0001;
%         mask(t)=0;
%         betaeff(t)=beta;            % no masks used
%     else
%         mask(t)=mask(t-1);
%         betaeff(t)=betaeff(t-1);
%     end
    
    

%   random number of people who get the mask

     
         mask(t)=0;
         maskscount=1000;
         peopleWithMaskDis=0;
         
         while(maskscount>0) %if maskscount = 0, then the dispenser will no longer run
             maskDis(t)=rand()<0.65;  %random number of people who get mask

             if(maskDis(t)==1)
                 mask(t)=1;
                 peopleWithMaskDis= peopleWithMaskDis+1; %the people who get a mask from the dispenser increases by 1
                 maskscount=maskscount - 1;   %the number of masks availiable decreases when someone takes a mask
                 touchDis=rand(size(A,1),1); %random number of people who touch dispenser
                    if touchDis< 0.001 %max rate virus spreads when touching surfaces      
                        NewI= NewI + touchDis;
                    end
             end
             t=t+1;
         end
       
   
     if mask(t)==1
                NewE=rand(N,1)<1-(1-betaeff(t)).^(NI); % mask used
                NewE=NewE + maskDis(t);
           else
                NewE=rand(N,1)<1-(1-beta).^(NI); % no masks used
     end
   
% number of people who get mask will be added to the exposed      
% how many people get the masks code with random number (done)
% people who get infected by touching the machine (done) 
% people who wouldn't actually take the mask #
% number of masks left in dispenser (done)