%% Preamble
% Program: MainFileNew
% Author: Jonathan Myers (Modified From S. Tay et al. 2010 Nature)
% Date: July 21st, 2019
% Purpose: To make the Tay MATLAB MainFile script more streamlined and
% easier to understand and modify.
% Arguments: None.
% Calls: ParametersNew, StatusChangeNew (ModelNew), and AllCellPlottingNew.
% Returns: None. Generates 'last.mat' file of simulation history and
% variables.

%% Setup
clear;              %reset all
clc;                %clear comand window
starttime=clock;    %current time

% Not Needed: Vestigial code from likely other code. 
%rand('twister', sum(1000*clock));

TNF=10; % TNF Dose

ANa=2; % ANa=2 - # IKBa alleles

AN=2; % AN=2 - # A20 alleles

ANR=2; % ANR=2 - # Reporter gene alleles

% Set AN=0 to study A20 knockout

N=10; % Number of cell to be simulated

%% Simulation Time Points
% Various time protocols can be studied within this frame.
% 10h randomization of initial conditions
t000=10*3600;
% 10h equilibrium waiting time
t00=10*3600;
% 1 step, time when TNF is being introduced into the system (in seconds)
t0=50*60;
% 2 step, length of TNF stimulation
tw1=5*60;
% 3 step length of first break   (White breaks: 3600s, 6000s, 12000s, our break 170*60s)
te1=100*60;
% 4 step length of second TNF stimulation
tw2=5*60;
% 5 step length of second break
te2=100*60;
% 6 step length of third TNF stimulation
tw3=5*60;
% 7 step length of third break
te3=100*60;

%% Numerical Simulation Parameters
tt=1000;  % time interval for ODEs solving 
                     
YYY=0;  % matrix of average, all variables y0(i)(t) 
NFKB=0; % total nuclear NF-kB 
GGa=0;  % status of Ikba gene
GG=0;   % status of A20 gene
GGT=0;  % status of TNF 
GGR=0;  % status of reporter genes
Bb=0;   % number of active receptors
MM=0;
NFF=0;

%% Start
for i=1:N % Simulate each cell (i)
    
i % Cell Number
    
%Call Parameters
[NF0,NF1,NF2,M0,M1,M2,k4,ka20,AB,kv,q1,q2,c1,c3,c4,c5,k1,k2,k3,a1,a2,a3,...
c1a,c5a,c6a,i1,i1a,e1a,e2a,dt,tp,KN,KNN,ka,ki,kb,kf,Tdeg,q1r,q2r,q2rr,...
c1r,c1rr,c3r]=Parameters;
    
%% Randomize Initial Total TNF Receptors And NF-kB Conditions
NF=round(NF0*exp(NF2+randn*NF1))                %NF-kB level
while NF > 10*NF0
NF=round(NF0*exp(NF2+randn*NF1))
end
     
%NF=NF0;     %uncomment to remove extrinsic noise 
   
% Lognormal distribution with Median=NF0, Mean=NF0*Exp(NF1^2/2), %
% Variance = NF0^2 * (Exp (NF1^2 -1) * Exp(NF1^2)
   

M=round(M0*exp(M2+randn*M1))            % number of TNFR1 receptors
while M > 10*M0
M=round(M0*exp(M2+randn*M1))
end
     
%M=M0;         %uncomment to remove extrinsic noise  
   
% Lognormal distribution with Median=M0, Mean=M0*Exp(M1^2/2), %
% Variance = M0^2 * (Exp (M1^2 -1) * Exp(M1^2)

%% Accumulator Setup
y0=zeros(1,19);     %initial conditions set to zero and next: 

y0(14)=NF;          %NF-kB is given in cytoplasmic complex(IkBa|NFkB), standard = 10^5
y0(2)=2*10^5;       %initial IKKn, total IKK kept constant  
y0(11)=0.14*y0(14); %free cytoplasmic IkBa protein 
y0(12)=0.06*y0(14); %free nuclear IkBa protein
y0(13)=10;          %IkBa mRNA
y0(10)=10;          %10 A20 mRNA
y0(9)=10000;        %10000 A20 protein
   
y0(10)=AB*y0(10);                   
y0(9)=AB*y0(9); 

Ga=0;               % initial status of IkBa promoter
G=0;                % initial status of A20 promoter
GT=0;               % initial status of TNF promoter
GR=0;               % initial status of reporter gene promoter
B=0;                % initial number of active receptors
yy0=y0;             % initial conditions y0(i)   

%% Randomization of Remaining Initial Conditions
realtime=0;                     %simulated time   
phase=round(rand*t000/dt)*dt;   %random initial time (dt -simulation time step -10s)  
tspan=[0:dt:tt];                %time for which the solution is derived to find the switching time, tt=1h

while (realtime<phase)
    [T0,Y0]=ode23tb(@ModelNew,tspan,yy0,[],Ga,G,GR,B,M); 
    Yact=Y0(:,8);               %amount of NF-kBn  
    Yin=Y0(:,12);               %amount of IkBan 
    TR=Y0(:,16);                %TNF level  
    Gax=Ga;Gx=G;GRx=GR;Bx=B;
    [mk,Ga,G,GR,B]=StatusChangeNew(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M) ;  %function determining the change of gene status, calls statuschange

    tc=T0(mk);                  %time when the status changes 
    yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration  
    realtime=realtime+tc;
end;

if (realtime>phase)
    nn=(realtime-phase)/dt;
    yy0=Y0(mk-nn,:);
    Ga=Gax;G=Gx;GR=GRx;B=Bx;
end;                            %status before the last change it occured outside of the time interval

clear Yact Yin Y0 T0 nn mk phase tc;

%####################################################
%###### 0 step - waiting for "equilibrium"    #######
%####################################################

realtime=0;

while (realtime<t00)
    [T0,Y0]=ode23tb(@ModelNew,tspan,yy0,[],Ga,G,GR,B,M); 
    Yact=Y0(:,8);               %amount of NF-kBn  
    Yin=Y0(:,12);               %amount of IkBan  
    TR=Y0(:,16);                %TNF level  
    Gax=Ga;Gx=G;GRx=GR;Bx=B;
    [mk,Ga,G,GR,B]=StatusChangeNew(AN,ANa,ANR,TR,Gax,Gx,GRx,Bx,Yact,Yin,M); %function determining the change of gene status, calls statuschange

    tc=T0(mk);                  %time when the status changes 
    yy0=Y0(mk,:);               %transfer of initial conditions to the next iteration  
    realtime=realtime+tc;
end;

if (realtime>t00)
    nn=(realtime-t00)/dt;
    yy0=Y0(mk-nn,:);
    Ga=Gax;G=Gx;GR=GRx;B=Bx;
end;                            %status before the last change it occured ouside of the time interval

clear Yact Yin Y0 T0 nn mk tc; 
