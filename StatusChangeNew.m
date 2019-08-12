%% Preamble
% Program: StatusChangeNew
% Author: Jonathan Myers (Modified From S. Tay et al. 2010 Nature)
% Date: July 9th, 2019
% Purpose: To make the Tay MATLAB StatusChange script more streamlined and
% easier to understand and modify. Uses Gillespie Algorithm.
% Arguments: Gene Copy Number, Extracellular TNF-alpha Level, Gene
% Activaition States, Active TNFR1 Receptor Number, Nuclear NF-kB Level, Nuclear
% IkBa Level, and Total TNFR1 Receptor Number.
% Calls: Parameters.
% Returns: Reaction Time Index, Gene Activation States, and Active TNFR1
% Receptor Number.

%% StatusChangeNew
function [mk,Ga,G,GR,B]=StatusChangeNew(AN,ANa,ANR,TRx,Gax,Gx,GRx,Bx,Yact,Yin,M)
%Initalize Gene States and Active Receptor Number
Ga=Gax;
G=Gx;
GR=GRx;
B=Bx;

%Call Parameters
[NF0,NF1,NF2,M0,M1,M2,k4,ka20,AB,kv,q1,q2,c1,c3,c4,c5,k1,k2,k3,a1,a2,a3,...
    c1a,c5a,c6a,i1,i1a,e1a,e2a,dt,tp,KN,KNN,ka,ki,kb,kf,Tdeg,q1r,q2r,q2rr,...
    c1r,c1rr,c3r]=Parameters;
%% Pick Time Interval
%Individual Reaction Propensities
p1a=(ANa-Gax)*q1*Yact;  % risk of NF-kB association to IKBa site at time mk 
p2a=Gax*(q2*Yin);       % risk of  NF-kB dissociation from IkBa at time mk

p1=(AN-Gx)*q1*Yact;     % risk of NF-kB association to A20 site 
p2=Gx*(q2*Yin);         % risk of  NF-kB dissociation from A20 site

p1r=(ANR-GRx)*q1r*Yact; % risk of NF-kB association to reporter gene 
p2r=GRx*(q2r*Yin+q2rr); % risk of  NF-kB dissociation from reporter gene

p3=(M-Bx)*kb*TRx;       % risk of TNFR1-TNF binding   
p4=Bx*kf;               % risk of TNFR1 inactivation

%Find "Total Reaction Propensity"
ro=(p1a+p2a+p1+p2+p1r+p2r+p3+p4);

%Using Cumulative Distribution Function of Exponential Probability
%Distribution, Find "Total Reaction Propensity Probability"
roint=dt*cumtrapz(ro);        % propensity function integrated using trapezoidal method
fd=1-exp(-roint);             % Distribution of the switching time

%Generate P1 to Determine Switching Time in Stochastic Simulation
%If P1 Larger than "Total Reaction Propensity Probability," Then Do Not
%Change Stochastic Reaction States
%If P1 Smaller than "Total Reaction Propensity Probability," Then Change
%Stochastic Reaction States

r=rand;
mk=length(fd); %mk = index of last element of fd array

if (fd(end)>=r)  
a=abs(fd-r);
mk=find((a-min(a))==0);         % mk = index (time) of next reaction 

clear a fd ro roint;      %Clear No Longer Needed Variables    

%% Pick Reaction Occurance
%Determine Reaction State to Change Based on (0,1) Range Normalization of
%Reaction Propensities Using P2

rnumber=rand;

ss= p1a(mk) + p2a(mk) +p1(mk) + p2(mk) + p1r(mk) + p2r(mk) + p3(mk) + p4;
p1a=p1a(mk)/ss;
p2a=p2a(mk)/ss;
p1=p1(mk)/ss;
p2=p2(mk)/ss;
p1r=p1r(mk)/ss;
p2r=p2r(mk)/ss;
p3=p3(mk)/ss;
p4=p4/ss;

if (rnumber<p1a) 
Ga=Ga+1;end   %IKBa activates
if (rnumber>=p1a)&&(rnumber<p1a+p2a)                 
Ga=Ga-1;end  %IKBa inactivates

if (rnumber>=p1a+p2a)&&(rnumber<p1a+p2a+p1)                    
G=G+1;end     %A20 activates
if (rnumber>=p1a+p2a+p1)&&(rnumber<p1a+p2a+p1+p2)             
G=G-1;end    %A20 inactivates

if (rnumber>=p1a+p2a+p1+p2)&&(rnumber<p1a+p2a+p1+p2+p1r)             
GR=GR+1; end    %reporter gene activates
if (rnumber>=p1a+p2a+p1+p2+p1r)&&(rnumber<p1a+p2a+p1+p2+p1r+p2r)      
GR=GR-1; end    %reporter gene inactivates

if (rnumber>=p1a+p2a+p1+p2+p1r+p2r)&&(rnumber<p1a+p2a+p1+p2+p1r+p2r+p3)          
B=B+1;end   %receptor activation
if (rnumber>=p1a+p2a+p1+p2+p1r+p2r+p3)&&(rnumber<p1a+p2a+p1+p2+p1r+p2r+p3+p4)   
B=B-1;end      %receptor deactivation 

end
end
