%% Preamble
% Program: AverageCellPlottingNew
% Author: Jonathan Myers (Modified From S. Tay et al. 2010 Nature)
% Date: July 21st, 2019
% Purpose: To make the Tay MATLAB AverageCellPlotting script more streamlined and
% easier to understand and modify.
% Arguments: None. Uses variables T, YYY, GGa, GG, GGR, Bb that should be
% in Workspace memory.
% Calls: None.
% Returns: None. Generates Plots of Molecular Concentrations and Gene States.

%% Figure 1
figure(1)
set(gcf,'Color',[1,1,1])

subplot(4,2,1);
plot(T,YYY(:,16));
grid on;
title('  TNF activity');

subplot(4,2,2);
plot(T,YYY(:,1));
grid on;
title('  Active IKKK, (IKKKp)');

subplot(4,2,3);
plot(T,KNN-YYY(:,2)-YYY(:,3)-YYY(:,4));
grid on;
title('  Intermediate Inactive IKK, (IKKii)');

subplot(4,2,4);
plot(T,YYY(:,2));
grid on;
title('  Neutral IKK, (IKKn)');

subplot(4,2,5);
plot(T,YYY(:,3));
grid on;
title('  Active IKK, (IKKa)');

subplot(4,2,6);
plot(T,YYY(:,4));
grid on;
title('  Inactive IKK, (IKKi)');

subplot(4,2,7);
plot(T,YYY(:,5));
grid on;
title('  phosphorylated IkBa');
xlabel('Time in min')

subplot(4,2,8);
plot(T,YYY(:,6));
grid on;
title('  phosphorylated IkBa|NFkB');
xlabel('Time in min')

%% Figure 2
figure(2)
set(gcf,'Color',[1,1,1])

subplot(4,2,1);
plot(T,YYY(:,3));
grid on;
title(' Active IKK, (IKKa)');

subplot(4,2,2);
plot(T,YYY(:,12)+YYY(:,15)+YYY(:,11)+YYY(:,14)+YYY(:,5)+YYY(:,6));
grid on;
title('  Total IkBa');

subplot(4,2,3)
plot(T,YYY(:,8));
grid on;
title('  Free nuclear NF-kB');

subplot(4,2,4);
plot(T,GGa);
grid on;
title('  Activity of IkBa gene, both alleles ');

subplot(4,2,5);
plot(T,GG);
grid on;
title('   Activity of A20 gene,both alleles'); 

subplot(4,2,6);
plot(T,YYY(:,13));
grid on;
title('  IkBa mRNA transcript');

subplot(4,2,7);
plot(T,YYY(:,10));
grid on;
title('  A20 mRNA transcript');
xlabel('Time in min');
subplot(4,2,8);

plot(T,YYY(:,9));
grid on;
title(' A20 protein');
xlabel('Time in min')

%% Figure 3
figure(3)
set(gcf,'Color',[1,1,1])

subplot(3,2,1);
plot(T,YYY(:,14));
grid on;
title(' cytoplasmic IkBa|NFkB');

subplot(3,2,2);
plot(T,YYY(:,15));
grid on;
title(' nuclear IkBa|NFkB');

subplot(3,2,3);
plot(T,YYY(:,5)+YYY(:,6));
grid on;
title('  total phosphorylated IkBa');

subplot(3,2,4);
plot(T,YYY(:,3));
grid on;
title('  Active IKK, (IKKa)');

subplot(3,2,5);
plot(T,YYY(:,11));
grid on;
title('  Free Cytoplasmic IkBa');
xlabel('Time in min')

subplot(3,2,6);
plot(T,YYY(:,12));
grid on;
title('  Free nuclear IkBa');
xlabel('Time in min')

%% Figure 4
figure(4)
set(gcf,'Color',[1,1,1])
plot(T,Bb);
grid on;
title('number of active receptors');
xlabel('Time in min')

%% Figure 5
figure(5)
set(gcf,'Color',[1,1,1])

subplot(2,1,1);
plot(T,YYY(:,3),'k');
grid on;
hold on;
title('Aktive IKK (IKKa)');

subplot(2,1,2);
plot(T,(YYY(:,8) + YYY(:,15))/10^5,'k');
grid on;
hold on;
title('nuclear NF-kB');
xlabel('Time in min')

%% Figure 6
figure(6)
set(gcf,'Color',[1,1,1])

subplot(2,1,1);
plot(T,GGR);
grid on;
title('Activity of Reporter gene');

subplot(2,1,2);
plot(T,YYY(:,17));
grid on;
title('mRNA of Reporter gene');
xlabel('Time in min')
