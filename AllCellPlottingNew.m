%% Preamble
% Program: AllCellPlottingNew
% Author: Jonathan Myers (Modified From S. Tay et al. 2010 Nature)
% Date: July 9th, 2019
% Purpose: To make the Tay MATLAB AllCellPlotting script more streamlined and
% easier to understand and modify.
% Arguments: None. Uses variables T, XXX, XGa, XG, XGR, XB that should be
% in Workspace memory.
% Calls: None.
% Returns: None. Generates Plots of Molecular Concentrations and Gene States.

%% Figure 1
figure(1)
set(gcf,'Color',[1,1,1])

subplot(4,2,1);
plot(T,XXX(:,:,16));
grid on;
title('TNF activity');

subplot(4,2,2);
plot(T,XXX(:,:,1));
grid on;
title('Active IKKK, (IKKKp)');

subplot(4,2,3);
plot(T,KNN-XXX(:,:,2)-XXX(:,:,3)-XXX(:,:,4));
grid on;
title('Intermediate Inactive IKK, (IKKii)');

subplot(4,2,4);
plot(T,XXX(:,:,2));
grid on;
title('Neutral IKK, (IKKn)');

subplot(4,2,5);
plot(T,XXX(:,:,3));
grid on;
title('Active IKK, (IKKa)');

subplot(4,2,6);
plot(T,XXX(:,:,4));
grid on;
title('[Inactive IKK, (IKKi)');

subplot(4,2,7);
plot(T,XXX(:,:,5));
grid on;
title('Phosphorylated IkBa');
xlabel('Time (min)')

subplot(4,2,8);
plot(T,XXX(:,:,6));
grid on;
title('Phosphorylated IkBa|NFkB');
xlabel('Time (min)')

%% Figure 2
figure(2)
set(gcf,'Color',[1,1,1])

subplot(4,2,1);
plot(T,XXX(:,:,3));
grid on;
title('Active IKK, (IKKa)');

subplot(4,2,2);
plot(T,XXX(:,:,12)+XXX(:,:,15)+XXX(:,:,11)+XXX(:,:,14)+XXX(:,:,5)+XXX(:,:,6));
grid on;
title(' Total IkBa');

subplot(4,2,3)
plot(T,XXX(:,:,8));
grid on;
title('Free nuclear NF-kB');

subplot(4,2,4);
plot(T,XGa);
grid on;
title('Activity of IkBa gene');

subplot(4,2,5);
plot(T,XG);
grid on;
title('Activity of A20 gene'); 

subplot(4,2,6);
plot(T,XXX(:,:,13));
grid on;
title('IkBa mRNA transcript');

subplot(4,2,7);
plot(T,XXX(:,:,10));
grid on;
title('A20 mRNA transcript');
xlabel('Time (min)')

subplot(4,2,8);
plot(T,XXX(:,:,9));
grid on;
title('A20 protein');
xlabel('Time (min)')

%% Figure 3
figure(3)
set(gcf,'Color',[1,1,1])

subplot(3,2,1);
plot(T,XXX(:,:,14));
grid on;
title('Cytoplasmic IkBa|NFkB');

subplot(3,2,2);
plot(T,XXX(:,:,15));
grid on;
title('Nuclear IkBa|NFkB');

subplot(3,2,3);
plot(T,XXX(:,:,5)+XXX(:,:,6));
grid on;
title('Total phosphorylated IkBa');
xlabel('Time (hours)')

subplot(3,2,4);
plot(T,XXX(:,:,3));
grid on;
title('Active IKK, (IKKa)');

subplot(3,2,5);
plot(T,XXX(:,:,11));
grid on;
title(' Free Cytoplasmic IkBa');
xlabel('Time (min)')

subplot(3,2,6);
plot(T,XXX(:,:,12));
grid on;
title('Free nuclear IkBa');
xlabel('Time (min)')

%% Figure 4
figure(4)
set(gcf,'Color',[1,1,1])
plot(T,XB+0.001);
grid on;
title('Number of active receptors');
xlabel('Time (min)')

%% Figure 5
figure(5)
set(gcf,'Color',[1,1,1])

subplot(2,1,1);
plot(T,XXX(:,:,8) + XXX(:,:,15),'LineWidth',1);
grid off;
hold on;
title('Nuclear NF-kB','FontSize',10);
ylabel('Number of Molecules','FontSize',10);
text(0.95,0.90,'A','FontSize',15,'Units','normalized');

subplot(2,1,2);
plot(T,XXX(:,:,3),'LineWidth',1);
grid off;
hold on;
title('Active IKK (IKKa)','FontSize',10);
xlabel('Time (min)','FontSize',12);
ylabel('Number of Molecules','FontSize',10);
text(0.95,0.90,'B','FontSize',15,'Units','normalized');

%% Figure 6
figure(6)
set(gcf,'Color',[1,1,1])

subplot(2,1,1);
plot(T,XGR,'LineWidth',1);
grid off;
title('Activity of reporter gene','FontSize',10);
ylabel('States','FontSize',10);
text(0.95,0.90,'A','FontSize',15,'Units','normalized');

subplot(2,1,2);
plot(T,XXX(:,:,17),'LineWidth',1);
grid off;
title('Reporter mRNA','FontSize',10);
ylabel('Number of Molecules','FontSize',10);
xlabel('Time (min)','FontSize',10);
text(0.95,0.90,'B','FontSize',15,'Units','normalized');

%% Abstract Figure
% figure(7)
% set(gcf,'Color',[1,1,1])
% 
% subplot(2,1,1);
% plot(T,XXX(:,:,8) + XXX(:,:,15),'LineWidth',1);
% grid off;
% hold on;
% title('Translocation of Transcription Factors','FontSize',12);
% xlabel('Time (min)','FontSize',10);
% ylabel('Number of Molecules','FontSize',10);
