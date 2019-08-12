% %##########################################################################
% %##################      PLOTS ALL CELLS TRAJECTORIES     #################
% %##########################################################################


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

% %##########################################################################

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

%##########################################################################

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


%##########################################################################

figure(4)
set(gcf,'Color',[1,1,1])
plot(T,XB+0.001);
grid on;
title('Number of active receptors');
xlabel('Time (min)')

%##########################################################################

figure(5)
subplot(2,1,1);
plot(T,XXX(:,:,8) + XXX(:,:,15));
%plot(T,(XXX(:,:,8) + XXX(:,:,15))/10^5);
grid on;
hold on;
title('Nuclear NF-kB');
ylabel('Number of Molecules');

subplot(2,1,2);
plot(T,XXX(:,:,3));
grid on;
hold on;
title('Active IKK (IKKa)');
xlabel('Time (min)')
ylabel('Number of Molecules');

%##########################################################################

figure(6)
set(gcf,'Color',[1,1,1])

subplot(2,1,1);
plot(T,XGR);
grid on;
title('Activity of reporter gene');
ylabel('States');

subplot(2,1,2);
plot(T,XXX(:,:,17));
grid on;
title('Reporter mRNA');
ylabel('Number of Molecules');

