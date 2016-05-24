% figure 2(c,d,e,f)
%
% form two paramer sets create figures of ODE solutions
%
% 1- asymprotically divergent
% 2- stochastically divergent

clear;

load data_para_series2;
para_sto=para_series(:,215);
% asymptotially divergent

load data_para_series3;
para_det=para_series(:,99);%38);
% transitnely divergent

it_max=size(feat_series,2);

x=1:it_max;

para_det1=[1 transpose(para_det)];
para_det2=[2 transpose(para_det)];
para_det3=[3 transpose(para_det)];
% different MOI 

ini_cond=zeros(9,1);
% start from various initial conditions

[t1 x1]=ode45(@model_final,[0 1000],ini_cond,[],para_det1);
[t2 x2]=ode45(@model_final,[0 1000],ini_cond,[],para_det2);
[t3 x3]=ode45(@model_final,[0 1000],ini_cond,[],para_det3);

% figure 2(B)
% CI-Q phase plane
figure;
hold on;
plot(x1(:,6),x1(:,9),x2(:,6),x2(:,9),'--',x3(:,6),x3(:,9),'-.','LineWidth',3);
plot([100 100 0],[0 100 100],'k:','LineWidth',3);
set(gca,'box','on');
legend('$\mathcal M=1$','$\mathcal M=2$','$\mathcal M=3$','Threshold');
h=legend;
set(h,'Interpreter','latex');
legend('boxoff');
plot(x1(end,6),x1(end,9),'.',x2(end,6),x2(end,9),'.',x3(end,6),x3(end,9),'.','Markersize',45);
%plot(t1,x1(:,9),t2,x2(:,9),'--',t3,x3(:,9),'-.','LineWidth',3);
set(gca,'Fontsize',30,'LineWidth',2);
ylim([0 170]);
xlim([0 550]);
ylabel('Q (nM)');
xlabel('cI (nM)');
legend('boxoff');
title('Asymptotically divergent','Interpreter','latex');
annotation('textbox','String','B','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.02 0.88 0.1 0.1]);
arrow([x1(floor(end/8)-1,6),x1(floor(end/8)-1,9)],[x1(floor(end/8),6),x1(floor(end/8),9)],25,'FaceColor','b','Edgecolor','b');
arrow([x1(floor(end/4)-1,6),x1(floor(end/4)-1,9)],[x1(floor(end/4),6),x1(floor(end/4),9)],25,'FaceColor','b','Edgecolor','b');
arrow([x2(floor(3*end/4)-1,6),x2(floor(3*end/4)-1,9)],[x2(floor(3*end/4),6),x2(floor(3*end/4),9)],25,'FaceColor','g','Edgecolor','g');
arrow([x2(floor(end/1.6)-1,6),x2(floor(end/1.6)-1,9)],[x2(floor(end/1.6),6),x2(floor(end/1.6),9)],25,'FaceColor','g','Edgecolor','g');
arrow([x3(floor(3*end/4)-1,6),x3(floor(3*end/4)-1,9)],[x3(floor(3*end/4),6),x3(floor(3*end/4),9)],25,'FaceColor','r','Edgecolor','r');
arrow([x3(floor(end/1.6)-1,6),x3(floor(end/1.6)-1,9)],[x3(floor(end/1.6),6),x3(floor(end/1.6),9)],25,'FaceColor','r','Edgecolor','r');
hold off;

% figure 2(C)
% Q dynamics
figure;
plot(t1,x1(:,9),t2,x2(:,9),'--',t3,x3(:,9),'-.',[0 500],[100 100],'k:','LineWidth',3);
set(gca,'Fontsize',30,'LineWidth',2);
ylim([0 170]);
xlim([0 500]);
ylabel('Q (nM)');
xlabel('t (min)');
legend('boxoff');
title('Asymptotically divergent','Interpreter','latex');
annotation('textbox','String','C','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.02 0.88 0.1 0.1]);



para_sto1=[1 transpose(para_sto)];
para_sto2=[2 transpose(para_sto)];
para_sto3=[3 transpose(para_sto)];

[ts1 xs1]=ode45(@model_final,[0 1000],ini_cond,[],para_sto1);
[ts2 xs2]=ode45(@model_final,[0 1000],ini_cond,[],para_sto2);
[ts3 xs3]=ode45(@model_final,[0 1000],ini_cond,[],para_sto3);

% figure 2(E)
figure;
hold on;
plot(xs1(:,6),xs1(:,9),xs2(:,6),xs2(:,9),'--',xs3(:,6),xs3(:,9),'-.','LineWidth',3);
plot([100 100 0],[0 100 100],'k:','LineWidth',3);
set(gca,'box','on');
legend('$\mathcal M=1$','$\mathcal M=2$','$\mathcal M=3$','Threshold');
h=legend;
set(h,'Interpreter','latex');
legend('boxoff');
plot(xs1(end,6),xs1(end,9),'.',xs2(end,6),xs2(end,9),'.',xs3(end,6),xs3(end,9),'.','Markersize',45);
%plot(ts1,xs1(:,9),ts2,xs2(:,9),'--',ts3,xs3(:,9),'-.','LineWidth',3);
set(gca,'Fontsize',30,'LineWidth',2); 
ylim([0 170]);
xlim([0 550]);
ylabel('Q (nM)');
xlabel('cI (nM)');
title('Transiently divergent','Interpreter','latex');
annotation('textbox','String','E','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.02 0.88 0.1 0.1]);
arrow([xs1(floor(7*end/8)-1,6),xs1(floor(7*end/8)-1,9)],[xs1(floor(7*end/8),6),xs1(floor(7*end/8),9)],25,'FaceColor','b','Edgecolor','b');
arrow([xs1(floor(5*end/8)-1,6),xs1(floor(5*end/8)-1,9)],[xs1(floor(5*end/8),6),xs1(floor(5*end/8),9)],25,'FaceColor','b','Edgecolor','b');
arrow([xs2(floor(3*end/8)-1,6),xs2(floor(3*end/8)-1,9)],[xs2(floor(3*end/8),6),xs2(floor(3*end/8),9)],25,'FaceColor','g','Edgecolor','g');
arrow([xs2(floor(end/2)-1,6),xs2(floor(end/2)-1,9)],[xs2(floor(end/2),6),xs2(floor(end/2),9)],25,'FaceColor','g','Edgecolor','g');
arrow([xs3(floor(3*end/8)-1,6),xs3(floor(3*end/8)-1,9)],[xs3(floor(3*end/8),6),xs3(floor(3*end/8),9)],25,'FaceColor','r','Edgecolor','r');
arrow([xs3(floor(end/2)-1,6),xs3(floor(end/2)-1,9)],[xs3(floor(end/2),6),xs3(floor(end/2),9)],25,'FaceColor','r','Edgecolor','r');
hold off;

% figure 2(F)
figure;
plot(ts1,xs1(:,9),ts2,xs2(:,9),'--',ts3,xs3(:,9),'-.',[0 500],[100 100],'k:','LineWidth',3);
set(gca,'Fontsize',30,'LineWidth',2); 
ylim([0 170]);
xlim([0 500]);
ylabel('Q (nM)');
xlabel('t (min)');
title('Transiently divergent','Interpreter','latex');
annotation('textbox','String','F','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.02 0.88 0.1 0.1]);
 
