% figure S1
%
% create heatmap of transcription probability
% 
% based on Shea & Ackers
%
x=0:0.005:3;
y=0:0.005:3;

[xx yy]=meshgrid(10.^x,10.^y);

alpha_x=0.06;
% shea 85
beta_x=0.66;
% shea 85
delta_x=0.9;
% arkin 98
alpha_y=0.84;
alpha_z=0.8;
alpha_q=0.75;
delta_aq=2;
alpha_n=0.6;

gamma_x=0.01;
gamma_y=0.06;%4;
gamma_z=0.10;
gamma_q=0.01;%0.01;
gamma_n=0.02;

gamma_m=0.1;

c_d_x=1/20;
% Darling 2000
c_d_y=1/(1/5.8);%1/100;
% Jana et al. JMB 97
c_d_z=1/20;

c_t_z=1/20;%1/10;

c_p_aq=1/5;

sigma=0.5;

zeta=0.1;

cx1=[-1+sqrt(1+8.*c_d_x.*10.^x)]./4./c_d_x;
cy1=[-1+sqrt(1+8.*c_d_y.*10.^y)]./4./c_d_y;
% equilibrium monomer concentration

cx2=cx1.^2.*c_d_x;
cy2=cy1.^2.*c_d_y;

cz1=zeros(length(x));

for i=1:length(x)
    while cz1(i)+2*c_d_z*cz1(i)^2+4*c_t_z*c_d_z^2*cz1(i)^4<10.^x(i)
        cz1(i)=cz1(i)+0.001;
    end
end
% fine monomer concentration

cz4=c_t_z*c_d_z^2.*(cz1-0.001).^4;


zz1=zeros(length(xx),length(yy));
zz2=zeros(length(xx),length(yy));
zz3=zeros(length(xx),length(yy));
zz4=zeros(length(xx),length(yy));

zz5=zeros(xx,1);

for i=1:length(x)
    for j=1:length(y)
        [zz1(i,j) zz2(i,j) zz3(i,j)]=pr_prm(cx2(i),cy2(j));
        zz4(i,j)=pl(cx2(i),cy2(j));
    end
    
    zz5(i)=pre(cz4(i));
end

zz1(end,1)=0;
zz1(1,end)=1;
zz2(end,1)=0;
zz2(1,end)=1;
zz3(end,1)=0;
zz3(1,end)=1;
zz4(end,1)=0;
zz4(1,end)=1;
%zz5(end,1)=0;
%zz5(1,end)=1;


[c,ch]=contourf(xx,yy,transpose(zz1),100);
set(ch,'edgecolor','none');
set(gca,'FontSize',30,'Linewidth',2);
set(gca,'Xscale','log');
set(gca,'Yscale','log');
set(gca, 'XTick',[1 10 100 1000]);
xlabel('CI (nM)');
ylabel('CRO (nM)');
title('f^{basal}_{RM}');
xlim([1 1000]);
ylim([1 1000]);
colorbar;
annotation('textbox','String','A','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.02 0.88 0.1 0.1]);

figure;
[ch,ch]=contourf(xx,yy,transpose(zz2),100);
set(ch,'edgecolor','none');
set(gca,'FontSize',30,'Linewidth',2);
set(gca,'Xscale','log');
set(gca,'Yscale','log');
set(gca, 'XTick',[1 10 100 1000]);
xlabel('CI (nM)');
ylabel('CRO (nM)');
title('f^{activated}_{RM}');
xlim([1 1000]);
ylim([1 1000]);
colorbar;
annotation('textbox','String','B','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.02 0.88 0.1 0.1]);

figure;
[ch, ch]=contourf(xx,yy,transpose(zz3),100);
set(ch,'edgecolor','none');
set(gca,'FontSize',30,'Linewidth',2);
set(gca,'Xscale','log');
set(gca,'Yscale','log');
set(gca, 'XTick',[1 10 100 1000]);
xlabel('CI (nM)');
ylabel('CRO (nM)');
title('f_{R}');
xlim([1 1000]);
ylim([1 1000]);
colorbar;
annotation('textbox','String','C','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.02 0.88 0.1 0.1]);

figure;
plot(10.^x,zz5,'LineWidth',3);
set(gca,'FontSize',30,'Linewidth',2);
set(gca,'Xscale','log');
set(gca, 'XTick',[1 10 100 1000]);
xlabel('CII (nM)');
ylabel('f_{RE}');
annotation('textbox','String','D','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.02 0.88 0.1 0.1]);
xlim([1 1000]);


