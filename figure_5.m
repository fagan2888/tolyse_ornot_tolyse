% figure 5
%
% based on experimental data from Lanying
%

beta0=[1 1];

xxx=0:0.01:15;

sty{1}='bo';
sty{2}='r+';
sty{3}='gx';
sty{4}='ms';
sty{5}='cd';

x_ref=[0.7 0.9 1.1 1.3 1.5];
%x_ref=[1.5 1.3 1.1 0.9 0.7];
% cell length
% proportional to V

f_ref=[0.5532    0.7298    0.8304    0.8441    0.8967; 0.5241    0.7150    0.7246    0.7367    0.7595; 0.3701    0.5206    0.6190    0.7604    0.7730; 0.3260    0.4464    0.5988    0.5933    0.6814;    0.2143    0.2497    0.3953    0.4992    0.6631];
% fraction of lysogeny
% taken from Zeng et al 2010 Cell

for i=1:5
    for j=1:5
        xmv_ref((i-1)*5+j)=j./x_ref(i);
        xf_lyg1((i-1)*5+j)=f_ref(i,j);
        xf_lyg2((i-1)*5+j)=f_ref(i,j)^(1/j);
        xf_lyg3((i-1)*5+j)=f_ref(i,1)^(1/j);
        
        mv_ref(i,j)=j./x_ref(i);
        f_lyg1(i,j)=f_ref(i,j);
        f_lyg2(i,j)=f_ref(i,j)^(1/j);
        
    end
end

beta3=nlinfit(xmv_ref,xf_lyg1,@hill_ftn,beta0);
beta4=nlinfit(xmv_ref,xf_lyg2,@hill_ftn,beta0);

yyy3=hill_ftn(beta3,xxx);
yyy4=hill_ftn(beta4,xxx);

%figure('units','normalized','position',[0.02 .1 .96 .7]);
figure('units','normalized','position',[0.05 .1 .5 .7]);
hold on;

%subplot(1,3,1);
hold on;
for i=1:5
    plot(mv_ref(:,i),f_lyg1(:,i),sty{i},'MarkerSize',12,'LineWidth',4);
end
legend('$\mathcal{M}=1$','$\mathcal{M}=2$','$\mathcal{M}=3$','$\mathcal{M}=4$','$\mathcal{M}=5$','Location','Southeast');
h=legend;
set(h,'Interpreter','Latex');
legend('Boxoff');
set(gca,'Box','on');
ax=get(gca,'Position');
ax(1)=ax(1)+0.05; %or wathever
ax(2)=ax(2)+0.05;
ax(3)=ax(3);
ax(4)=ax(4);
set(gca,'Position',ax); 
plot(xxx,yyy3,'k','LineWidth',2);
set(gca,'FontSize',30,'LineWidth',2);
xlim([0 8]);
ylim([0 1]);
%title(tit{i});
xlabel('$\mathcal{M}/V$','Interpreter','Latex');
ylabel('$P_{lysg}$','Interpreter','Latex');
annotation('textbox','String','A','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.01 0.9 0.1 0.1]);


figure('units','normalized','position',[0.05 .1 .5 .7]);
%subplot(1,3,2);
hold on;
for i=1:5
    plot(mv_ref(:,i),f_lyg2(:,i),sty{i},'MarkerSize',12,'LineWidth',4);
end
legend('$\mathcal{M}=1$','$\mathcal{M}=2$','$\mathcal{M}=3$','$\mathcal{M}=4$','$\mathcal{M}=5$','Location','Southeast');
h=legend;
set(h,'Interpreter','Latex');
legend('Boxoff');
set(gca,'Box','on');
ax=get(gca,'Position');
ax(1)=ax(1)+0.1; %or wathever
ax(2)=ax(2)+0.05;
ax(3)=ax(3)-0.05;
ax(4)=ax(4);
set(gca,'Position',ax); 
plot(xxx,yyy4,'k','LineWidth',2);
set(gca,'FontSize',30,'LineWidth',2);
xlim([0 8]);
ylim([0 1]);
%title(tit{i});
xlabel('$\mathcal{M}/V$','Interpreter','Latex');
ylabel('$P_{lysg}^{1/\mathcal{M}}$','Interpreter','Latex');
annotation('textbox','String','C','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.01 0.9 0.1 0.1]);


figure('units','normalized','position',[0.05 .1 .5 .7]);
%subplot(1,3,3);
fit_hill=zeros(100,1);
fit_hill2=zeros(100,1);

rmad=zeros(length(xxx),25);
rmad_min=zeros(25,1);

for h=1:100
    alpha=h/100;
    for i=1:5
        for j=1:5
            xmv_ref((i-1)*5+j)=j^alpha./x_ref(i);
            xf_lyg1((i-1)*5+j)=f_ref(i,j);
        end
    end
    
    bb=nlinfit(xmv_ref,xf_lyg1,@hill_ftn,beta0);
   
    r21=0;
    r22=0;
    r_rma=0;
    
    ymv_ref=hill_ftn(bb,xmv_ref);
    
    y_rma=hill_ftn(bb,xxx);
    
    for i=1:25
        
        for k=1:length(xxx)
            rmad(k,i)=(xf_lyg1(i)-y_rma(k))^2+(xmv_ref(i)-xxx(k))^2;
        end
        
        r_rma=r_rma+min(rmad(:,i));
        
        r21=r21+(xf_lyg1(i)-ymv_ref(i))^2;
        r22=r22+(xf_lyg1(i)-mean(xf_lyg1))^2;
    end
    
    fit_hill(h)=1-r21/r22;
    fit_hill2(h)=1-r_rma/r22
end

[max_r2 max_ind]=max(fit_hill)

%alpha=max_ind/1000;
%alpha=0.6;
alpha=0.5;

hold on;
for i=1:5
    for j=1:5
        xmv_ref2((i-1)*5+j)=j^alpha./x_ref(i);
        xf_lyg1((i-1)*5+j)=f_ref(i,j);
        
        plot(j^alpha./x_ref(i),f_ref(i,j),sty{j},'MarkerSize',12,'LineWidth',4);        
    end
end
legend('$\mathcal{M}=1$','$\mathcal{M}=2$','$\mathcal{M}=3$','$\mathcal{M}=4$','$\mathcal{M}=5$','Location','Southeast');
h=legend;
set(h,'Interpreter','Latex');
legend('Boxoff');

bb=nlinfit(xmv_ref2,xf_lyg1,@hill_ftn,beta0);

yyy5=hill_ftn(bb,xxx);


plot(xxx,yyy5,'k','LineWidth',2);
max_r2
hold off;

ax=get(gca,'Position');
ax(1)=ax(1)+0.1; %or wathever
ax(2)=ax(2)+0.05;
ax(3)=ax(3)-0.05;
ax(4)=ax(4);
set(gca,'Position',ax);
set(gca,'Box','on');
%set(h,'Position',ax); 
%plot(xxx,yyy3,'k','LineWidth',2);
set(gca,'FontSize',30,'LineWidth',2);
xlim([0 4]);
ylim([0 1]);
%title(tit{i});
xlabel('$\mathcal{M}^{\epsilon}/V$','Interpreter','Latex');
ylabel('$P_{lysg}$','Interpreter','Latex');
annotation('textbox','String','D','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.01 0.9 0.1 0.1]);



figure('units','normalized','position',[0.05 .1 .5 .7]);

col{1}='b';
col{2}='r';
col{3}='g';
col{4}='m';
col{5}='c';

x=0.5:0.1:2;
hold on;
for i=1:5
    %y=hill_ftn(beta3,x.*i);
    y=hill_ftn(beta3,x);
    plot(x.*i,y.^(1/i),col{i},'LineWidth',3);
    %plot(x.*i,y.^(i),col{i},'LineWidth',3);
    
    %xf_2(:,i)=f_ref(:,1).^(1/i);
    
    %plot(mv_ref(:,i),f_ref(:,1).^(1/i),sty{i},'MarkerSize',12,'LineWidth',4);
end
legend('$\mathcal{M}=1$','$\mathcal{M}=2$','$\mathcal{M}=3$','$\mathcal{M}=4$','$\mathcal{M}=5$','Location','Southeast');
h=legend;
set(h,'Interpreter','Latex');
legend('Boxoff');

%xf_2l=[xf_2(:,1); xf_2(:,2); xf_2(:,3); xf_2(:,4); xf_2(:,5)];

% b6=nlinfit(xmv_ref,xf_lyg3,@hill_ftn,beta0);
% yyy6=hill_ftn(b6,xxx);
% plot(xxx,yyy6,'k','LineWidth',2);

hold off;
%set(h,'Position',ax); 
%plot(xxx,yyy3,'k','LineWidth',2);
ax=get(gca,'Position');
ax(1)=ax(1)+0.1; %or wathever
ax(2)=ax(2)+0.05;
ax(3)=ax(3)-0.05;
ax(4)=ax(4);
set(gca,'Position',ax); 
set(gca,'FontSize',30,'LineWidth',2);
set(gca,'Box','on');
xlim([0 8]);
ylim([0 1]);
%title(tit{i});
xlabel('$\mathcal{M}$/V','Interpreter','Latex');
ylabel('$P_{lysg}^{1/\mathcal{M}}$','Interpreter','Latex');
annotation('textbox','String','B','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.01 0.9 0.1 0.1]);

