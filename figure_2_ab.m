% figure 2(a,b)
%
% phase plane ODE solutions
%
%
% 1- asymptotically divergent
% 2- trnasiently divergent

clear;

load data_para_series2;
para_sto=para_series(:,215);
% transiently div

load data_para_series3;
para_det=para_series(:,99);%38);
% asymptotically div

it_max=size(feat_series,2);

x=1:it_max;

para_det1=[1 transpose(para_det)];
% moi=1

ini_cond=zeros(9,1);

% figure 2(A)
hold on;

for i=0:50:200;
   
    for j=[0 300];
        
        ini_cond=zeros(9,1);
        ini_cond(6)=j;
        ini_cond(9)=i;
        [t x]=ode45(@model_final,[0 1000],ini_cond,[],para_det1);
        % change initial conditions
        
        plot(x(:,6),x(:,9),'k','LineWidth',1);
        plot(x(end,6),x(end,9),'k.','Markersize',40);
        
        for k=[80]
            arrow([x(k-1,6) x(k-1,9)],[x(k,6) x(k,9)],12,'FaceColor','k','Edgecolor','k');
        end
        % draw little arrows on the solution
        
        clear t x;
    end
 end

for i=0:30:300;
   
    for j=[0 200];
        
        ini_cond=zeros(9,1);
        ini_cond(6)=i;
        ini_cond(9)=j;
        [t x]=ode45(@model_final,[0 1000],ini_cond,[],para_det1);

        plot(x(:,6),x(:,9),'k','LineWidth',1);
        plot(x(end,6),x(end,9),'k.','Markersize',40);

        for k=[80]
            arrow([x(k-1,6) x(k-1,9)],[x(k,6) x(k,9)],12,'FaceColor','k','Edgecolor','k');
        end

        clear t x;
    end
end

set(gca,'box','on');
set(gca,'FontSize',30,'LineWidth',2);
ylim([0 200]);
xlim([0 300]);
ylabel('Q (nM)');
xlabel('cI (nM)');
title('Asymptotically divergent, $\mathcal M=1$','Interpreter','latex');
annotation('textbox','String','A','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.02 0.88 0.1 0.1]);
hold off; 



% figure 2(B)
figure;

para_det1=[1 transpose(para_sto)];
 
ini_cond=zeros(9,1);
 
hold on;

for i=0:50:200;
   
    for j=[0 300];
        
        ini_cond=zeros(9,1);
        ini_cond(6)=j;
        ini_cond(9)=i;
        [t x]=ode45(@model_final,[0 1000],ini_cond,[],para_det1);

        plot(x(:,6),x(:,9),'k','LineWidth',1);
        plot(x(end,6),x(end,9),'k.','Markersize',40);
        
        for k=[60]
            arrow([x(k-1,6) x(k-1,9)],[x(k,6) x(k,9)],12,'FaceColor','k','Edgecolor','k');
        end
        clear t x;
    end
 end

for i=0:30:300;
   
    for j=[0 200];
        
        ini_cond=zeros(9,1);
        ini_cond(6)=i;
        ini_cond(9)=j;
        [t x]=ode45(@model_final,[0 1000],ini_cond,[],para_det1);

        plot(x(:,6),x(:,9),'k','LineWidth',1);
        plot(x(end,6),x(end,9),'k.','Markersize',40);

        for k=[60]
            arrow([x(k-1,6) x(k-1,9)],[x(k,6) x(k,9)],12,'FaceColor','k','Edgecolor','k');
        end

        clear t x;
    end
end
set(gca,'box','on');
set(gca,'FontSize',30,'LineWidth',2);
ylim([0 200]);
xlim([0 300]);
ylabel('Q (nM)');
xlabel('cI (nM)');
title('Transiently divergent, $\mathcal M=1$','Interpreter','latex');
annotation('textbox','String','D','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.02 0.88 0.1 0.1]);
hold off; 
 
