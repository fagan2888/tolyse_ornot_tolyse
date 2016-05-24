% create figure 3
%
% indiviaul realization of stochstic dynamics
%
% transiently divergent: MOI=1, 2
%
clear;

th_ci=100;
th_q=100;

% transiently divergent
load para_sto_series3;

% moi=1
hold on;
for i=1:50
    
    t_ind=ceil(4*min(ciexit(th_ci,i,1),qexit(th_q,i,1)))+1; 
    % time for first passage
    
    xx=zeros(t_ind,2);
    
    % only draw before first passage
    for j=2:t_ind-1
        xx(j,1)=ci_dyn(j-1,i,1);
        xx(j,2)=q_dyn(j-1,i,1);
    end
    
    % if lysogenic
    if ciexit(th_ci,i,1)<qexit(th_q,i,1)
        xx(t_ind,1)=th_ci;
        xx(t_ind,2)=q_dyn(t_ind-1,i,1);
        
        plot(xx(:,1),xx(:,2),'b');
        
    end
    
    % if lytic
    if ciexit(th_ci,i,1)>qexit(th_q,i,1)
        xx(t_ind,1)=ci_dyn(t_ind-1,i,1);
        xx(t_ind,2)=th_q;
        
        plot(xx(:,1),xx(:,2),'r');
        
    end
end
hold off;     
set(gca,'Fontsize',25,'LineWidth',2);
title('$\mathcal M=1$','Interpreter','latex');
% transiently driven
xlabel('cI (nM)');
ylabel('Q (nM)');
set(gca,'box','on');
set(gca, 'XTick', 0:20:100);
set(gca, 'YTick', 0:20:100);
annotation('textbox','String','A','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0 0.89 0.1 0.1]);
annotation('textbox','String','Lysogeny 10%','Fontsize',20,'LineStyle','none','FontWeight','bold','Position',[0.64 0.33 0.3 0.1],'color','b');
annotation('textbox','String','Lysis 90%','Fontsize',20,'LineStyle','none','FontWeight','bold','Position',[0.53 0.81 0.3 0.1],'Color','r');


% moi=2
figure;

hold on;
for i=1:50
    
    t_ind=ceil(4*min(ciexit(th_ci,i,2),qexit(th_q,i,2)))+1; 

    xx=zeros(t_ind,2);
    
    for j=2:t_ind-1
        xx(j,1)=ci_dyn(j-1,i,2);
        xx(j,2)=q_dyn(j-1,i,2);
    end
    
    if ciexit(th_ci,i,2)<qexit(th_q,i,2)
        xx(t_ind,1)=th_ci;
        xx(t_ind,2)=q_dyn(t_ind-1,i,2);
        
        plot(xx(:,1),xx(:,2),'b');
        
    end
    
    if ciexit(th_ci,i,2)>qexit(th_q,i,2)
        xx(t_ind,1)=ci_dyn(t_ind-1,i,2);
        xx(t_ind,2)=th_q;
        
        plot(xx(:,1),xx(:,2),'r');
        
    end
end
hold off;     
set(gca,'Fontsize',25,'Linewidth',2);
title('$\mathcal M=2$','Interpreter','latex');
xlabel('cI (nM)');
ylabel('Q (nM)');
set(gca,'box','on');
set(gca, 'XTick', 0:20:100);
set(gca, 'YTick', 0:20:100);
annotation('textbox','String','B','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0 0.89 0.1 0.1]);
annotation('textbox','String','Lysogeny 45%','Fontsize',20,'LineStyle','none','FontWeight','bold','Position',[0.6 0.25 0.3 0.1],'color','b');
annotation('textbox','String','Lysis 55%','Fontsize',20,'LineStyle','none','FontWeight','bold','Position',[0.18 0.88 0.3 0.1],'Color','r');

