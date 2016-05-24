clear;

load golding_data;

plot(m_v(:,1),frac_lysg(:,1),'x',m_v(:,2),frac_lysg(:,2),'o',m_v(:,3),frac_lysg(:,3),'*',m_v(:,4),frac_lysg(:,4),'s',m_v(:,5),frac_lysg(:,5),'^','MarkerSize',15,'LineWidth',4);
set(gca,'FontSize',30,'LineWidth',2);
xlabel('$\mathcal M / V$','Interpreter','latex');
ylabel('Fraction of lysogeny');
h=legend('$\mathcal M=1$','$\mathcal M=2$','$\mathcal M=3$','$\mathcal M=4$','$\mathcal M=5$');
set(h,'Interpreter','latex','Location','Southeast');
legend('boxoff');
title('Zeng et. al. (2010)');
ylim([0 1]);
xlim([0.5 5]);
annotation('textbox','String','A','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.01 0.89 0.1 0.1]);

% fail=20%

for i=1:5
% observed moi    
    for j=1:6
    % actual moi    
    % j= m-1     

        mprob(i,j)=factorial(i)/factorial(max(i-j+1,0))/factorial(j-1)*(0.8)^(j-1)*(0.2)^(i-j+1)*[i+1>=j];
        % probability of MOI=j-1 when observed MOI=i;
    end
end

for i=1:5
    
    frac_lysg2(i,1)=(frac_lysg(i,1))./mprob(1,2);
    frac_lysg2(i,2)=(frac_lysg(i,2)-mprob(2,2)*frac_lysg2(i,1))./mprob(2,3);
    frac_lysg2(i,3)=(frac_lysg(i,3)-mprob(3,2)*frac_lysg2(i,1)-mprob(3,3)*frac_lysg2(i,2))./mprob(3,4);
    frac_lysg2(i,4)=(frac_lysg(i,4)-mprob(4,2)*frac_lysg2(i,1)-mprob(4,3)*frac_lysg2(i,2)-mprob(4,4)*frac_lysg2(i,3))./mprob(4,5);    
    frac_lysg2(i,5)=(frac_lysg(i,5)-mprob(5,2)*frac_lysg2(i,1)-mprob(5,3)*frac_lysg2(i,2)-mprob(5,4)*frac_lysg2(i,3)-mprob(5,5)*frac_lysg2(i,4))./mprob(5,6);

end

figure;

plot(m_v(:,1),frac_lysg2(:,1),'x',m_v(:,2),frac_lysg2(:,2),'o',m_v(:,3),frac_lysg2(:,3),'*',m_v(:,4),frac_lysg2(:,4),'s',m_v(:,5),frac_lysg2(:,5),'^','MarkerSize',15,'LineWidth',4);
set(gca,'FontSize',30,'LineWidth',2);
xlabel('$\mathcal M / V$','Interpreter','latex');
ylabel('Fraction of lysogeny');
h=legend('$\mathcal M=1$','$\mathcal M=2$','$\mathcal M=3$','$\mathcal M=4$','$\mathcal M=5$');
set(h,'Interpreter','latex','Location','Southeast');
legend('boxoff');
title('Zeng et. al. (2010) with 20% Failure');
ylim([0 1]);
xlim([0.5 5]);
annotation('textbox','String','B','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.01 0.89 0.1 0.1]);


figure;

%n_file=7;
n_file=1;

n_step=10;

p_lysg=zeros(n_step,5);
p_lyss=zeros(n_step,5);
p_und=zeros(n_step,5);

for iii=1:n_file

%loadfile=strcat('mv_sto',num2str(iii));    
%load(loadfile);
load mv_det_5min_2;

%th_ci=200;
%th_q=150;

th_ci=100;
th_q=100;
% max for sto

%th_ci=100;
%th_q=100;
% max for det

[a1 a2 a3]=size(ciexit);

mv=zeros(n_step,5);

for i=1:a3
    mv(:,i)=i.*(0.5+1.5/n_step/2:1.5/n_step:2);
end

%p_error=zeros(n_step,5,4);

for i=1:a2    
    vol_ind=ceil((1/vol(i)-0.5)*n_step/1.5);
    
    for j=1:a3 
        if ciexit(th_ci,i,j)<qexit(th_q,i,j)
            p_lysg(vol_ind,j)=p_lysg(vol_ind,j)+1;
        elseif ciexit(th_ci,i,j)>qexit(th_q,i,j)
            p_lyss(vol_ind,j)=p_lyss(vol_ind,j)+1;
        else
            p_und(vol_ind,j)=p_und(vol_ind,j)+1;
        end
        
%         if ciexit(th_ci-5,i,j)<qexit(th_q,i,j)
%             p_error(vol_ind,j,1)=p_error(vol_ind,j,1)+1;
%         end
%         
%         if ciexit(th_ci+5,i,j)<qexit(th_q,i,j)
%             p_error(vol_ind,j,2)=p_error(vol_ind,j,2)+1;
%         end
%         
%         if ciexit(th_ci,i,j)<qexit(th_q-5,i,j)
%             p_error(vol_ind,j,3)=p_error(vol_ind,j,3)+1;
%         end
%         
%         if ciexit(th_ci,i,j)<qexit(th_q+5,i,j)
%             p_error(vol_ind,j,4)=p_error(vol_ind,j,4)+1;
%         end
        
    end
end

end

for i=1:n_step
    for j=1:a3

        summ=p_lysg(i,j)+p_lyss(i,j)+p_und(i,j);
        
        p_lysg(i,j)=p_lysg(i,j)/summ;
        p_lyss(i,j)=p_lyss(i,j)/summ;
        p_und(i,j)=p_und(i,j)/summ;
        
%         for k=1:4
%             p_error(i,j,k)=p_error(i,j,k)/summ;
%         end
        
    end
end

% er_l=zeros(n_step,5);
% er_h=zeros(n_step,5);
% 
% for i=1:n_step
%     for j=1:5
%         
%         er_l(i,j)=abs( min ( min(p_error(i,j,:)-p_lysg(i,j)),0));
%         er_h(i,j)=max(max(p_error(i,j,:)-p_lysg(i,j)),0);
%         
%     end
% end

% errorbar(mv,p_lysg,er_l,er_h,'x','Markersize',15,'LineWidth',3);
plot(mv(:,1),p_lysg(:,1),'x',mv(:,2),p_lysg(:,2),'o',mv(:,3),p_lysg(:,3),'*',mv(:,4),p_lysg(:,4),'s',mv(:,5),p_lysg(:,5),'^','MarkerSize',15,'LineWidth',4);
set(gca,'FontSize',30,'LineWidth',2);
xlabel('$\mathcal M / V$','Interpreter','latex');
ylabel('Fraction of lysogeny');
h=legend('$\mathcal M=1$','$\mathcal M=2$','$\mathcal M=3$','$\mathcal M=4$','$\mathcal M=5$');
set(h,'Interpreter','latex','Location','Southeast');
legend('boxoff');
%title('Transiently driven');
title('Asymptotically divergent 5min')
ylim([0 1]);
xlim([0.5 5]);
annotation('textbox','String','C','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.01 0.89 0.1 0.1]);



figure;

clear;

n_file=5;

n_step=10;

p_lysg=zeros(n_step,5);
p_lyss=zeros(n_step,5);
p_und=zeros(n_step,5);

for iii=2:n_file%1:n_file
    
loadfile=strcat('mv_det_15min_',num2str(iii));      
load(loadfile);
%load mv_det1;

%load mv_det_15min_2;

%th_ci=200;
%th_q=150;

th_ci=100;
th_q=100;
% max for sto

%th_ci=100;
%th_q=100;
% max for det

[a1 a2 a3]=size(ciexit);

mv=zeros(n_step,5);

for i=1:a3
    mv(:,i)=i.*(0.5+1.5/n_step/2:1.5/n_step:2);
end

%p_error=zeros(n_step,5,4);

for i=1:a2    
    vol_ind=ceil((1/vol(i)-0.5)*n_step/1.5);
    
    for j=1:a3 
        if ciexit(th_ci,i,j)<qexit(th_q,i,j)
            p_lysg(vol_ind,j)=p_lysg(vol_ind,j)+1;
        elseif ciexit(th_ci,i,j)>qexit(th_q,i,j)
            p_lyss(vol_ind,j)=p_lyss(vol_ind,j)+1;
        else
            p_und(vol_ind,j)=p_und(vol_ind,j)+1;
        end
        
%         if ciexit(th_ci-5,i,j)<qexit(th_q,i,j)
%             p_error(vol_ind,j,1)=p_error(vol_ind,j,1)+1;
%         end
%         
%         if ciexit(th_ci+5,i,j)<qexit(th_q,i,j)
%             p_error(vol_ind,j,2)=p_error(vol_ind,j,2)+1;
%         end
%         
%         if ciexit(th_ci,i,j)<qexit(th_q-5,i,j)
%             p_error(vol_ind,j,3)=p_error(vol_ind,j,3)+1;
%         end
%         
%         if ciexit(th_ci,i,j)<qexit(th_q+5,i,j)
%             p_error(vol_ind,j,4)=p_error(vol_ind,j,4)+1;
%         end
        
    end
end

end

for i=1:n_step
    for j=1:a3

        summ=p_lysg(i,j)+p_lyss(i,j)+p_und(i,j);
        
        p_lysg(i,j)=p_lysg(i,j)/summ;
        p_lyss(i,j)=p_lyss(i,j)/summ;
        p_und(i,j)=p_und(i,j)/summ;
        
%         for k=1:4
%             p_error(i,j,k)=p_error(i,j,k)/summ;
%         end
        
    end
end

% er_l=zeros(n_step,5);
% er_h=zeros(n_step,5);
% 
% for i=1:n_step
%     for j=1:5
%         
%         er_l(i,j)=abs( min ( min(p_error(i,j,:)-p_lysg(i,j)),0));
%         er_h(i,j)=max(max(p_error(i,j,:)-p_lysg(i,j)),0);
%         
%     end
% end

% errorbar(mv,p_lysg,er_l,er_h,'x','Markersize',15,'LineWidth',3);
plot(mv(:,1),p_lysg(:,1),'x',mv(:,2),p_lysg(:,2),'o',mv(:,3),p_lysg(:,3),'*',mv(:,4),p_lysg(:,4),'s',mv(:,5),p_lysg(:,5),'^','MarkerSize',15,'LineWidth',4);
set(gca,'FontSize',30,'LineWidth',2);
xlabel('$\mathcal M / V$','Interpreter','latex');
ylabel('Fraction of lysogeny');
h=legend('$\mathcal M=1$','$\mathcal M=2$','$\mathcal M=3$','$\mathcal M=4$','$\mathcal M=5$');
set(h,'Interpreter','latex','Location','Southeast');
legend('boxoff');
title('Asymptotically divergent 15min');
ylim([0 1]);
xlim([0.5 5]);
annotation('textbox','String','D','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.01 0.89 0.1 0.1]);

