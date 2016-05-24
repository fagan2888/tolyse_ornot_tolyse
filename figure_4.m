% figure 4
%
% (A) asymprotically
% (B) transiently 
%
clear;

n_file=7;

n_step=10;

p_lysg=zeros(n_step,5);
p_lyss=zeros(n_step,5);
p_und=zeros(n_step,5);

for iii=1:n_file

% data file
%
% varying V as wells as MOI
loadfile=strcat('mv_sto',num2str(iii));    
% transiently div

load(loadfile);

th_ci=100;
th_q=100;
% max for sto

[a1 a2 a3]=size(ciexit);

mv=zeros(n_step,5);

for i=1:a3
    mv(:,i)=i.*(0.5+1.5/n_step/2:1.5/n_step:2);
end
% M/V

for i=1:a2    
    vol_ind=ceil((1/vol(i)-0.5)*n_step/1.5);
    
    % number of cells at each MOI that goes lysogeny and lysis
    for j=1:a3 
        if ciexit(th_ci,i,j)<qexit(th_q,i,j)
            p_lysg(vol_ind,j)=p_lysg(vol_ind,j)+1;
        elseif ciexit(th_ci,i,j)>qexit(th_q,i,j)
            p_lyss(vol_ind,j)=p_lyss(vol_ind,j)+1;
        else
            p_und(vol_ind,j)=p_und(vol_ind,j)+1;
        end        
    end
end

end

for i=1:n_step
    for j=1:a3

        summ=p_lysg(i,j)+p_lyss(i,j)+p_und(i,j);
        
        p_lysg(i,j)=p_lysg(i,j)/summ;
        p_lyss(i,j)=p_lyss(i,j)/summ;
        p_und(i,j)=p_und(i,j)/summ;
        
    end
end
plot(mv(:,1),p_lysg(:,1),'x',mv(:,2),p_lysg(:,2),'o',mv(:,3),p_lysg(:,3),'*',mv(:,4),p_lysg(:,4),'s',mv(:,5),p_lysg(:,5),'^','MarkerSize',15,'LineWidth',4);
set(gca,'FontSize',30,'LineWidth',2);
xlabel('$\mathcal M / V$','Interpreter','latex');
ylabel('$P_lysg$','Interpreter','latex');
legend('$\mathcal{M}=1$','$\mathcal{M}=2$','$\mathcal{M}=3$','$\mathcal{M}=4$','$\mathcal{M}=5$','Location','Southeast');
h=legend;
set(h,'Interpreter','Latex');
legend('boxoff');
title('Transiently divergent');
ylim([0 1]);
xlim([0.5 5]);
annotation('textbox','String','B','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.01 0.89 0.1 0.1]);



figure;

clear;

n_file=5;

n_step=10;

p_lysg=zeros(n_step,5);
p_lyss=zeros(n_step,5);
p_und=zeros(n_step,5);

for iii=1:n_file

loadfile=strcat('mv_det',num2str(iii));      
% asymptotically divergent

load(loadfile);

th_ci=100;
th_q=100;
% thresholds

[a1 a2 a3]=size(ciexit);

mv=zeros(n_step,5);

for i=1:a3
    mv(:,i)=i.*(0.5+1.5/n_step/2:1.5/n_step:2);
end

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
        
    end
end

end

for i=1:n_step
    for j=1:a3

        summ=p_lysg(i,j)+p_lyss(i,j)+p_und(i,j);
        
        p_lysg(i,j)=p_lysg(i,j)/summ;
        p_lyss(i,j)=p_lyss(i,j)/summ;
        p_und(i,j)=p_und(i,j)/summ;
        
    end
end

plot(mv(:,1),p_lysg(:,1),'x',mv(:,2),p_lysg(:,2),'o',mv(:,3),p_lysg(:,3),'*',mv(:,4),p_lysg(:,4),'s',mv(:,5),p_lysg(:,5),'^','MarkerSize',15,'LineWidth',4);
set(gca,'FontSize',30,'LineWidth',2);
xlabel('$\mathcal M / V$','Interpreter','latex');
ylabel('$P_{lysg}$','Interpreter','latex');
legend('$\mathcal{M}=1$','$\mathcal{M}=2$','$\mathcal{M}=3$','$\mathcal{M}=4$','$\mathcal{M}=5$','Location','Southeast');
h=legend;
set(h,'Interpreter','Latex');
legend('boxoff');
title('Asymptotically divergent');
ylim([0 1]);
xlim([0.5 5]);
annotation('textbox','String','A','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.01 0.89 0.1 0.1]);

