% figure 6
%
% calculate M/V collapse based on dosage compsneation
%
%
clear;

comp_coeff=0.5;
% coefficient for partial gene dosage compensation

cith=100;
qth=120;
% threshold

i_max=10;
i_min=1;

sty{1}='bo';
sty{2}='r+';
sty{3}='gx';
sty{4}='ms';
sty{5}='cd';

st{1}='b.';
st{2}='r.';
st{3}='g.';
st{4}='m.';
st{5}='c.';

frac_lysg=zeros(6,5);
m_v=zeros(6,5,2);
% third index
%   1: compensated
%   2: m/v

erx=zeros(6,5,2);
% 1: uncompensated
% 2: compensated

ery=zeros(6,5);

tot_num=zeros(6);
% number with given V

num_lysg=zeros(6,5);
% 1st index=1/volume
% 2nd index=moi

for i=i_min:i_max;
    load_file=strcat('stoc_data\asymp_',num2str(comp_coeff),'_',num2str(i),'.mat');
    %load_file=strcat('stoc_data\trans_',num2str(comp_coeff),'_',num2str(i),'.mat');
    load(load_file);

    % bin volume into 6 bins
    for j=1:1000
        if 1/vol(j)<.75
            tot_num(1,1)=tot_num(1,1)+1;
            for k=1:5
                if ciexit(cith,j,k)<qexit(qth,j,k)
                    num_lysg(1,k)=num_lysg(1,k)+1;
                end
            end
        elseif 1/vol(j)<1
            tot_num(2,1)=tot_num(2,1)+1;
            for k=1:5
                if ciexit(cith,j,k)<qexit(qth,j,k)
                    num_lysg(2,k)=num_lysg(2,k)+1;
                end
            end
        elseif 1/vol(j)<1.25
            tot_num(3,1)=tot_num(3,1)+1;
            for k=1:5
                if ciexit(cith,j,k)<qexit(qth,j,k)
                    num_lysg(3,k)=num_lysg(3,k)+1;
                end
            end
        elseif 1/vol(j)<1.5
            tot_num(4,1)=tot_num(4,1)+1;
            for k=1:5
                if ciexit(cith,j,k)<qexit(qth,j,k)
                    num_lysg(4,k)=num_lysg(4,k)+1;
                end
            end
        elseif 1/vol(j)<1.75
            tot_num(5,1)=tot_num(5,1)+1;
            for k=1:5
                if ciexit(cith,j,k)<qexit(qth,j,k)
                    num_lysg(5,k)=num_lysg(5,k)+1;
                end
            end
        else
            tot_num(6,1)=tot_num(6,1)+1;
            for k=1:5
                if ciexit(cith,j,k)<qexit(qth,j,k)
                    num_lysg(6,k)=num_lysg(6,k)+1;
                end
            end
        end
    end 
end

inv_v=0.625:0.25:2;%((2+1.75)/2):-0.25:0.5;

for i=1:6
    for j=1:5
        
        frac_lysg(i,j)=num_lysg(i,j)/tot_num(i);
        
        ery(i,j)=sqrt(frac_lysg(i,j))/sqrt(tot_num(i));
        
        m_v(i,j,1)=j^comp_coeff*inv_v(i);%/inv_v(i);
        m_v(i,j,2)=j*inv_v(i);
        
        erx(i,j,1)=j^comp_coeff*0.125/2;
        erx(i,j,2)=j*0.125/2;
    end
end

%tit_text=strcat('\lambda=',num2str(comp_coeff));

figure('units','normalized','position',[0.05 .1 .5 .7]);

hold on;
for i=1:5
    %errorbar(m_v(:,i,1),frac_lysg(:,i),ery(:,i),st{i});
    %herrorbar(m_v(:,i,1),frac_lysg(:,i),erx(:,i,1),st{i});
    plot(m_v(:,i,1),frac_lysg(:,i),sty{i},'MarkerSize',12,'LineWidth',3);
end
legend('$\mathcal{M}=1$','$\mathcal{M}=2$','$\mathcal{M}=3$','$\mathcal{M}=4$','$\mathcal{M}=5$','Location','Southeast');
for i=1:5
    %errorbar(m_v(:,i,1),frac_lysg(:,i),ery(:,i),st{i});
    %herrorbar(m_v(:,i,1),frac_lysg(:,i),erx(:,i,1),st{i});
    %%plot(m_v(:,i,1),frac_lysg(:,i),sty{i},'MarkerSize',12,'LineWidth',3);
end
h=legend;
set(h,'Interpreter','Latex');
legend('Boxoff');
ax=get(gca,'Position');
ax(1)=ax(1)+0.1; %or wathever
ax(2)=ax(2)+0.05;
ax(3)=ax(3)-0.05;
ax(4)=ax(4);
set(gca,'Position',ax); 
set(gca,'FontSize',30,'LineWidth',2);
set(gca,'Box','on');
ylim([0 1]);
xlim([0 4]);
%title(tit_text);
xlabel('$\mathcal{M}^\epsilon$/V','Interpreter','Latex');
ylabel('$P_{lysg}$','Interpreter','Latex');
annotation('textbox','String','B','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.01 0.9 0.1 0.1]);


figure('units','normalized','position',[0.05 .1 .5 .7]);

hold on;
for i=1:5   
    %errorbar(m_v(:,i,2),frac_lysg(:,i),ery(:,i),st{i});
    %herrorbar(m_v(:,i,2),frac_lysg(:,i),erx(:,i,2),st{i});
    plot(m_v(:,i,2),frac_lysg(:,i),sty{i},'MarkerSize',12,'LineWidth',3);
end
legend('$\mathcal{M}=1$','$\mathcal{M}=2$','$\mathcal{M}=3$','$\mathcal{M}=4$','$\mathcal{M}=5$','Location','Southeast');
for i=1:5   
    %errorbar(m_v(:,i,2),frac_lysg(:,i),ery(:,i),st{i});
    %herrorbar(m_v(:,i,2),frac_lysg(:,i),erx(:,i,2),st{i});
    %%plot(m_v(:,i,2),frac_lysg(:,i),sty{i},'MarkerSize',12,'LineWidth',3);
end
h=legend;
set(h,'Interpreter','Latex');
legend('Boxoff');
ax=get(gca,'Position');
ax(1)=ax(1)+0.1; %or wathever
ax(2)=ax(2)+0.05;
ax(3)=ax(3)-0.05;
ax(4)=ax(4);
set(gca,'Position',ax); 
set(gca,'FontSize',30,'LineWidth',2);
ylim([0 1]);
%title(tit_text);
xlabel('$\mathcal{M}$/V','Interpreter','latex');
ylabel('$P_{lysg}$','Interpreter','Latex');
set(gca,'Box','on');
annotation('textbox','String','A','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.01 0.9 0.1 0.1]);

