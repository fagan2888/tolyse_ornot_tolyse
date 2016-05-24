% figure S2
%
% draw heatmap of mean decision time for transicently driven case
%
% moi=1
%
% exclude undecided case
%
clear;

moi=1;

load para_sto_series3;
% transiently div
% moi=1

[a1 a2 a3]=size(ciexit);

t_dec=zeros(300,200);

f_lysg=zeros(300,200);

for i=1:300
    for j=1:200
        
        t_sum=0;
        t_num=0;
        
        for k=1:a2
            %if (ciexit(i,k,1)<1000 || qexit(j,k,1)<1000)
                t_sum=t_sum+min(ciexit(i,k,moi),qexit(j,k,moi));
                t_num=t_num+1;
            %end
            
            if ciexit(i,k,moi)<qexit(j,k,moi)
                f_lysg(i,j)=f_lysg(i,j)+1/a2;
            end
            % fraction of lysogeny
        end
        
        if t_num>0
            t_dec(i,j)=t_sum/t_num;
        end
        % average decision time
        
    end
    
    fprintf('%d / 500\n',i);
    
end

xgrid=1:300;
ygrid=1:200;

[xx yy]=meshgrid(xgrid,ygrid);

hold on;
contourf(xx,yy,transpose(log10(t_dec)),100,'LineStyle','none');
plot(100,100,'k.','Markersize',60);
set(gca,'FontSize',30);
xlabel('cI threshold (nM)');
ylabel('Q threshold (nM)');
colorbar;
title('Average decision time (min)');
colorbar('YTick',[0.699, 1, 1.30, 1.69, 2 ,2.30, 2.69],'YTickLabel',{'5', '10', '20', '50', '100', '200', '500'});
annotation('textbox','String','B','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.01 0.89 0.1 0.1]);

f_lysg(1,end)=0;
f_lysg(end,1)=1;

figure;
hold on;
contourf(xx,yy,transpose(f_lysg),100,'LineStyle','none');
plot(100,100,'k.','Markersize',60);
hold off;
set(gca,'FontSize',30);
xlabel('cI threshold (nM)');
ylabel('Q threshold (nM)');
colorbar;
title('Fraction of Lysogeny');
annotation('textbox','String','A','FontSize',40,'LineStyle','none','FontWeight','bold','Position',[0.01 0.89 0.1 0.1]);
