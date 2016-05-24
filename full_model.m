% stochastic version of cI, cro, cII, Q and antiQ model
%
% single event update by Gillespie algorithm
%
% This function uses 'pr_prm' & 'pre'
%
% [ci_exit q_exit t_dyn c_dyn]=stoc_full_model(para)
%
% input: para (array 1X24)
%
% ci_exit(i) = exit time for ci threshold i(nM)
% q_exit(j) = exit time for q threhdhold j(nM)
%
% t_dyn(:) = sampled 10 times / min
%
% c_dyn(:, 1-9) = number of mRNA and protein sampled 10 times /min
%
% c_dyn(1)=cI mRNA
% c_dyn(2)=cro mRNA
% c_dyn(3)=cII mRNA
% c_dyn(4)=Q mRNA
% c_dyn(5)=antiQ mRNA
% c_dyn(6)=cI
% c_dyn(7)=cro
% c_dyn(8)=cII
% c_dyn(9)=Q
%
%
function [ci_exit q_exit t_dyn c_dyn marker]=full_model(para)

    % redefine all the input parameters
    t_max=para(22);
    
    ci_max=para(23);
    q_max=para(24);
    
    x=zeros(9,1);
    c=zeros(9,1);
    
    N=para(1);
    
    alpha_x=para(2);
    beta_x=para(3);
    delta_x=para(4);
    
    alpha_y=para(5);
    alpha_z=para(6);
    
    alpha_q=para(7);
    
    delta_aq=para(8);
    % transcription rate
    
    gamma_x=para(9);
    gamma_y=para(10);
    gamma_z=para(11);
    gamma_q=para(12);
    % protein degradation
    
    gamma_m=para(13);
    % mRNA degradation
    
    c_d_x=para(14);
    c_d_y=para(15);
    c_d_z=para(16);
    % dimerization
    
    c_t_z=para(17);
    
    c_p_aq=para(18);
    
    sigma=para(19);
    
    zeta=para(20);
    
    vol=para(21);
    
    reaciton=zeros(21,1);
    % there are total of 21 possible reactions
    
    C=10^-9;   
    % reference concentration
    % nM
    
    marker1=0;
    marker2=0;
    marker3=0;
    t=0;
    
    for i=1:9
        c_series(1,i)=0;
        t_series(1)=0;
    end
    k=1;
    
    while (marker1==0 && marker2==0 && marker3==0)
        
        cx1=[-1+sqrt(1+8*c_d_x*c(6))]/4/c_d_x;
        cy1=[-1+sqrt(1+8*c_d_y*c(7))]/4/c_d_y;
        % equilibrium monomer concentration
        
        cx2=cx1^2*c_d_x;
        cy2=cy1^2*c_d_y;
        
        cz1=0;
        
        while cz1+2*c_d_z*cz1^2+4*c_t_z*c_d_z^2*cz1^4<c(8)
            cz1=cz1+0.001;
        end
        % determine equilibrium cII monomer concentration
        
        cz4=c_t_z*c_d_z^2*(cz1-0.001)^4;
    
        [cI_bas cI_act cro]=pr_prm(cx2,cy2);
        [cI_pre]=pre(cz4);

        reaction(1)=alpha_x*binornd(N,cI_bas);%N*cI_bas*alpha_x;%/vol;
        % cI basal transcription
        reaction(2)=beta_x*binornd(N,cI_act);%*beta_x;%/vol;
        % cI activated transcription
        reaction(3)=delta_x*binornd(N,cI_pre);%*delta_x;%/vol;
        % cI transcription by Pre

        reaction(4)=gamma_m*c(1)*vol;
        % cI mRNA degradation

        reaction(5)=binornd(N,cro)*alpha_y;%/vol;
        % cro transcription

        reaction(6)=gamma_m*c(2)*vol;
        % cro mRNA degradation

        reaction(7)=binornd(N,cro)*alpha_z;%/vol;
        % cII trancription

        reaction(8)=gamma_m*c(3)*vol;
        % cII mRNA degradation

        reaction(9)=binornd(N,cro)*alpha_q;%/vol;
        % Q transcription
       
        reaction(10)=gamma_m*c(4)*vol;
        % Q mRNA degradation
        
        reaction(11)=delta_aq*binornd(N,c_p_aq*cz4/(1+c_p_aq*cz4));
        % antiQ transcription
 
        reaction(12)=gamma_m*c(5)*vol;
        
        reaction(13)=zeta*c(4)*c(5)*vol;
        
        reaction(14)=sigma*c(1)*vol;
        reaction(15)=sigma*c(2)*vol;
        reaction(16)=sigma*c(3)*vol;
        reaction(17)=sigma*c(4)*vol;
        % translation
        % assume # ribosmes linearly scale with vol
        %
        % if translation is a two step process,
        % reaction rate=(# ribosomes is fixed)X(mRNA concentration)
        %
  
        reaction(18)=gamma_x*c(6)*vol;
        reaction(19)=gamma_y*c(7)*vol;
        reaction(20)=gamma_z*c(8)*vol;
        reaction(21)=gamma_q*c(9)*vol;
        % protein degradation
        
        % Gillespie alrogorithm
        tot_reaction=sum(reaction);

        % calculate the time of next event
        % exponential distribution
        time_update=-1/tot_reaction*log(rand);
        t=t+time_update;
        
        for i=1:length(reaction)
            reaction(i)=reaction(i)/tot_reaction;
        end

        p=rand;

        for i=1:length(reaction)-1
            cum_react(i)=sum(reaction(1:i));
        end


        if p<=cum_react(1)
            x(1)=x(1)+1;
        elseif p<=cum_react(2)
            x(1)=x(1)+1;
        elseif p<=cum_react(3)
            x(1)=x(1)+1;
        elseif p<=cum_react(4)
            x(1)=x(1)-1;
        elseif p<=cum_react(5)
            x(2)=x(2)+1;
        elseif p<=cum_react(6)
            x(2)=x(2)-1;
        elseif p<=cum_react(7)
            x(3)=x(3)+1;
        elseif p<=cum_react(8)
            x(3)=x(3)-1;
        elseif p<=cum_react(9)
            x(4)=x(4)+1;
        elseif p<=cum_react(10)
            x(4)=x(4)-1;
        elseif p<=cum_react(11)
            x(5)=x(5)+1;
        elseif p<=cum_react(12)
            x(5)=x(5)-1;            
        elseif p<=cum_react(13)    
            x(5)=x(5)-1;
            x(4)=x(4)-1;
        elseif p<=cum_react(14)
            x(6)=x(6)+1;
        elseif p<=cum_react(15)
            x(7)=x(7)+1;
        elseif p<=cum_react(16)
            x(8)=x(8)+1;
        elseif p<=cum_react(17)
            x(9)=x(9)+1;
        elseif p<=cum_react(18)
            x(6)=x(6)-1;
        elseif p<=cum_react(19)
            x(7)=x(7)-1;
        elseif p<=cum_react(20)
            x(8)=x(8)-1;
        else 
            x(9)=x(9)-1;
        end

        for i=1:9
            c(i)=x(i)/vol;
            c_series(k+1,i)=c(i);
        end

        t_series(k+1)=t;
        
        k=k+1;
        
        if t>t_max
            marker1=1;
        end
        if c(6)>ci_max
            marker2=1;
        end
        if c(9)>q_max
            marker3=1;
        end
        
    end

    ci_exit=t_max.*ones(round(ci_max),1);
    q_exit=t_max.*ones(round(q_max),1);
    
    i=2;

    t_dyn=0:0.25:t_max;
    c_dyn=zeros(length(t_dyn),9);

    % calculate exit time and truncated time series
    while i<length(t_series)+1

        ind_ci_exit=floor(c_series(i,6));
        ind_q_exit=floor(c_series(i,9));

        if (ind_ci_exit>0 && ind_ci_exit<ci_max+1)
            if ci_exit(ind_ci_exit)>t_series(i)
                for j=1:ind_ci_exit
                    ci_exit(j)=min(t_series(i),ci_exit(j));
                end
            end
        end
        
        if (ind_q_exit>0 && ind_q_exit<q_max+1)
            if q_exit(ind_q_exit)>t_series(i)
                for j=1:ind_q_exit
                    q_exit(j)=min(t_series(i),q_exit(j));
                end
            end
        end
        
        c_dyn(ceil(t_series(i)*4),:)=c_series(i,:);

        i=i+1;
    end
    
    %if marker1==1
    %    c_dyn(end,:)=[];
    %end
  
    for i=2:length(t_dyn)
        if sum(c_dyn(i,:))==0
            c_dyn(i,:)=c_dyn(i-1,:);
        end
    end
    
    marker=[marker1 marker2 marker3];

end
