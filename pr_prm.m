% calculate probability of states which can transcribe mRNA
%
% free energy is based on Shea Ackers 1985
%
% [cI_basal cI_act cro]=pr_prm(x,y)
%
% x: cI dimer concentration (nM)
% y: cro dimer concentration (nM)
%
function [cI_basal cI_act cro]=pr_prm(x,y)

    R=1.986*10^-3;  % kcal K^-1 mol^-1
    T=298;  % K

    G=-1.*[0 11.7 10.1 10.1 10.8 10.8 12.1 11.5 12.5 23.7 21.8 22.2 21.6 22.9 22.9 24.0 22.5 20.9 20.9 23.8 20.9 22.2 22.6 21.6 23.2 24.6 22.3 22.3 33.8 33.7 35.8 32.6 33.0 31.7 33.0 34.6 35.2 33.1 34.0 32.4]; % kcal
    % Gibbs free energy

    n_cI=[0 1 1 1 0 0 0 0 0 2 2 2 0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 0 3 0 2 2 2 1 1 1 2 0 1 1]; 
    % # of cI dimers bound

    n_cro=[0 0 0 0 1 1 1 0 0 0 0 0 2 2 2 0 1 1 1 1 1 1 0 0 0 1 1 1 0 3 1 1 1 2 2 2 0 2 1 1];
    % # of cro dimers bound

    n_RNAP=[0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 2 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1];
    % # of RNAP bound

    C=10^-9;
    z=C*30;
    
    index_cI_basal=[8 16 25 27 28 38 39];
    index_cI_act=[24 37 40];
    index_cro=[9 16 23 26];
    % configuration of operators with possible transcription
    
    f=zeros(length(G),1);
    
    for k=1:length(G)
    
        f(k)=exp(-G(k)/R/T)*(C*x)^n_cI(k)*(C*y)^n_cro(k)*z^n_RNAP(k);
        % C*x: nM->M
    end

    f_sum=sum(f(:));

    f(:)=f(:)./f_sum;

    cI_basal=0;
    cI_act=0;
    cro=0;
    
    for k=1:length(index_cI_basal)

        cI_basal=cI_basal+f(index_cI_basal(k));

    end    
    for k=1:length(index_cI_act)

        cI_act=cI_act+f(index_cI_act(k));

    end
    for k=1:length(index_cro)

        cro=cro+f(index_cro(k));

    end  
end
