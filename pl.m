% calculate probability of states which can transcribe mRNA
%
% free energy is based on Arkin, Ross, and McAdams 98
%
% [n_pl]=pl(x,y)
%
% x: cI dimer concentration (nM)
% y: cro dimer concentration (nM)
%
function [n_pl]=pl(x,y)

    R=1.986*10^-3;  % kcal K^-1 mol^-1
    T=298;  % K

    G=-1.*[0 10.9 12.1 11.7 10.1 12.5 22.9 20.9 22.8 23.7];

    n_cI=[0 0 0 1 1 0 0 1 1 2];
    n_cro=[0 1 1 0 0 0 2 1 1 0];
    n_RNAP=[0 0 0 0 0 1 0 0 0 0];
    
    C=10^-9;
    RNAP=C*30;

    f=zeros(length(G),1);
    
    for k=1:length(G)
    
        f(k)=exp(-G(k)/R/T)*(C*x)^n_cI(k)*(C*y)^n_cro(k)*RNAP^n_RNAP(k);
        % C*z -> M

    end

    f_sum=sum(f(:));

    f(:)=f(:)./f_sum;
    
    n_pl=f(6);
    % only 4th state can lead to cI transcription
    %
    % state 2 also allows transcription, but the transcription initiation
    % rate is more than 100 fold smaller.
end
