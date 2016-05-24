% calculate probability of states which can transcribe mRNA
%
% free energy is based on Arkin, Ross, and McAdams 98
%
% [cI_pe]=pre(z)
%
% z: cII dimer concentration (nM)
%
function [cI_pre]=pre(z)

    R=1.986*10^-3;  % kcal K^-1 mol^-1
    T=298;  % K

    G=-1.*[0 9.9 9.7 21.5];

    n_cII=[0 0 1 1];
    n_RNAP=[0 1 0 1];
    
    C=10^-9;
    RNAP=C*30;

    f=zeros(length(G),1);
    
    for k=1:length(G)
    
        f(k)=exp(-G(k)/R/T)*(C*z)^n_cII(k)*RNAP^n_RNAP(k);
        % C*z -> M

    end

    f_sum=sum(f(:));

    f(:)=f(:)./f_sum;
    
    cI_pre=f(4);
    % only 4th state can lead to cI transcription
    %
    % state 2 also allows transcription, but the transcription initiation
    % rate is more than 100 fold smaller.
end
