% Dynamics of cI, cro, cII, Q and int at the protein and mRNA level
%
% a version with N and cII tetramerization
%
% no dynamics for N
% degradation occurs in monomer, dimer, tetramer
% total protein concentration
%
% dxdt=model1(t, x, para)
% 
% x(1)=cI mRNA
% x(2)=cro mRNA
% x(3)=cII mRNA
% x(4)=aQ mRNA
% x(5)=Q mRNA
% x(6)=cI
% x(7)=cro
% x(8)=cII
% x(9)=Q
%
function dxdt=model_final(t,x,para)
    
    dxdt=zeros(9,1);
    
    N=para(1);
    % copy number
    
    alpha_x=para(2);
    beta_x=para(3);
    delta_x=para(4);
    alpha_y=para(5);
    alpha_z=para(6);
    alpha_q=para(7);
    delta_aq=para(8);
    % transcription rate
    % alpha: basal
    % beta: cI activation
    % delta: cII activation
    
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
    % dimerization rate
    
    c_t_z=para(17);
    % tetramerization rate
    
    c_p_aq=para(18);
    % binding constant for cII tetramer
    
    sigma=para(19);
    % translation
    
    zeta=para(20);
    % binding between Q and antiQ mRNA
    
    C=10^-9;

    
    % monomer + dimer = total protein
    % ignore proteins bound to operators
    %
    % x1+2*c_d_x*x1^2=X_tot
    %
    % x1=[-1+sqrt(1+8 c_d_x X_tot)]/4 c_d_x
    
    x1=[-1+sqrt(1+8*c_d_x*x(6))]/4/c_d_x;
    y1=[-1+sqrt(1+8*c_d_y*x(7))]/4/c_d_y;
 
    x2=c_d_x*x1^2;
    y2=c_d_y*y1^2;
   
    z1=0;
    
    z_step=0.001;
    
    while 4*z1^4*c_d_z^2*c_t_z+2*z1^2*c_d_z+z1<x(8)
        z1=z1+z_step;
    end
    
    z2=c_d_z*(z1-z_step)^2;
    z4=c_t_z*z2^2;
    z1=x(9)-2*z2-4*z4;
    
    [cI_bas cI_act cro]=pr_prm(x2,y2);
    [cI_pre]=pre(z4);
    [n_pl]=pl(x2,y2);
    
    dxdt(1)=N*cI_bas*alpha_x+N*cI_act*beta_x+N*cI_pre*delta_x-gamma_m*x(1);
    % cI mRNA
    dxdt(2)=N*cro*alpha_y-gamma_m*x(2);
    % cro mRNA
    dxdt(3)=N*cro*alpha_z-gamma_m*x(3);
    % cII mRNA
    dxdt(4)=N*cro*alpha_q-gamma_m*x(4)-zeta*x(4)*x(5);
    % Q mRNA
    dxdt(5)=N*delta_aq*c_p_aq*z4/(1+c_p_aq*z4)-gamma_m*x(5)-zeta*x(4)*x(5);
    % aQ mRNA 
    %
    % antisense mRNA also decrease with degradation
    %
    dxdt(6)=sigma*x(1)-gamma_x*x(6);%*x1;
    % cI
    dxdt(7)=sigma*x(2)-gamma_y*x(7);%y1;
    % cro
    dxdt(8)=sigma*x(3)-gamma_z*x(8);%z1;
    % cII
    dxdt(9)=sigma*x(4)-gamma_q*x(9);
    % Q
end