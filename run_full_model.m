% main_stoc_full_model
%
%

clear;

it_max=100;
it_rep=1;

tic;

ci_max=110;%300;%500;
q_max=110;%150;%200;

t_max=1000;

% alpha_x=0.06;
% % shea 85
% beta_x=0.66;
% % shea 85
% delta_x=0.9;
% % arkin 98
% alpha_y=0.84;
% alpha_z=0.8;
% alpha_q=0.75;
% delta_aq=2;
%  
% gamma_x=0.01;
% gamma_y=0.06;%4;
% gamma_z=0.10;
% gamma_q=0.01;%0.01;
% 
% gamma_m=0.1;
% 
% c_d_x=1/20;
% % Darling 2000
% c_d_y=1/(1/5.8);%1/100;
% % Jana et al. JMB 97
% c_d_z=1/20;
% 
% c_t_z=1/20;%1/10;
% 
% c_p_aq=1/5;
% 
% sigma=0.5;
% 
% zeta=0.1;
% 
% transcp=[alpha_x beta_x delta_x alpha_y alpha_z alpha_q delta_aq];
% gamma=[gamma_x gamma_y gamma_z gamma_q gamma_m];
% dim=[c_d_x c_d_y c_d_z];
% tet=c_t_z;
% 
% para_model=[transcp gamma dim tet c_p_aq sigma zeta];

load('data_para_series2');
para_model=transpose(para_series(:,215));
% stochastic
% no seperation of steady states

%load data_para_series3;
%para_model=transpose(para_series(:,99));
% deterministic
% separation of steady states

for it_ind=1:it_rep

    ciexit=zeros(ci_max,it_max,5);
    qexit=zeros(q_max,it_max,5);

    ci_dyn=zeros(4*t_max+1,it_max,5);
    q_dyn=zeros(4*t_max+1,it_max,5);

    vol=zeros(it_max,1);

    for i=1:it_max

        %inv_vol=0.5+1.5*rand;
        %vol(i)=1/inv_vol;
        %vol(i)=0.5+1.5*rand;

        N=1;
        vol=0.5;
        para1=[N para_model vol t_max ci_max q_max];
        N=2;
        vol=1;
        para2=[N para_model vol t_max ci_max q_max];
        N=3;
        vol=1.5;
        para3=[N para_model vol t_max ci_max q_max];
        N=4;
        vol=2;
        para4=[N para_model vol t_max ci_max q_max];
        N=5;
        vol=2.5;
        para5=[N para_model vol t_max ci_max q_max];

        [ci_exit1 q_exit1 t_dyn1 c_dyn1 m1]=full_model(para1);
        [ci_exit2 q_exit2 t_dyn2 c_dyn2 m2]=full_model(para2);
        [ci_exit3 q_exit3 t_dyn3 c_dyn3 m3]=full_model(para3);
        [ci_exit4 q_exit4 t_dyn4 c_dyn4 m4]=full_model(para4);
        [ci_exit5 q_exit5 t_dyn5 c_dyn5 m5]=full_model(para5);

        for j=1:ci_max
            ciexit(j,i,1)=ci_exit1(j);
            ciexit(j,i,2)=ci_exit2(j);
            ciexit(j,i,3)=ci_exit3(j);
            ciexit(j,i,4)=ci_exit4(j);
            ciexit(j,i,5)=ci_exit5(j);
        end
        for j=1:q_max
            qexit(j,i,1)=q_exit1(j);
            qexit(j,i,2)=q_exit2(j);
            qexit(j,i,3)=q_exit3(j);
            qexit(j,i,4)=q_exit4(j);
            qexit(j,i,5)=q_exit5(j);
        end

        for j=1:4*t_max+1
            q_dyn(j,i,1)=c_dyn1(j,9);
            q_dyn(j,i,2)=c_dyn2(j,9);
            q_dyn(j,i,3)=c_dyn3(j,9);
            q_dyn(j,i,4)=c_dyn4(j,9);
            q_dyn(j,i,5)=c_dyn5(j,9);

            ci_dyn(j,i,1)=c_dyn1(j,6);
            ci_dyn(j,i,2)=c_dyn2(j,6);
            ci_dyn(j,i,3)=c_dyn3(j,6);
            ci_dyn(j,i,4)=c_dyn4(j,6);
            ci_dyn(j,i,5)=c_dyn5(j,6);
        end

        fprintf('%d / ',i);
        fprintf('  %d   \n',it_max)

    end
    
    save('test_data1');

end