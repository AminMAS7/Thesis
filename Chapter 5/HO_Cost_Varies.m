
clc
clear all
close all



Mainseed = 7;
rng(Mainseed)
%%%%%%%%%%%%%%%%%%%%%% THz
load("Molecular_absorption.mat","BW","K_abs");

c = physconst('LightSpeed');
k_B = 1.381e-23;  % Boltzmann constant
T0 = 290; % noise temperature (Kelvin)

General_Struct = {};
General_Struct.c = c;
% THz:
General_Struct.f_THz = 0.4e12;
General_Struct.k_a = K_abs(BW == General_Struct.f_THz);
General_Struct.THz_BW = 0.8e9;
General_Struct.Ther_Noise_THz_Orig = k_B*T0*General_Struct.THz_BW;
Pmax_THz_Tot_dBm = 25; % [dBm]
General_Struct.Pmax_THz = db2pow(Pmax_THz_Tot_dBm-30);
% RF:
General_Struct.f_RF = 8e9;
General_Struct.RF_BW = 100e6;
General_Struct.Ther_Noise_RF_Orig = k_B*T0*General_Struct.RF_BW;
Pmax_RF_Tot_dBm = 40; % [dBm]
General_Struct.Pmax_RF = db2pow(Pmax_RF_Tot_dBm-30);
% MBN:
General_Struct.RF_BW_Norm = General_Struct.RF_BW/General_Struct.THz_BW;


Traj_data = {};

v = 50; % Users velocity

Traj_data.v = v;


Traj_data.Trajec_Time = 0.2;
Traj_data.tau_eps = 0.1;
N = ceil(Traj_data.Trajec_Time/Traj_data.tau_eps);
Traj_data.N = N;


Traj_data.N_Hori = ceil( 0.1 / Traj_data.tau_eps);
General_Struct.T_NUA = 20;


tau_eps_GT = 1e-2;
N_GT = ceil((Traj_data.tau_eps)/tau_eps_GT);

RF_Coh_tau = 1e-2;

Traj_data.N_RF_Coh = RF_Coh_tau/tau_eps_GT;

Traj_data.tau_eps_GT = tau_eps_GT;
Traj_data.N_GT = N_GT;


H_min = 10;
% RF Antenna Gain
General_Struct.G_Tx_RF = db2pow(10); % [dB]
General_Struct.G_Rx_RF = db2pow(8); % [dB]

%%%%%%%%%%% THz Antenna Gain
General_Struct.G_Tx_THz = db2pow(15); % [dB]
General_Struct.G_Rx_THz = db2pow(8); % [dB]

% Main Parameters loading:
K = 15; % Number of users
B_T = 5; % Number of TBSs
B_R = 3; % Number of RBSs

M = 504;
M = K * ceil(M/K);
NumAnt_per_SubArr_THz = M/K;


M_hat = 80;
M_hat = K * ceil(M_hat/K);
NumAnt_per_SubArr_RF = M_hat/K;
General_Struct.RF_PL_expon = 2.5;
General_Struct.RicianFac = 100;

General_Struct.THz_Hyb = 1;
General_Struct.NFC_flg_THz = 1;
General_Struct.Part_Conn_THz = 0;
General_Struct.Full_Conn_THz = (~General_Struct.Part_Conn_THz);


General_Struct.RF_Hyb = 1;
General_Struct.NFC_flg_RF = 1;
General_Struct.Part_Conn_RF = 0;
General_Struct.Full_Conn_RF = (~General_Struct.Part_Conn_RF);


General_Struct.err_flg = 0;


% Positions loading:
Channel_info = {};
Channel_info.L_x_max = 350;
Channel_info.W_y_max = 250;
Channel_info.W_min_THz = 35;

%%% Parameters Loading:
eta_BL = 1e-3;

Clus_Siz_THz = 2;
Clus_Siz_RF = 2;
Clus_Siz_MBN = 4;

Rth = 0.7e9;


n = 1;
eps_conv = 1e-3;



% Channel Loading:
seed = Mainseed;


THz_Hyb = General_Struct.THz_Hyb;
RF_Hyb = General_Struct.RF_Hyb;
Pmax_THz = General_Struct.Pmax_THz;
Pmax_RF = General_Struct.Pmax_RF;
THz_BW = General_Struct.THz_BW;
RF_BW_Norm = General_Struct.RF_BW_Norm;
RF_BW = General_Struct.RF_BW;

Rth_vec = Rth*ones(1,K);
General_Struct.S_HO = 1;


Ch_err = 0;

General_Struct.Mobil_flg = 0;

General_Struct.MRT_Init_We = 1;

HO_Pen_1 = 1;
HO_Pen_2 = 1;

[THz_BS_Loc,RF_BS_Loc,All_BS_Loc] = BS_Pos(B_T,B_R,Channel_info);

Users_in_X = rand(K,1);

L_x_max_Users = 135;
Users_Init_x = (L_x_max_Users)*Users_in_X;
[Users_x_mat,Users_Init_y] = Users_pos(K,THz_BS_Loc,N,Traj_data.tau_eps,v,Users_Init_x,Channel_info);

Max_Min_Init = 1;
General_Struct.MaxMin_flg = 1;

MNT = 50*6;
Seed_Vec = 1:MNT;

eta_Vec = linspace(0,0.95,10);


for mnt = Seed_Vec
    seed = Seed_Vec(mnt);
    mnt

    rng(seed)
    Users_in_X = rand(K,1);
    Users_Init_x = (L_x_max_Users)*Users_in_X;

     [F_HO_AveCh_err_LogCo(:,mnt), ~, F_HO_AveRate_err_LogCo(:,mnt), ~, NuHO_MBN_LogCo(:,mnt), ...
    F_HO_AveCh_err_P1(:,mnt), ~, F_HO_AveRate_err_P1(:,mnt), ~, NuHO_MBN_P1(:,mnt),...
    F_HO_AveCh_err_NoCo(:,mnt), ~, F_HO_AveRate_err_NoCo(:,mnt), ~, NuHO_MBN_NoCo(:,mnt),...
    F_HO_AveCh_err_B1(:,mnt), ~,F_HO_AveRate_err_B1(:,mnt), ~, NuHO_MBN_B1(:,mnt),...
    F_HO_AveCh_err_TOn_LogCo(:,mnt), ~,F_HO_AveRate_err_TOn_LogCo(:,mnt), ~, NuHO_TOn_LogCo(:,mnt), ...
    F_HO_AveCh_err_TOn_P1(:,mnt), ~,F_HO_AveRate_err_TOn_P1(:,mnt), ~, NuHO_TOn_P1(:,mnt),...
    F_HO_AveCh_err_TOn_NoCo(:,mnt), ~,F_HO_AveRate_err_TOn_NoCo(:,mnt), ~, NuHO_TOn_NoCo(:,mnt)] = Inner_Func(Traj_data,K,B_R,M_hat,Channel_info,...
        H_min,General_Struct,B_T,M,...
    eps_conv,Pmax_THz,Rth_vec,THz_BW,...
    Clus_Siz_THz,Clus_Siz_RF,Clus_Siz_MBN,THz_Hyb,eta_BL,...
    Pmax_RF,RF_BW_Norm,RF_Hyb,RF_BW,Max_Min_Init,HO_Pen_1,HO_Pen_2,...
    Ch_err,Users_Init_x,THz_BS_Loc,RF_BS_Loc,All_BS_Loc,eta_Vec,seed);

end

%%% MBN:
F_HO_AveCh_err_LogCo_Mean = mean(F_HO_AveCh_err_LogCo,2);
F_HO_AveRate_err_LogCo_Mean = mean(F_HO_AveRate_err_LogCo,2);

NuHO_MBN_LogCo_Mean = mean(NuHO_MBN_LogCo,2);


F_HO_AveCh_err_P1_Mean = mean(F_HO_AveCh_err_P1,2);
F_HO_AveRate_err_P1_Mean = mean(F_HO_AveRate_err_P1,2);

NuHO_MBN_P1_Mean = mean(NuHO_MBN_P1,2);


F_HO_AveCh_err_NoCo_Mean = mean(F_HO_AveCh_err_NoCo,2);
F_HO_AveRate_err_NoCo_Mean = mean(F_HO_AveRate_err_NoCo,2);

NuHO_MBN_NoCo_Mean = mean(NuHO_MBN_NoCo,2);


F_HO_AveCh_err_B1_Mean = mean(F_HO_AveCh_err_B1,2);
F_HO_AveRate_err_B1_Mean = mean(F_HO_AveRate_err_B1,2);

NuHO_MBN_B1_Mean = mean(NuHO_MBN_B1,2);


%%%% THz Only:
F_HO_AveCh_err_TOn_LogCo_Mean = mean(F_HO_AveCh_err_TOn_LogCo,2);
F_HO_AveRate_err_TOn_LogCo_Mean = mean(F_HO_AveRate_err_TOn_LogCo,2);

NuHO_TOn_LogCo_Mean = mean(NuHO_TOn_LogCo,2);


F_HO_AveCh_err_TOn_P1_Mean = mean(F_HO_AveCh_err_TOn_P1,2);
F_HO_AveRate_err_TOn_P1_Mean = mean(F_HO_AveRate_err_TOn_P1,2);

NuHO_TOn_P1_Mean = mean(NuHO_TOn_P1,2);


F_HO_AveCh_err_TOn_NoCo_Mean = mean(F_HO_AveCh_err_TOn_NoCo,2);
F_HO_AveRate_err_TOn_NoCo_Mean = mean(F_HO_AveRate_err_TOn_NoCo,2);

NuHO_TOn_NoCo_Mean = mean(NuHO_TOn_NoCo,2);

Color1 = [0, 0, 1]; % blue
Color2 = [0, 0.5, 0]; % green
Color3 = [1, 0, 0]; % red 
Color4 = [0.25, 0.25, 0.25]; % black
Color7 = [0, 0.75, 0.75]; % cyan
Color5 = [0.75, 0.75, 0];
Color6 = [0.75, 0, 0.75]; % purple
Color8 = [0.8500, 0.3250, 0.0980]; % dark orange


Gb_factor = (1/1e9);
figure(2)
box on
plot(eta_Vec,Gb_factor*F_HO_AveCh_err_LogCo_Mean,'-','Color',Color1,'MarkerIndices',1:2:length(eta_Vec),'LineWidth',0.8)
hold on
plot(eta_Vec,Gb_factor*F_HO_AveCh_err_P1_Mean,'-s','Color',Color2,'MarkerIndices',1:2:length(eta_Vec),'LineWidth',0.8)
hold on
plot(eta_Vec,Gb_factor*F_HO_AveCh_err_NoCo_Mean,'-','Color',Color4,'MarkerIndices',1:2:length(eta_Vec),'LineWidth',0.8)
hold on
plot(eta_Vec,Gb_factor*F_HO_AveCh_err_B1_Mean,'-','Color',Color3,'MarkerIndices',1:2:length(eta_Vec),'LineWidth',0.8)
hold on
plot(eta_Vec,Gb_factor*F_HO_AveCh_err_TOn_LogCo_Mean,'--*','Color',Color1,'MarkerIndices',1:2:length(eta_Vec),'LineWidth',0.8)
hold on
plot(eta_Vec,Gb_factor*F_HO_AveCh_err_TOn_P1_Mean,'--*','Color',Color2,'MarkerIndices',1:2:length(eta_Vec),'LineWidth',0.8)
hold on
plot(eta_Vec,Gb_factor*F_HO_AveCh_err_TOn_NoCo_Mean,'--*','Color',Color4,'MarkerIndices',1:2:length(eta_Vec),'LineWidth',0.8)

grid on

legend({'MBN-Algo1-Cost','MBN-Algo1-MO','MBN-Algo1(No Cost)','MBN-B1(No Cost)','THzOn-Algo1-Cost','THzOn-Algo1-MO','THzOn-Algo1(No Cost)'},'FontSize',9,'Interpreter','latex')
ylabel('Average HO-aware Sum Rate [Gbps]', 'FontSize', 10)
xlabel('HO Cost, \eta^T = \eta^M', 'FontSize', 10)


figure(4)
box on
plot(eta_Vec,1*NuHO_MBN_LogCo_Mean,'-','Color',Color1,'MarkerIndices',1:2:length(eta_Vec),'LineWidth',0.8)
hold on
plot(eta_Vec,1*NuHO_MBN_P1_Mean,'-s','Color',Color2,'MarkerIndices',1:2:length(eta_Vec),'LineWidth',0.8)
hold on
plot(eta_Vec,1*NuHO_MBN_NoCo_Mean,'-','Color',Color4,'MarkerIndices',1:2:length(eta_Vec),'LineWidth',0.8)
hold on
plot(eta_Vec,1*NuHO_MBN_B1_Mean,'-','Color',Color3,'MarkerIndices',1:2:length(eta_Vec),'LineWidth',0.8)
hold on
plot(eta_Vec,1*NuHO_TOn_LogCo_Mean,'--*','Color',Color1,'MarkerIndices',1:2:length(eta_Vec),'LineWidth',0.8)
hold on
plot(eta_Vec,1*NuHO_TOn_P1_Mean,'--*','Color',Color2,'MarkerIndices',1:2:length(eta_Vec),'LineWidth',0.8)
hold on
plot(eta_Vec,1*NuHO_TOn_NoCo_Mean,'--*','Color',Color4,'MarkerIndices',1:2:length(eta_Vec),'LineWidth',0.8)

grid on

legend({'MBN-Algo1-Cost','MBN-Algo1-MO','MBN-Algo1(No Cost)','MBN-B1(No Cost)','THzOn-Algo1-Cost','THzOn-Algo1-MO','THzOn-Algo1(No Cost)'},'FontSize',9,'Interpreter','latex')
ylabel('Average Number of HOs', 'FontSize', 10)
xlabel('HO Cost, \eta^T = \eta^M', 'FontSize', 10)



function [F_HO_AveCh_err_LogCo, F_HO_AveCh_GT_LogCo, F_HO_AveRate_err_LogCo, F_HO_AveRate_GT_LogCo, NuHO_MBN_LogCo, ...
    F_HO_AveCh_err_P1, F_HO_AveCh_GT_P1, F_HO_AveRate_err_P1, F_HO_AveRate_GT_P1, NuHO_MBN_P1,...
    F_HO_AveCh_err_NoCo, F_HO_AveCh_GT_NoCo, F_HO_AveRate_err_NoCo, F_HO_AveRate_GT_NoCo, NuHO_MBN_NoCo,...
    F_HO_AveCh_err_B1, F_HO_AveCh_GT_B1,F_HO_AveRate_err_B1, F_HO_AveRate_GT_B1, NuHO_MBN_B1,...
    F_HO_AveCh_err_TOn_LogCo, F_HO_AveCh_GT_TOn_LogCo,F_HO_AveRate_err_TOn_LogCo, F_HO_AveRate_GT_TOn_LogCo, NuHO_TOn_LogCo, ...
    F_HO_AveCh_err_TOn_P1, F_HO_AveCh_GT_TOn_P1,F_HO_AveRate_err_TOn_P1, F_HO_AveRate_GT_TOn_P1, NuHO_TOn_P1,...
    F_HO_AveCh_err_TOn_NoCo, F_HO_AveCh_GT_TOn_NoCo,F_HO_AveRate_err_TOn_NoCo, F_HO_AveRate_GT_TOn_NoCo, NuHO_TOn_NoCo] = Inner_Func(Traj_data,K,B_R,M_hat,Channel_info,...
        H_min,General_Struct,B_T,M,...
    eps_conv,Pmax_THz,Rth_vec,THz_BW,...
    Clus_Siz_THz,Clus_Siz_RF,Clus_Siz_MBN,THz_Hyb,eta_BL,...
    Pmax_RF,RF_BW_Norm,RF_Hyb,RF_BW,Max_Min_Init,HO_Pen_1,HO_Pen_2,...
    Ch_err,Users_Init_x,THz_BS_Loc,RF_BS_Loc,All_BS_Loc,eta_Vec,seed)


for i = 1:length(eta_Vec)
    eta_THz = eta_Vec(i);
    eta_RF = eta_Vec(i);

    Func_Type = 'WithCost';
    [F_HO_AveCh_err_LogCo(i),F_HO_AveRate_err_LogCo(i),F_HO_AveCh_GT_LogCo(i),F_HO_AveRate_GT_LogCo(i),NuHO_MBN_LogCo(i),~,~,~,~,~] = MBN_Mobile_Func(Traj_data,K,B_R,M_hat,Channel_info,...
        H_min,General_Struct,B_T,M,...
    eps_conv,Pmax_THz,Rth_vec,THz_BW,...
    Clus_Siz_THz,Clus_Siz_RF,Clus_Siz_MBN,THz_Hyb,Func_Type,eta_BL,...
    Pmax_RF,RF_BW_Norm,RF_Hyb,RF_BW,Max_Min_Init,eta_THz,eta_RF,HO_Pen_1,HO_Pen_2,Ch_err,Users_Init_x,THz_BS_Loc,RF_BS_Loc,seed);



    %%%%% THz only:

    Func_Type = 'THzOn-WithCost';
    [F_HO_AveCh_err_TOn_LogCo(i),F_HO_AveRate_err_TOn_LogCo(i),F_HO_AveCh_GT_TOn_LogCo(i),F_HO_AveRate_GT_TOn_LogCo(i),NuHO_TOn_LogCo(i),~,~,~] = THzOn_Mobile_Func(Traj_data,K,B_R,Channel_info,...
        H_min,General_Struct,B_T,M,eps_conv,Pmax_THz,Rth_vec,THz_BW,Clus_Siz_MBN,THz_Hyb,Func_Type,eta_BL,...
    Max_Min_Init,eta_THz,HO_Pen_1,HO_Pen_2,Ch_err,Users_Init_x,All_BS_Loc,seed);


    if i == 1
        Func_Type = 'Pen-1';
        [F_HO_AveCh_err_P1(i),F_HO_AveRate_err_P1(i),F_HO_AveCh_GT_P1(i),F_HO_AveRate_GT_P1(i),NumHO_MBN_P1,W_MAT_P1,U_MAT_P1,alpha_mat_P1,beta_mat_P1,~] = MBN_Mobile_Func(Traj_data,K,B_R,M_hat,Channel_info,...
            H_min,General_Struct,B_T,M,...
        eps_conv,Pmax_THz,Rth_vec,THz_BW,...
        Clus_Siz_THz,Clus_Siz_RF,Clus_Siz_MBN,THz_Hyb,Func_Type,eta_BL,...
        Pmax_RF,RF_BW_Norm,RF_Hyb,RF_BW,Max_Min_Init,eta_THz,eta_RF,HO_Pen_1,HO_Pen_2,Ch_err,Users_Init_x,THz_BS_Loc,RF_BS_Loc,seed);
        NuHO_MBN_P1(i) = NumHO_MBN_P1;
    else
        [F_HO_AveCh_err_P1(i),F_HO_AveRate_err_P1(i),F_HO_AveCh_GT_P1(i),F_HO_AveRate_GT_P1(i)] = MBN_Mobile_Func_HO(Traj_data,K,B_R,M_hat,Channel_info,...
            H_min,General_Struct,B_T,M,THz_BW,W_MAT_P1,U_MAT_P1,alpha_mat_P1,beta_mat_P1,...
            THz_Hyb,eta_BL,RF_Hyb,RF_BW,eta_THz,eta_RF,Ch_err,Users_Init_x,THz_BS_Loc,RF_BS_Loc,Rth_vec,seed);
        NuHO_MBN_P1(i) = NumHO_MBN_P1;
    end
    
   
    if i == 1
        Func_Type = 'NoCost';
        [F_HO_AveCh_err_NoCo(i),F_HO_AveRate_err_NoCo(i),F_HO_AveCh_GT_NoCo(i),F_HO_AveRate_GT_NoCo(i),NumHO_MBN_NoCo,W_MAT_NoCo,U_MAT_NoCo,alpha_mat_NoCo,beta_mat_NoCo,~] = MBN_Mobile_Func(Traj_data,K,B_R,M_hat,Channel_info,...
            H_min,General_Struct,B_T,M,...
        eps_conv,Pmax_THz,Rth_vec,THz_BW,...
        Clus_Siz_THz,Clus_Siz_RF,Clus_Siz_MBN,THz_Hyb,Func_Type,eta_BL,...
        Pmax_RF,RF_BW_Norm,RF_Hyb,RF_BW,Max_Min_Init,eta_THz,eta_RF,HO_Pen_1,HO_Pen_2,Ch_err,Users_Init_x,THz_BS_Loc,RF_BS_Loc,seed);
        NuHO_MBN_NoCo(i) = NumHO_MBN_NoCo;
        else
        [F_HO_AveCh_err_NoCo(i),F_HO_AveRate_err_NoCo(i),F_HO_AveCh_GT_NoCo(i),F_HO_AveRate_GT_NoCo(i)] = MBN_Mobile_Func_HO(Traj_data,K,B_R,M_hat,Channel_info,...
            H_min,General_Struct,B_T,M,THz_BW,W_MAT_NoCo,U_MAT_NoCo,alpha_mat_NoCo,beta_mat_NoCo,...
            THz_Hyb,eta_BL,RF_Hyb,RF_BW,eta_THz,eta_RF,Ch_err,Users_Init_x,THz_BS_Loc,RF_BS_Loc,Rth_vec,seed);
        NuHO_MBN_NoCo(i) = NumHO_MBN_NoCo;

    end

    if i == 1
        Func_Type = 'MBN-B1';
        [F_HO_AveCh_err_B1(i),F_HO_AveRate_err_B1(i),F_HO_AveCh_GT_B1(i),F_HO_AveRate_GT_B1(i),NumHO_MBN_B1,W_MAT_B1,U_MAT_B1,alpha_mat_B1,beta_mat_B1,~] = MBN_Mobile_Func(Traj_data,K,B_R,M_hat,Channel_info,...
            H_min,General_Struct,B_T,M,...
        eps_conv,Pmax_THz,Rth_vec,THz_BW,...
        Clus_Siz_THz,Clus_Siz_RF,Clus_Siz_MBN,THz_Hyb,Func_Type,eta_BL,...
        Pmax_RF,RF_BW_Norm,RF_Hyb,RF_BW,Max_Min_Init,eta_THz,eta_RF,HO_Pen_1,HO_Pen_2,Ch_err,Users_Init_x,THz_BS_Loc,RF_BS_Loc,seed);
        NuHO_MBN_B1(i) = NumHO_MBN_B1;
    else
        [F_HO_AveCh_err_B1(i),F_HO_AveRate_err_B1(i),F_HO_AveCh_GT_B1(i),F_HO_AveRate_GT_B1(i)] = MBN_Mobile_Func_HO(Traj_data,K,B_R,M_hat,Channel_info,...
            H_min,General_Struct,B_T,M,THz_BW,W_MAT_B1,U_MAT_B1,alpha_mat_B1,beta_mat_B1,...
            THz_Hyb,eta_BL,RF_Hyb,RF_BW,eta_THz,eta_RF,Ch_err,Users_Init_x,THz_BS_Loc,RF_BS_Loc,Rth_vec,seed);
        NuHO_MBN_B1(i) = NumHO_MBN_B1;

    end

    if i == 1
        Func_Type = 'THzOn-Pen-1';
        [F_HO_AveCh_err_TOn_P1(i),F_HO_AveRate_err_TOn_P1(i),F_HO_AveCh_GT_TOn_P1(i),F_HO_AveRate_GT_TOn_P1(i),NumHO_TOn_P1,W_MAT_TOn_P1,alpha_mat_TOn_P1,~] = THzOn_Mobile_Func(Traj_data,K,B_R,Channel_info,...
            H_min,General_Struct,B_T,M,eps_conv,Pmax_THz,Rth_vec,THz_BW,Clus_Siz_MBN,THz_Hyb,Func_Type,eta_BL,...
        Max_Min_Init,eta_THz,HO_Pen_1,HO_Pen_2,Ch_err,Users_Init_x,All_BS_Loc,seed);
        NuHO_TOn_P1(i) = NumHO_TOn_P1;
    else
        [F_HO_AveCh_err_TOn_P1(i),F_HO_AveRate_err_TOn_P1(i),F_HO_AveCh_GT_TOn_P1(i),F_HO_AveRate_GT_TOn_P1(i)] = THzOn_Mobile_Func_HO(Traj_data,K,B_R,Channel_info,...
            H_min,General_Struct,B_T,M,THz_BW,W_MAT_TOn_P1,alpha_mat_TOn_P1,THz_Hyb,eta_BL,...
        eta_THz,Ch_err,Users_Init_x,All_BS_Loc,Rth_vec,seed);
        NuHO_TOn_P1(i) = NumHO_TOn_P1;

    end


    if i == 1
        Func_Type = 'THzOn-NoCost';
        [F_HO_AveCh_err_TOn_NoCo(i),F_HO_AveRate_err_TOn_NoCo(i),F_HO_AveCh_GT_TOn_NoCo(i),F_HO_AveRate_GT_TOn_NoCo(i),NumHO_TOn_NoCo,W_MAT_TOn_NoCo,alpha_mat_TOn_NoCo,~] = THzOn_Mobile_Func(Traj_data,K,B_R,Channel_info,...
            H_min,General_Struct,B_T,M,eps_conv,Pmax_THz,Rth_vec,THz_BW,Clus_Siz_MBN,THz_Hyb,Func_Type,eta_BL,...
        Max_Min_Init,eta_THz,HO_Pen_1,HO_Pen_2,Ch_err,Users_Init_x,All_BS_Loc,seed);
        NuHO_TOn_NoCo(i) = NumHO_TOn_NoCo;
    else
        [F_HO_AveCh_err_TOn_NoCo(i),F_HO_AveRate_err_TOn_NoCo(i),F_HO_AveCh_GT_TOn_NoCo(i),F_HO_AveRate_GT_TOn_NoCo(i)] = THzOn_Mobile_Func_HO(Traj_data,K,B_R,Channel_info,...
            H_min,General_Struct,B_T,M,THz_BW,W_MAT_TOn_NoCo,alpha_mat_TOn_NoCo,THz_Hyb,eta_BL,...
        eta_THz,Ch_err,Users_Init_x,All_BS_Loc,Rth_vec,seed);
        NuHO_TOn_NoCo(i) = NumHO_TOn_NoCo;

    end



end

end



