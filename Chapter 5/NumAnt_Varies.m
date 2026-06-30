
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

v = 40; % Users velocity


Traj_data.Trajec_Time = 0.3;
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
K = 8; % Number of users
B_T = 4; % Number of TBSs
B_R = 2; % Number of RBSs

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
Channel_info.W_min_THz = 30;

%%% Parameters Loading:
eta_BL = 0e-3;

THz_Ph_fac = 0.1;
RF_Ph_fac = 0.1;
Clus_Siz_THz = 20;
Clus_Siz_RF = 20;
Clus_Siz_MBN = 40;

Rth = 0.1e9;
Rth_vec = Rth*ones(1,K);

Log_Util = 0;
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

General_Struct.S_HO = 1;

Ch_err = 0;

General_Struct.Mobil_flg = 0;

General_Struct.MRT_Init_We = 10;


eta_THz = 0;
eta_RF = 0;


[THz_BS_Loc,RF_BS_Loc,All_BS_Loc] = BS_Pos(B_T,B_R,Channel_info);
L_x_max_Users = 135;


M_vec = [50,100,200,300,400,500,600,700,800,1000];

M_hat_vec = [20,50,100,150,200,250,300,350,400,500];

v_vec = [0,40];

x_vec = [0.5, 1];


MNT = 50*6;
Seed_Vec = 1:MNT;



T_Inp = 30;
T_MaxMin = 15;
Max_Min_Init = 0;
General_Struct.MaxMin_flg = 1;

for mnt = Seed_Vec
    seed = Seed_Vec(mnt);
    mnt
    [SumRate_HO_MBN_NF_NF(:,:,:,mnt),SumRate_HO_THz_NF_NF(:,:,:,mnt),SumRate_HO_RF_NF_NF(:,:,:,mnt),Infeas_flg_NF_NF(:,:,:,mnt),Tmax_NF_NF(:,:,:,mnt),...
        SumRate_HO_MBN_FF_NF(:,:,:,mnt),SumRate_HO_THz_FF_NF(:,:,:,mnt),SumRate_HO_RF_FF_NF(:,:,:,mnt),Infeas_flg_FF_NF(:,:,:,mnt),Tmax_FF_NF(:,:,:,mnt)]...
    = Inner_Func(Traj_data,K,B_R,Channel_info,...
        H_min,General_Struct,B_T,eps_conv,Pmax_THz,THz_BW,Clus_Siz_THz,Clus_Siz_RF,Clus_Siz_MBN,eta_BL,...
    Pmax_RF,RF_BW_Norm,RF_BW,Max_Min_Init,eta_THz,eta_RF,Rth_vec,L_x_max_Users,THz_BS_Loc,RF_BS_Loc,...
    T_Inp,T_MaxMin,v_vec,x_vec,M_vec,M_hat_vec,seed);
end

SumRate_HO_MBN_NF_NF_Mean = mean(SumRate_HO_MBN_NF_NF,4);
SumRate_HO_THz_NF_NF_Mean = mean(SumRate_HO_THz_NF_NF,4);
SumRate_HO_RF_NF_NF_Mean = mean(SumRate_HO_RF_NF_NF,4);

SumRate_HO_MBN_FF_NF_Mean = mean(SumRate_HO_MBN_FF_NF,4);
SumRate_HO_THz_FF_NF_Mean = mean(SumRate_HO_THz_FF_NF,4);
SumRate_HO_RF_FF_NF_Mean = mean(SumRate_HO_RF_FF_NF,4);


Color1 = [0, 0, 1]; % blue
Color2 = [0, 0.5, 0]; % green
Color3 = [1, 0, 0]; % red 
Color4 = [0.25, 0.25, 0.25]; % black
Color7 = [0, 0.75, 0.75]; % cyan
Color5 = [0.75, 0.75, 0];
Color6 = [0.75, 0, 0.75]; % purple
Color8 = [0.8500, 0.3250, 0.0980]; % dark orange


Gb_factor = (1/1e9);
figure(1)
box on

plot(M_vec, Gb_factor * squeeze(SumRate_HO_MBN_NF_NF_Mean(1,1,:)), '-', 'Color', Color1, 'MarkerIndices', 1:1:length(M_vec), 'LineWidth', 0.8)
hold on
plot(M_vec, Gb_factor * squeeze(SumRate_HO_MBN_NF_NF_Mean(1,2,:)), '-*', 'Color', Color2, 'MarkerIndices', 1:1:length(M_vec), 'LineWidth', 0.8)
hold on
plot(M_vec, Gb_factor * squeeze(SumRate_HO_MBN_NF_NF_Mean(2,1,:)), '-s', 'Color', Color4, 'MarkerIndices', 1:1:length(M_vec), 'LineWidth', 0.8)
hold on
plot(M_vec, Gb_factor * squeeze(SumRate_HO_MBN_NF_NF_Mean(2,2,:)), '-o', 'Color', Color3, 'MarkerIndices', 1:1:length(M_vec), 'LineWidth', 0.8)
hold on

plot(M_vec, Gb_factor * squeeze(SumRate_HO_MBN_FF_NF_Mean(1,1,:)), '--', 'Color', Color1, 'MarkerIndices', 1:1:length(M_vec), 'LineWidth', 0.8)
hold on
plot(M_vec, Gb_factor * squeeze(SumRate_HO_MBN_FF_NF_Mean(1,2,:)), '--*', 'Color', Color2, 'MarkerIndices', 1:1:length(M_vec), 'LineWidth', 0.8)
hold on
plot(M_vec, Gb_factor * squeeze(SumRate_HO_MBN_FF_NF_Mean(2,1,:)), '--s', 'Color', Color4, 'MarkerIndices', 1:1:length(M_vec), 'LineWidth', 0.8)
hold on
plot(M_vec, Gb_factor * squeeze(SumRate_HO_MBN_FF_NF_Mean(2,2,:)), '--o', 'Color', Color3, 'MarkerIndices', 1:1:length(M_vec), 'LineWidth', 0.8)


ylabel('Sum Rate [Gbps]')
xlabel('Number of THz Antennas', 'FontSize', 9.5)

grid on

ax1 = gca;


set(ax1, 'XTick', M_vec, 'XTickLabel', M_vec,'FontSize', 9);

ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', ...
           'Color', 'none', ...
           'XColor', 'k', ...
           'YColor', 'none','FontSize', 9);

xlabel(ax2, 'Number of UMB Antennas', 'FontSize', 9)

linkaxes([ax1, ax2], 'x');



set(ax2, 'XTick', M_vec, 'XTickLabel', M_hat_vec);

uistack(ax1, 'top');
legend(ax1, {'NF-NF, v=0[m/s], x=\lambda/2', ...
             'NF-NF, v=0[m/s], x=\lambda', ...
             'NF-NF, v=40[m/s], x=\lambda/2', ...
             'NF-NF, v=40[m/s], x=\lambda', ...
             'FF-NF, v=0[m/s], x=\lambda/2', ...
             'FF-NF, v=0[m/s], x=\lambda', ...
             'FF-NF, v=40[m/s], x=\lambda/2', ...
             'FF-NF, v=40[m/s], x=\lambda'}, ...
            'FontSize', 8, 'Location', 'northwest');

xlim([min(M_vec), max(M_vec)])



function [SumRate_HO_MBN_NF_NF,SumRate_HO_THz_NF_NF,SumRate_HO_RF_NF_NF,Infeas_flg_NF_NF,Tmax_NF_NF,...
        SumRate_HO_MBN_FF_NF,SumRate_HO_THz_FF_NF,SumRate_HO_RF_FF_NF,Infeas_flg_FF_NF,Tmax_FF_NF]...
    = Inner_Func(Traj_data,K,B_R,Channel_info,...
        H_min,General_Struct,B_T,eps_conv,Pmax_THz,THz_BW,Clus_Siz_THz,Clus_Siz_RF,Clus_Siz_MBN,eta_BL,...
    Pmax_RF,RF_BW_Norm,RF_BW,Max_Min_Init,eta_THz,eta_RF,Rth_vec,L_x_max_Users,THz_BS_Loc,RF_BS_Loc,...
    T_Inp,T_MaxMin,v_vec,x_vec,M_vec,M_hat_vec,seed)


        
for i = 1 : length(v_vec)
    
    v = v_vec(i);

    for j = 1 : length(x_vec)
        x = x_vec(j);

        for k = 1:length(M_vec)
            k;
            M = M_vec(k);
            M_hat = M_hat_vec(k);  
            
            
            rng(seed)
            Users_in_X = rand(K,1);
            Users_Init_x = (L_x_max_Users)*Users_in_X;
            
            
            N = Traj_data.N;
            
            tau_eps_GT = Traj_data.tau_eps_GT;
            N_GT = Traj_data.N_GT;
            N_RF_Coh = Traj_data.N_RF_Coh;
            
            
            
            rng(seed)
            
            RF_Rayleigh_MAT_TEMP = 1/sqrt(2) * (randn(M_hat,B_R,K,ceil(N_GT/N_RF_Coh),N) + 1j*randn(M_hat,B_R,K,ceil(N_GT/N_RF_Coh),N));
            RF_Ray_idx = repelem(1:ceil(N_GT/N_RF_Coh), N_RF_Coh);
            RF_Rayleigh_mat_TEMP = RF_Rayleigh_MAT_TEMP(:,:,:,RF_Ray_idx,:);
            
            
            alpha_n_1 = ones(B_T,K);
            beta_n_1 = ones(B_R,K);
            
            n=1;
            UA_flg = 1;
            
            
            
            RF_Rayleigh_mat_GT = RF_Rayleigh_mat_TEMP(:,:,:,:,n);
            
            [Users_x_mat_GT,Users_Init_y] = Users_pos(K,THz_BS_Loc,N_GT,tau_eps_GT,v,Users_Init_x,Channel_info);
            
            [THz_User_Dist_GT,AoD_THz_Users_Angle_GT] = BS_Users_Dis_ang(Users_x_mat_GT,Users_Init_y,THz_BS_Loc,H_min);
            
            Blockage_Prob_MBN = exp(-THz_User_Dist_GT(:,:,1)*eta_BL);
            rng(seed*n)
            BL_mat = rand(size(THz_User_Dist_GT(:,:,1))) < Blockage_Prob_MBN;
            
            
            %%%%% NFC:
            
            NFC_flg = 1;
            [PL_THz_GT,PL_tilde_THz_GT,~,H_NF,~,~] = THz_Channel(x,v,K,B_T,M,THz_User_Dist_GT,AoD_THz_Users_Angle_GT,Traj_data,General_Struct,NFC_flg);
            
            F_Analog_NF = exp(1j*angle(H_NF))/sqrt(M);
            [Eff_H_mat_NF_NF,Eff_H_tilde_mat_NF_NF,~] = Eff_Channel_THz_V4(F_Analog_NF,H_NF,PL_THz_GT,PL_tilde_THz_GT,B_T,K,M,General_Struct,Traj_data);
            
            H_mat_NF_NF = Eff_H_mat_NF_NF(:,:,:,1);
            H_tilde_mat_NF_NF = Eff_H_tilde_mat_NF_NF(:,:,:,1);
            
            
            
            [RF_User_Dist_GT,AoD_RF_Users_Angle_GT] = BS_Users_Dis_ang(Users_x_mat_GT,Users_Init_y,RF_BS_Loc,H_min);
            
            NFC_flg = 1;
            [PL_RF_GT,~,G_NF,~] = RF_Channel(x,v,K,B_R,M_hat,RF_User_Dist_GT,AoD_RF_Users_Angle_GT,RF_Rayleigh_mat_GT,Traj_data,General_Struct,NFC_flg);
            
            
            Q_Analog_NF = exp(1j*angle(G_NF))/sqrt(M_hat);
            [Eff_G_mat_NF_NF,~] = Eff_Channel_RF_V4(Q_Analog_NF,G_NF,PL_RF_GT,B_R,M_hat,K,General_Struct,Traj_data);
            
            G_mat_NF_NF = Eff_G_mat_NF_NF(:,:,:,1);
            
            
            
            [SumRate_HO_MBN_NF_NF(i,j,k),~,SumRate_HO_THz_NF_NF(i,j,k),SumRate_HO_RF_NF_NF(i,j,k),~,~,Infeas_flg_NF_NF(i,j,k),Tmax_NF_NF(i,j,k),~,~,~,~,~,~] = ...
            MBN_Opt_Stat_Func(T_Inp,T_MaxMin,General_Struct,Traj_data,eps_conv,Pmax_THz,H_mat_NF_NF,H_tilde_mat_NF_NF,B_T,K,M,Rth_vec,THz_BW,Clus_Siz_THz,UA_flg,...
            Pmax_RF,G_mat_NF_NF,M_hat,B_R,RF_BW_Norm,Clus_Siz_RF,Clus_Siz_MBN,RF_BW,Max_Min_Init,BL_mat,eta_THz,eta_RF,alpha_n_1,beta_n_1);
            
            %%%%% FFC:
            
            NFC_flg = 0;
            [~,~,~,H_FF,~,~] = THz_Channel(x,v,K,B_T,M,THz_User_Dist_GT,AoD_THz_Users_Angle_GT,Traj_data,General_Struct,NFC_flg);
            
            F_Analog_FF = exp(1j*angle(H_FF))/sqrt(M);
            [Eff_H_mat_FF_NF,Eff_H_tilde_mat_FF_NF,~] = Eff_Channel_THz_V4(F_Analog_FF,H_NF,PL_THz_GT,PL_tilde_THz_GT,B_T,K,M,General_Struct,Traj_data);
            
            H_mat_FF_NF = Eff_H_mat_FF_NF(:,:,:,1);
            H_tilde_mat_FF_NF = Eff_H_tilde_mat_FF_NF(:,:,:,1);
            
            
            
            NFC_flg = 0;
            [~,~,G_FF,~] = RF_Channel(x,v,K,B_R,M_hat,RF_User_Dist_GT,AoD_RF_Users_Angle_GT,RF_Rayleigh_mat_GT,Traj_data,General_Struct,NFC_flg);
            
            
            Q_Analog_FF = exp(1j*angle(G_FF))/sqrt(M_hat);
            [Eff_G_mat_FF_NF,~] = Eff_Channel_RF_V4(Q_Analog_FF,G_NF,PL_RF_GT,B_R,M_hat,K,General_Struct,Traj_data);
            
            G_mat_FF_NF = Eff_G_mat_FF_NF(:,:,:,1);
            
            
            
            [SumRate_HO_MBN_FF_NF(i,j,k),~,SumRate_HO_THz_FF_NF(i,j,k),SumRate_HO_RF_FF_NF(i,j,k),~,~,Infeas_flg_FF_NF(i,j,k),Tmax_FF_NF(i,j,k),~,~,~,~,~,~] = ...
            MBN_Opt_Stat_Func(T_Inp,T_MaxMin,General_Struct,Traj_data,eps_conv,Pmax_THz,H_mat_FF_NF,H_tilde_mat_FF_NF,B_T,K,M,Rth_vec,THz_BW,Clus_Siz_THz,UA_flg,...
            Pmax_RF,G_mat_FF_NF,M_hat,B_R,RF_BW_Norm,Clus_Siz_RF,Clus_Siz_MBN,RF_BW,Max_Min_Init,BL_mat,eta_THz,eta_RF,alpha_n_1,beta_n_1);
    
        end
    end
end

end



