% Author: Mohammad Amin Saeidi
% Date: August 2023
%
% Usage and Redistribution Notice:
% - This code is provided for [educational/research/personal] purposes only.
% - Any use, modification, redistribution, or other use of this code
%   requires explicit permission from the author.
% - To request permission, please send an email to [amin96a@yorku.ca or m.amin.saeidi96@gmail.com].

clc
clear all
close all
 
% Loading Parameters:
Params.RF_Tx_gain = 1; % Transmitter antenna gain for RBS
Params.RF_Rx_gain = 1; % Reciever antenna gain for RBS
Params.fcRF = 2*10^9;  % RF carrier frequency (2 GHz)
Params.P_RF = 2;      % RF transmission power [Watt]
Params.THz_Tx_gain = 316.228; % Transmitter antenna gain for TBS
Params.THz_Rx_gain = 316.228; % Reciever antenna gain for TBS
%Params.fcTHz = 0.8375e12; % Terahertz carrier frequency (1 THz)
Params.fcTHz = 1e12;
Params.P_THz = 0.2;  % THz transmission power [Watt]
Params.RF_Bandwidth = 40e6; % RF bandwidth (40 MHz)
% load("Molecular_absorption.mat","BW","K_abs");
% kf = K_abs(BW == Params.fcTHz);
kf = Kav1(Params.fcTHz);
Params.kf = kf; % Molecular absorption coefficient at f_THz = 1 THz 
Params.MissAll_Probability = 0.006; % Probabely the random variable D ( Comes from the alignment between the user and interfering TBSs)


k_B = 1.381e-23;  % Boltzmann constant
T0 = 290; % noise temperature (Kelvin)
Thermal_Noise_RF = ( Params.RF_Bandwidth )*k_B*T0;
Params.Thermal_Noise_RF = Thermal_Noise_RF;
% Monte Carlo and Positions loading
R_max = 400;  % Maximum radius of the considered circle
Params.R_max = R_max;
A = pi*R_max^2; % Area of considered environment
Params.Area = A;
MNT = 1e6; % Number of Monte-Carlo iterations
T = MNT;


%%%%%%%%% CAPEX and OPEX Parameters Considering small cells %%%%%%%%%
RF_CAPEX = 33e3; % RF BS cost
THz_CAPEX = 38e3; % THz BS cost
Hybrid_CAPEX = 48e3; % Hybrid THz/RF BS cost
RF_OPEX = 2.8e3;
THz_OPEX = 2.7e3;
Hybrid_OPEX = 3.2e3;
RF_Cost = RF_CAPEX + RF_OPEX;
THz_Cost = THz_CAPEX + THz_OPEX;
Hybrid_Cost = Hybrid_CAPEX + Hybrid_OPEX;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Power_flg = 0;
Xu = 0;
Yu = 0;
idx = 0;
N_vec = [1:1:200];
tic
THz_BW_vec = [1e9 10e9];
THz_idx = 1;
for THz_BW = THz_BW_vec
    Params.THz_Bandwidth = THz_BW;
    Thermal_Noise_THz = ( Params.THz_Bandwidth )*k_B*T0;
    Params.Thermal_Noise_THz = Thermal_Noise_THz;
    idx = 0;
    for N = N_vec
        idx = idx +1
        Cost_SA = N * RF_Cost + N* THz_Cost;
        Cost_Intg = N * Hybrid_Cost;

        [Rate_SA,Sum_rate_SA,Sum_SE_SA,SE_SA,RF_ASSOC_SA,THz_ASSOC_SA,Rate_Intg,Sum_rate_Intg,Sum_SE_Intg,SE_Intg,...
            RF_ASSOC_Intg,THz_ASSOC_Intg] = MonteCarlo_func(N,T,Params,R_max);
        % Stand-alone Average rates:

        Rate_SA_mat(THz_idx,idx) = Rate_SA;
        Sum_Rate_SA_mat(THz_idx,idx) = Sum_rate_SA;
        SE_SA_mat(THz_idx,idx) = SE_SA;
        Sum_SE_SA_mat(THz_idx,idx) = Sum_SE_SA;
        DCE_Rate_Based_SA(THz_idx,idx) = Rate_SA/Cost_SA;
        DCE_SE_Based_SA(THz_idx,idx) = Sum_SE_SA/Cost_SA;

        RF_ASSOC_SA_mat(THz_idx,idx) = RF_ASSOC_SA;
        THz_ASSOC_SA_mat(THz_idx,idx) = THz_ASSOC_SA;

        % Integrated Average rates:
        Rate_Intg_mat(THz_idx,idx) = Rate_Intg;
        Sum_Rate_Intg_mat(THz_idx,idx) = Sum_rate_Intg;
        SE_Intg_mat(THz_idx,idx) = SE_Intg;
        Sum_SE_Intg_mat(THz_idx,idx) = Sum_SE_Intg;
        DCE_Rate_Based_Intg(THz_idx,idx) = Rate_Intg/Cost_Intg;
        DCE_SE_Based_Intg(THz_idx,idx) = Sum_SE_Intg/Cost_Intg;
    
        RF_ASSOC_Intg_mat(THz_idx,idx) = RF_ASSOC_Intg;
        THz_ASSOC_Intg_mat(THz_idx,idx) = THz_ASSOC_Intg;
    end
    THz_idx = THz_idx + 1
end
time = toc
Color1 = [0, 0, 1]; % blue
Color2 = [0, 0.5, 0]; % green
Color3 = [1, 0, 0]; % red 
Color7 = [0, 0.75, 0.75]; % cyan
Color5 = [0.75, 0.75, 0];
Color6 = [0.75, 0, 0.75]; % purple
Color4 = [0.25, 0.25, 0.25]; % black
Color8 = [0.8500, 0.3250, 0.0980]; % dark orange

figure(1)
box on
plot(N_vec , DCE_SE_Based_SA(1,:),'--','Color',Color1,'LineWidth',0.8);
hold on
plot(N_vec , DCE_SE_Based_Intg(1,:),'-','Color',Color1,'LineWidth',0.8);
hold on
plot(N_vec , DCE_SE_Based_SA(2,:),'--','Color',Color2,'LineWidth',0.8);
hold on
plot(N_vec , DCE_SE_Based_Intg(2,:),'-','Color',Color2,'LineWidth',0.8);
grid on

xlabel('Number of BSs (N_T = N_R = N_{Hyb})', 'FontSize', 12)
ylabel('Deployment Cost Efficiency [bit/s/Hz/$]', 'FontSize', 12)
legend({'SA-MBN, B_T=1 GHz','Int-MBN, B_T=1 GHz','SA-MBN, B_T=10 GHz','Int-MBN, B_T=10 GHz'},'FontSize',7)

figure(2)
box on
yyaxis left;
plot(N_vec , SE_SA_mat(1,:),'--','Color',Color1,'LineWidth',0.8);
hold on
plot(N_vec , SE_Intg_mat(1,:),'-','Color',Color1,'LineWidth',0.8);
hold on
line_fewer_markers(N_vec , SE_SA_mat(2,:),8,'--o','Color',Color2,'LineWidth',0.8,'MarkerSize',3.5,'MarkerEdgeColor',Color2,'MarkerFaceColor','w')
hold on
line_fewer_markers(N_vec , SE_Intg_mat(2,:),8,'-o','Color',Color2,'LineWidth',0.8,'MarkerSize',3.5,'MarkerEdgeColor',Color2,'MarkerFaceColor','w')
ylabel('Average spectral efficiency [bit/s/Hz]', 'FontSize', 12)
yyaxis right;
plot(N_vec , Rate_SA_mat(1,:)./1e9,'--','Color',Color3,'LineWidth',0.8);
hold on
plot(N_vec , Rate_Intg_mat(1,:)./1e9,'-','Color',Color3,'LineWidth',0.8);
hold on
line_fewer_markers(N_vec , Rate_SA_mat(2,:)./1e9,8,'--o','Color',Color3,'LineWidth',0.8,'MarkerSize',3.5,'MarkerEdgeColor',Color3,'MarkerFaceColor','w')
hold on
line_fewer_markers(N_vec , Rate_Intg_mat(2,:)./1e9,8,'-o','Color',Color3,'LineWidth',0.8,'MarkerSize',3.5,'MarkerEdgeColor',Color3,'MarkerFaceColor','w')
ylabel('Average rate of a typical user [Gbps]', 'FontSize', 12)
grid on
xlabel('Number of BSs (N_T = N_R = N_{Hyb})', 'FontSize', 12)
legend({'SE, SA-MBN, B_T=1 GHz','SE, Int-MBN, B_T=1 GHz','SE, SA-MBN, B_T=10 GHz','SE, Int-MBN, B_T=10 GHz',...
    'Rate, SA-MBN, B_T=1 GHz','Rate, Int-MBN, B_T=1 GHz','Rate, SA-MBN, B_T=10 GHz','Rate, Int-MBN, B_T=10 GHz'},'FontSize',6,'NumColumns',2)


figure(3)
box on
plot(N_vec , THz_ASSOC_SA_mat(1,:),'--','Color',Color1,'LineWidth',0.8);
hold on
plot(N_vec , THz_ASSOC_Intg_mat(1,:),'-','Color',Color1,'LineWidth',0.8);
hold on
plot(N_vec , THz_ASSOC_SA_mat(2,:),'--','Color',Color2,'LineWidth',0.8);
hold on
plot(N_vec , THz_ASSOC_Intg_mat(2,:),'-','Color',Color2,'LineWidth',0.8);

xlabel('Number of BSs (N_T = N_R = N_{Hyb})', 'FontSize', 12)
ylabel('THz Association = (1 - RF Association)', 'FontSize', 12)
legend({'SA-MBN, B_T=1 GHz','Int-MBN, B_T=1 GHz','SA-MBN, B_T=10 GHz','Int-MBN, B_T=10 GHz'},'FontSize',7)
grid on
