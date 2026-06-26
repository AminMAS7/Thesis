% Author: Mohammad Amin Saeidi
% Date: August 2023
%
% Usage and Redistribution Notice:
% - This code is provided for [educational/research/personal] purposes only.
% - Any use, modification, redistribution, or other use of this code
%   requires explicit permission from the author.
% - To request permission, please send an email to [amin96a@yorku.ca or m.amin.saeidi96@gmail.com].

clc
close all
clear all

% Loading Parameters:
Params.RF_Tx_gain = 1; % Transmitter antenna gain for RBS
Params.RF_Rx_gain = 1; % Reciever antenna gain for RBS
Params.fcRF = 2e9; % RF carrier frequency (2 GHz)
Params.P_RF = 2; % RF transmission power [Watt]
Params.THz_Tx_gain = 316.228; % Transmitter antenna gain for TBS
Params.THz_Rx_gain = 316.228; % Reciever antenna gain for TBS
Params.fcTHz = 1e12; % Terahertz carrier frequency (1 THz)
Params.P_THz = 0.2; % THz transmission power [Watt]
Params.RF_Bandwidth = 40e6; % RF bandwidth (40 MHz)
Params.THz_Bandwidth = 10e9; % THz bandwidth (10 GHz)
Kav1(1e12)
%load("Molecular_absorption.mat","BW","K_abs");
%kf = K_abs(BW == 0.8375e12);
%Params.kf = kf; % Molecular absorption coefficient at f_THz = 1 THz 

Params.MissAll_Probability = 0.006; % Probabely the random variable D ( Comes from the alignment between the user and interfering TBSs)
k_B = 1.381e-23;   % Boltzmann constant
T0 = 290;   % noise temperature (Kelvin)
Thermal_Noise_RF = ( Params.RF_Bandwidth )*k_B*T0;
Params.Thermal_Noise_RF = Thermal_Noise_RF;

Thermal_Noise_THz = ( Params.THz_Bandwidth )*k_B*T0;
Params.Thermal_Noise_THz = Thermal_Noise_THz;

% Monte Carlo and Positions loading
R_max = 400; % Maximum radius of the considered circle
A = pi*R_max^2; % Area of considered environment
Xu = 0; % Position of the user of interest
Yu = 0; % Position of the user of interest
T = 5e5; % Number of Monte-Carlo iterations
Nmax = 40; % Total number of BSs
N = Nmax;
% N_vec = [1:Nmax];
Kf_vec = linspace(0,0.3,100);
N_vec = [10,60];
tic

for j = [1,2]
    N = N_vec(j);
    idx = 0;
    for kf = Kf_vec
        Params.kf = kf;
        idx = idx +1
        parfor i = 1 : T
            % Loading BSs positions for Stand-alone deployment
            % THz:
            rTHzbs = R_max*sqrt(rand(N,1));
            thetaTb = 2*pi*rand(N,1);
            Xb = rTHzbs.*cos(thetaTb);
            Yb = rTHzbs.*sin(thetaTb);
            % RF:
            rRFzbs = R_max*sqrt(rand(N,1));
            thetaRb = 2*pi*rand(N,1);
            Xrb = rRFzbs.*cos(thetaRb);
            Yrb = rRFzbs.*sin(thetaRb);
            % Distacnces from THz BSs
            [Xmp_Tmat, Xp_Tmat] = meshgrid(Xu,Xb);
            [Ymp_Tmat, Yp_Tmat] = meshgrid(Yu,Yb);
            D_ue_Tbs = sqrt((Xmp_Tmat-Xp_Tmat).^2 + (Ymp_Tmat-Yp_Tmat).^2);
            % Distacnces from RF BSs
            [Xmp_Rmat, Xp_Rmat] = meshgrid(Xu,Xrb);
            [Ymp_Rmat, Yp_Rmat] = meshgrid(Yu,Yrb);
            D_ue_Rbs = sqrt((Xmp_Rmat-Xp_Rmat).^2 + (Ymp_Rmat-Yp_Rmat).^2);
            % Fading channel for RF BSs
            fadeRand = exprnd(1,N,1); % Channel fade RV based on number of users of interest (1) and number of RF BSs
    
            % Loading BSs positions for Integrated deployment
            rRF_THzbs = R_max*sqrt(rand(N,1));
            thetaRF_TH = 2*pi*rand(N,1);
            XbHyb = rRF_THzbs.*cos(thetaRF_TH);
            YbHyb = rRF_THzbs.*sin(thetaRF_TH);
    
            % Distacnces from RF_THz BSs
            [Xmp_Tmat_Integ, Xp_Tmat_Integ] = meshgrid(Xu,XbHyb);
            [Ymp_Tmat_Integ, Yp_Tmat_Integ] = meshgrid(Yu,YbHyb);
            D_ue_RF_THzbs = sqrt((Xmp_Tmat_Integ-Xp_Tmat_Integ).^2 + (Ymp_Tmat_Integ-Yp_Tmat_Integ).^2);
            % Fading channels for RF mode of Hybrid BSs
            fadeRand_Integ = exprnd(1,N,1);
    
            % Computing rates of both deployments
            % Stand-alone: 
            [Output_SA] = Stand_alone_Rate(D_ue_Rbs,D_ue_Tbs,Params,fadeRand);
    
            max_rate_Upper_SA(i) = Output_SA.max_rate_Upper;
            RFassoc_SA(i) = Output_SA.RFassoc;
            THassoc_SA(i) = Output_SA.THassoc;
            
            max_rate_T1_SA(i) = Output_SA.max_rate_T1;
            RFassoc_T1_SA(i) = Output_SA.RFassoc_T1;
            THassoc_T1_SA(i) = Output_SA.THassoc_T1;
            
            max_rate_SINR_SA(i) = Output_SA.max_rate_SINR;
            RFassoc_SINR_SA(i) = Output_SA.RFassoc_SINR;
            THassoc_SINR_SA(i) = Output_SA.THassoc_SINR;
            
            max_rate_Pow_SA(i) = Output_SA.max_rate_Pow;
            RFassoc_Pow_SA(i) = Output_SA.RFassoc_Pow;
            THassoc_Pow_SA(i) = Output_SA.THassoc_Pow;
    
            % Integrated:
    
            [Output_Intg] = Integ_Rate(D_ue_RF_THzbs,Params,fadeRand_Integ);
    
            max_rate_Upper_Intg(i) = Output_Intg.max_rate_Upper;
            RFassoc_Intg(i) = Output_Intg.RFassoc;
            THassoc_Intg(i) = Output_Intg.THassoc;
            
            max_rate_T1_Intg(i) = Output_Intg.max_rate_T1;
            RFassoc_T1_Intg(i) = Output_Intg.RFassoc_T1;
            THassoc_T1_Intg(i) = Output_Intg.THassoc_T1;
            
            max_rate_SINR_Intg(i) = Output_Intg.max_rate_SINR;
            RFassoc_SINR_Intg(i) = Output_Intg.RFassoc_SINR;
            THassoc_SINR_Intg(i) = Output_Intg.THassoc_SINR;
            
            max_rate_Pow_Intg(i) = Output_Intg.max_rate_Pow;
            RFassoc_Pow_Intg(i) = Output_Intg.RFassoc_Pow;
            THassoc_Pow_Intg(i) = Output_Intg.THassoc_Pow;
    
        end
        % Stand-alone Average rates:
        % Data Rate metric: 
        Net_max_rate_SA_mean(j,idx) = 2 * N * mean(max_rate_Upper_SA);
        RFassoc_SA_mean(j,idx) = mean(RFassoc_SA);
        THassoc_SA_mean(j,idx) = mean(THassoc_SA);
        % Test metric:
        Net_max_rate_SA_T1_mean(j,idx) = 2 * N * mean(max_rate_T1_SA);
        RFassoc_SA_T1_mean(j,idx) = mean(RFassoc_T1_SA);
        THassoc_SA_T1_mean(j,idx) = mean(THassoc_T1_SA);
        % SINR metric:
        Net_max_rate_SA_SINR_mean(j,idx) = 2 * N * mean(max_rate_SINR_SA);
        RFassoc_SA_SINR_mean(j,idx) = mean(RFassoc_SINR_SA);
        THassoc_SA_SINR_mean(j,idx) = mean(THassoc_SINR_SA);
        % Power metric:
        Net_max_rate_SA_Pow_mean(j,idx) = 2 * N * mean(max_rate_Pow_SA);
        RFassoc_SA_Pow_mean(j,idx) = mean(RFassoc_Pow_SA);
        THassoc_SA_Pow_mean(j,idx) = mean(THassoc_Pow_SA);
    
    
        % Integrated Average rates:
        % Data Rate metric: 
        Net_max_rate_Intg_mean(j,idx) = 2 * N * mean(max_rate_Upper_Intg);
        RFassoc_Intg_mean(j,idx) = mean(RFassoc_Intg);
        THassoc_Intg_mean(j,idx) = mean(THassoc_Intg);
        % Test metric:
        Net_max_rate_Intg_T1_mean(j,idx) = 2 * N * mean(max_rate_T1_Intg);
        RFassoc_Intg_T1_mean(j,idx) = mean(RFassoc_T1_Intg);
        THassoc_Intg_T1_mean(j,idx) = mean(THassoc_T1_Intg);
        % SINR metric:
        Net_max_rate_Intg_SINR_mean(j,idx) = 2 * N * mean(max_rate_SINR_Intg);
        RFassoc_Intg_SINR_mean(j,idx) = mean(RFassoc_SINR_Intg);
        THassoc_Intg_SINR_mean(j,idx) = mean(THassoc_SINR_Intg);
        % Power metric:
        Net_max_rate_Intg_Pow_mean(j,idx) = 2 * N * mean(max_rate_Pow_Intg);
        RFassoc_Intg_Pow_mean(j,idx) = mean(RFassoc_Pow_Intg);
        THassoc_Intg_Pow_mean(j,idx) = mean(THassoc_Pow_Intg);
    
    end
end




time = toc


Color1 = [0, 0, 1]; % blue
Color2 = [0, 0.5, 0]; % green
Color3 = [1, 0, 0]; % red 
Color4 = [0, 0.75, 0.75]; % cyan
Color5 = [0.75, 0.75, 0];
Color6 = [0.75, 0, 0.75]; % purple
Color7 = [0.25, 0.25, 0.25]; % black
Color8 = [0.8500, 0.3250, 0.0980]; % dark orange
Kf_vec = linspace(0,0.3,100);



figure(7)

%%%%%%%%%%%%%%%%%%%%%%% N = 60
%figure(2)
box on
plot(Kf_vec , Net_max_rate_SA_mean(2,:)./1e9,'--','Color',Color7,'LineWidth',0.3);
hold on
plot(Kf_vec , Net_max_rate_SA_T1_mean(2,:)./1e9,'--','Color',Color2,'LineWidth',0.3);
hold on
plot(Kf_vec , Net_max_rate_SA_Pow_mean(2,:)./1e9,'--','Color',Color1,'LineWidth',0.3);
hold on
plot(Kf_vec , Net_max_rate_SA_SINR_mean(2,:)./1e9,'--','Color',Color3,'LineWidth',0.3);
hold on
plot(Kf_vec , Net_max_rate_Intg_mean(2,:)./1e9,'-','Color',Color7,'LineWidth',0.3);
hold on
plot(Kf_vec , Net_max_rate_Intg_T1_mean(2,:)./1e9,'-','Color',Color2,'LineWidth',0.3);
hold on
plot(Kf_vec , Net_max_rate_Intg_Pow_mean(2,:)./1e9,'-','Color',Color1,'LineWidth',0.3);
hold on
plot(Kf_vec , Net_max_rate_Intg_SINR_mean(2,:)./1e9,'-','Color',Color3,'LineWidth',0.3);

%%%%%%%%%%%%%%%%%%%%%%% N = 10
hold on
line_fewer_markers(Kf_vec , Net_max_rate_SA_mean(1,:)./1e9,12,'--o','Color',Color7,'LineWidth',0.3,'MarkerSize',3.5,'MarkerEdgeColor',Color7,'MarkerFaceColor','w');
hold on
line_fewer_markers(Kf_vec , Net_max_rate_SA_T1_mean(1,:)./1e9,12,'--o','Color',Color2,'LineWidth',0.3,'MarkerSize',3.5,'MarkerEdgeColor',Color2,'MarkerFaceColor','w');
hold on
line_fewer_markers(Kf_vec , Net_max_rate_SA_Pow_mean(1,:)./1e9,12,'--o','Color',Color1,'LineWidth',0.3,'MarkerSize',3.5,'MarkerEdgeColor',Color1,'MarkerFaceColor','w');
hold on
line_fewer_markers(Kf_vec , Net_max_rate_SA_SINR_mean(1,:)./1e9,12,'--o','Color',Color3,'LineWidth',0.3,'MarkerSize',3.5,'MarkerEdgeColor',Color3,'MarkerFaceColor','w');
hold on
line_fewer_markers(Kf_vec , Net_max_rate_Intg_mean(1,:)./1e9,12,'-o','Color',Color7,'LineWidth',0.3,'MarkerSize',3.5,'MarkerEdgeColor',Color7,'MarkerFaceColor','w');
hold on
line_fewer_markers(Kf_vec , Net_max_rate_Intg_T1_mean(1,:)./1e9,12,'-o','Color',Color2,'LineWidth',0.3,'MarkerSize',3.5,'MarkerEdgeColor',Color2,'MarkerFaceColor','w');
hold on
line_fewer_markers(Kf_vec , Net_max_rate_Intg_Pow_mean(1,:)./1e9,12,'-o','Color',Color1,'LineWidth',0.3,'MarkerSize',3.5,'MarkerEdgeColor',Color1,'MarkerFaceColor','w');
hold on
line_fewer_markers(Kf_vec , Net_max_rate_Intg_SINR_mean(1,:)./1e9,12,'-o','Color',Color3,'LineWidth',0.3,'MarkerSize',3.5,'MarkerEdgeColor',Color3,'MarkerFaceColor','w');


ylim([0 1e5])


set(gca, 'YScale', 'log')
xlabel('Molecular absorption coefficient, k(f)', 'FontSize', 12)
ylabel('Average network data rate [Gbps]', 'FontSize', 12)
legend({'SA-MBN - Max rate Metric, N=60','SA-MBN - Biased Metric, N=60',...
    'SA-MBN - RSRP Metric, N=60','SA-MBN - SINR Metric, N=60',...
    'Int-MBN - Max rate Metric, N=60','Int-MBN - Biased Metric, N=60',...
    'Int-MBN - RSRP Metric, N=60','Int-MBN - SINR Metric, N=60',...
    'SA-MBN - Max rate Metric, N=10','SA-MBN - Biased Metric, N=10',...
    'SA-MBN - RSRP Metric, N=10','SA-MBN - SINR Metric, N=10',...
    'Int-MBN - Max rate Metric, N=10','Int-MBN - Biased Metric, N=10',...
    'Int-MBN - RSRP Metric, N=10','Int-MBN - SINR Metric, N=10'},'Fontsize',7.5)
grid on


x = 29:35;
x2 = Kf_vec(x);
p=axes('position',[.2 .66 .13 .22]);
box on % put box around new pair of axes
plot(Kf_vec , Net_max_rate_SA_mean(2,:)./1e9,'--','Color',Color7,'LineWidth',0.3);
hold on
plot(Kf_vec , Net_max_rate_SA_T1_mean(2,:)./1e9,'--','Color',Color2,'LineWidth',0.3);
hold on
plot(Kf_vec , Net_max_rate_SA_Pow_mean(2,:)./1e9,'--','Color',Color1,'LineWidth',0.3);
hold on
plot(Kf_vec , Net_max_rate_Intg_mean(2,:)./1e9,'-','Color',Color7,'LineWidth',0.3);
hold on
plot(Kf_vec , Net_max_rate_Intg_T1_mean(2,:)./1e9,'-','Color',Color2,'LineWidth',0.3);
set(gca, 'YScale', 'log')

x = 29:35;
x2 = Kf_vec(x);
p=axes('position',[.38 .66 .13 .22]);
box on % put box around new pair of axes
line_fewer_markers(Kf_vec , Net_max_rate_SA_mean(1,:)./1e9,12,'--o','Color',Color7,'LineWidth',0.3,'MarkerSize',3.5,'MarkerEdgeColor',Color7,'MarkerFaceColor','w');
hold on
line_fewer_markers(Kf_vec , Net_max_rate_SA_T1_mean(1,:)./1e9,12,'--o','Color',Color2,'LineWidth',0.3,'MarkerSize',3.5,'MarkerEdgeColor',Color2,'MarkerFaceColor','w');
hold on
line_fewer_markers(Kf_vec , Net_max_rate_SA_Pow_mean(1,:)./1e9,12,'--o','Color',Color1,'LineWidth',0.3,'MarkerSize',3.5,'MarkerEdgeColor',Color1,'MarkerFaceColor','w');
hold on
line_fewer_markers(Kf_vec , Net_max_rate_Intg_mean(1,:)./1e9,12,'-o','Color',Color7,'LineWidth',0.3,'MarkerSize',3.5,'MarkerEdgeColor',Color7,'MarkerFaceColor','w');
hold on
line_fewer_markers(Kf_vec , Net_max_rate_Intg_T1_mean(1,:)./1e9,12,'-o','Color',Color2,'LineWidth',0.3,'MarkerSize',3.5,'MarkerEdgeColor',Color2,'MarkerFaceColor','w');
set(gca, 'YScale', 'log')

