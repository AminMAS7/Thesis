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
Params.fcRF = 2*10^9; % RF carrier frequency (2 GHz)
Params.P_RF = 3; % RF transmission power [Watt]
Params.THz_Tx_gain = 316.228; % Transmitter antenna gain for TBS
Params.THz_Rx_gain = 316.228; % Reciever antenna gain for TBS
Params.fcTHz = 1e12; % Terahertz carrier frequency (1 THz)
Params.P_THz = 0.2; % THz transmission power [Watt]
Params.RF_Bandwidth = 40e6; % RF bandwidth (40 MHz)
%Params.THz_Bandwidth = 10e9; % THz bandwidth (10 GHz)
Params.kf = Kav1(Params.fcTHz); % Molecular absorption coefficient at f_THz = 1 THz 
Params.MissAll_Probability = 0.006; % Probabely the random variable D ( Comes from the alignment between the user and interfering TBSs)
k_B = 1.381e-23;   % Boltzmann constant
T0 = 290;   % noise temperature (Kelvin)
Thermal_Noise_RF = ( Params.RF_Bandwidth )*k_B*T0;
Params.Thermal_Noise_RF = Thermal_Noise_RF;

% Monte Carlo and Positions loading
R_max = 400; % Maximum radius of the considered circle
Params.R_max = R_max;
A = pi*R_max^2; % Area of considered environment
Params.Area = A;
T = 40e3; % Number of Monte-Carlo iterations
Nmax = 2000;

tic
R_th_vec = linspace(0.1,1000,100).*1e6;
R_th_vec = R_th_vec(1:51);

THz_BW_vec = [1e9, 10e9];
Intg_RF_Ass_mat = zeros(length(THz_BW_vec),length(R_th_vec));
Intg_THz_Ass_mat = zeros(length(THz_BW_vec),length(R_th_vec));

SA_RF_Ass_mat = zeros(length(THz_BW_vec),length(R_th_vec));
SA_THz_Ass_mat = zeros(length(THz_BW_vec),length(R_th_vec));

N_Intg_selected = zeros(length(THz_BW_vec),length(R_th_vec));
N_SA_selected = zeros(length(THz_BW_vec),length(R_th_vec));

mean_rate_SA_mat = zeros(length(THz_BW_vec),length(R_th_vec));
Num_RF = zeros(length(THz_BW_vec),length(R_th_vec));
Num_THz = zeros(length(THz_BW_vec),length(R_th_vec));
max_rate_Intg_mean_mat = zeros(length(THz_BW_vec),length(R_th_vec));
N_SA1_flg = zeros(length(THz_BW_vec),length(R_th_vec));

Power_flg = 0;

THz_BW_idx = 0;
for THz_Bandwidth = THz_BW_vec
    Params.THz_Bandwidth = THz_Bandwidth;
    THz_BW_idx = THz_BW_idx + 1
    R_idx = 0;
    for R_th = R_th_vec
        R_idx = R_idx + 1
        % Integrated Deployment:
        Intg_Stop_flg = 0;
        for NInteg = 2:Nmax
            NInteg;
            if Intg_Stop_flg == 1
                break
            end
            [max_rate_Intg_mean,RFassoc_Intg_mean,THassoc_Intg_mean] = Monte_Carlo_Integ(NInteg,T,Params,Power_flg);
            max_rate_Intg_mean_mat(THz_BW_idx,R_idx) = max_rate_Intg_mean;
            if max_rate_Intg_mean >= R_th
               Intg_Stop_flg = 1;
               N_Intg_selected(THz_BW_idx,R_idx) = NInteg;
               Intg_RF_Ass_mat(THz_BW_idx,R_idx) = RFassoc_Intg_mean;
               Intg_THz_Ass_mat(THz_BW_idx,R_idx) = THassoc_Intg_mean;
            end
        end
    
        % SA Deployment:
        SA_Stop_flg = 0;
    
    
        for N_SA = 2:Nmax
    
            if SA_Stop_flg == 1
                break
            end
    
            SA_idx = 0;
            rate_SA_mat = zeros(1,N_SA-1);
            for NR = 1:N_SA-1
                for NT = N_SA-1:-1:1
                    if NR + NT == N_SA
                        SA_idx = SA_idx + 1;
                        [max_rate_SA_mean,RFassoc_SA_mean,THassoc_SA_mean] = Monte_Carlo_SA(NR,NT,T,Params,Power_flg);
                        rate_SA_mat(SA_idx) = max_rate_SA_mean;
                        RFassoc_SA_mat(SA_idx) = RFassoc_SA_mean;
                        THassoc_SA_mat(SA_idx) = THassoc_SA_mean;
                    end
                end
            end
    
            mean_rate_SA_mat(THz_BW_idx,R_idx) = max(rate_SA_mat);
            if max(rate_SA_mat) >= R_th
                SA_Stop_flg = 1;
                N_SA_selected(THz_BW_idx,R_idx) = N_SA;
                [max_rate,max_rate_idx] = max(rate_SA_mat);
                Num_RF(THz_BW_idx,R_idx) = max_rate_idx;
                Num_THz(THz_BW_idx,R_idx) = N_SA-max_rate_idx;
            end
        end
    N_SA_selected = N_SA_selected
    N_Intg_selected = N_Intg_selected    
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


figure(1)
box on
plot(R_th_vec./1e6 , N_SA_Equal(1,:),'--','Color',Color1,'LineWidth',0.8);
hold on
plot(R_th_vec./1e6 , N_Intg_selected(1,:),'-','Color',Color1,'LineWidth',0.8);
hold on
plot(R_th_vec./1e6 , N_SA_selected(1,:),':','Color',Color1,'LineWidth',0.8);
hold on 
plot(R_th_vec./1e6 , N_SA_Equal(2,:),'--','Color',Color2,'LineWidth',0.8);
hold on
plot(R_th_vec./1e6 , N_SA_selected(2,:),':','Color',Color2,'LineWidth',0.8);
hold on 
plot(R_th_vec./1e6 , N_Intg_selected(2,:),'-','Color',Color2,'LineWidth',0.8);

xlabel('Target Rate Threshold [Mbps]', 'FontSize', 12)
ylabel('Number of Required BSs', 'FontSize', 12)
legend({'SA-MBN- B_T=1 GHz','Int-MBN- B_T=1 GHz','SA-MBN-FN- B_T=1 GHz','SA-MBN- B_T=10 GHz','SA-MBN-FN- B_T=10 GHz','Int-MBN- B_T=10 GHz'},'FontSize',7)
set(gca, 'YScale', 'log')
xlim([min(R_th_vec./1e6) max(R_th_vec./1e6)+10])
ylim([2 max(N_SA_Equal(1,:))])
grid on
