% The number of sub-bands varies:
clc
clear all
close all
% seed = 20;
% rng(seed)

TW_Tolerance = 0.006;
[w_I,w_e,f_I,f_e,k_f] = EdgeBands(TW_Tolerance);
Fractional_BW = 1e-2;
betta = 2;
Guard_Band = 1e8;
[S_hat_L,S_hat_U,w_T] = Sub_bands_Func(w_I,w_e,f_I,f_e,Fractional_BW,betta,Guard_Band);

% Loading parameters:
G_Rx = 25; %[dB]
G_Rx = db2pow(G_Rx);
G_Tx = 25; %[dB]
G_Tx = db2pow(G_Tx);
c = physconst('LightSpeed');
C = (G_Tx*G_Rx*c^2)/((4*pi)^2);
SysParams.C = C;
k_B = 1.381e-23;  % Boltzmann constant
T0 = 290; % noise temperature (Kelvin)
Ther_Noise_PSD = k_B*T0; % N0
alpha = -2; % Pathloss exponent
S = S_hat_L;
% Loading parameters:
G_Rx = 25; %[dB]
G_Rx = db2pow(G_Rx);
G_Tx = 25; %[dB]
G_Tx = db2pow(G_Tx);
c = physconst('LightSpeed');
C = (G_Tx*G_Rx*c^2)/((4*pi)^2);
SysParams.C = C;
k_B = 1.381e-23;  % Boltzmann constant
T0 = 290; % noise temperature (Kelvin)
Ther_Noise_PSD = k_B*T0; % N0
alpha = -2; % Pathloss exponent
SystemParams.alpha = alpha;
SystemParams.KF = k_f;


S = S_hat_L;
%S = 20;
B = 10;
N_Orig = 12;

% Generating users and BSs locations:


[w,F_S_vec] = BW_Freq_Vec_Func(S,f_I,w_T,w_e,w_I,Guard_Band);
SystemParams.F_S_vec = F_S_vec;
SystemParams.f_THz = median(F_S_vec);
SystemParams.c = c;

SystemParams.w = w;

OptParams.Init_scale = 0.1;
OptParams.FP_Iterations = 8000;
OptParams.T_main_Direct_FP = 100;
OptParams.eps_Power = 1e-3;
OptParams.eps_main_FP = 1e-4;
BiSecParam.mu_b_Min = 0;
BiSecParam.mu_b_Max = 1e4;
BiSecParam.errTol_Bisec = 1e-7;

OptParams.Outer_ADMM_Iterations = 8000;
OptParams.Outer_ADMM_Accuracy = 1e-3;
OptParams.Inner_ADMM_Iterations = 3000;
OptParams.Inner_ADMM_Accuracy = 1e-3;
OptParams.ADMM_Penalty_Fac = 2.2;

Assoc_params.Idle_BS_flg = 1;
Assoc_params.Gam_U = S;

SystemParams.InvAntennaDirec = 0.5;
SystemParams.Pmax = 1;
Assoc_params.Gam_L = 1;
MNT = 4;
HI_Vec = linspace(0,0.7,5); 

for mnt = 1 : MNT
    mnt
    seed = mnt+5;
    rng(seed)
    
    Ther_Noise = (Ther_Noise_PSD * w .* (F_S_vec).^2)/C;
    Ther_Noise = reshape(Ther_Noise,[1,S,1]);

    [WSR_main_FP_Rnd_Assoc,Mean_Conn_Order_Rnd_Assoc,WSR_EqPow,Mean_Conn_Order_Eq_Pow,...
    WSR_main_FP,Mean_Conn_Order_FP,WSR_Main_ADMM,Mean_Conn_Order_ADMM,WSR_main_FP_Sing_Conn,Mean_Conn_Order_FP_Sing_Conn] = ...
    InnerFunc_HWI_Varies(HI_Vec,Assoc_params,S,N_Orig,Ther_Noise,OptParams,BiSecParam,...
        SystemParams,B,seed);
    
    WSR_main_FP_Rnd_Assoc_MNT(mnt,:,:) = WSR_main_FP_Rnd_Assoc;
    Mean_Conn_Order_Rnd_Assoc_MNT(mnt,:,:) = Mean_Conn_Order_Rnd_Assoc;
    WSR_EqPow_MNT(mnt,:,:) = WSR_EqPow;
    Mean_Conn_Order_Eq_Pow_MNT(mnt,:,:) = Mean_Conn_Order_Eq_Pow;
    WSR_main_FP_MNT(mnt,:,:) = WSR_main_FP;
    Mean_Conn_Order_FP_MNT(mnt,:,:) = Mean_Conn_Order_FP;
    WSR_Main_ADMM_MNT(mnt,:,:) = WSR_Main_ADMM;
    Mean_Conn_Order_ADMM_MNT(mnt,:,:) = Mean_Conn_Order_ADMM;
    WSR_Main_Sing_Conn_MNT(mnt,:,:) = WSR_main_FP_Sing_Conn;
    Mean_Conn_Order_Sing_Conn_MNT(mnt,:,:) = Mean_Conn_Order_FP_Sing_Conn;

end

WSR_main_FP_Rnd_Assoc_AL = squeeze(mean(WSR_main_FP_Rnd_Assoc_MNT,1)); 
Mean_Conn_Order_Rnd_Assoc_AL = squeeze(mean(Mean_Conn_Order_Rnd_Assoc_MNT,1));
WSR_EqPow_AL = squeeze(mean(WSR_EqPow_MNT,1));
Mean_Conn_Order_Eq_Pow_AL = squeeze(mean(Mean_Conn_Order_Eq_Pow_MNT,1));
WSR_main_FP_AL = squeeze(mean(WSR_main_FP_MNT,1));
Mean_Conn_Order_FP_AL = squeeze(mean(Mean_Conn_Order_FP_MNT,1));
WSR_Main_ADMM_AL = squeeze(mean(WSR_Main_ADMM_MNT,1));
Mean_Conn_Order_ADMM_AL = squeeze(mean(Mean_Conn_Order_ADMM_MNT,1));
WSR_Main_Sing_Conn_AL = squeeze(mean(WSR_Main_Sing_Conn_MNT,1));
Mean_Conn_Order_Sing_Conn_AL = squeeze(mean(Mean_Conn_Order_Sing_Conn_MNT,1));




Color1 = [0, 0, 1]; % blue
Color2 = [0, 0.5, 0]; % green
Color3 = [1, 0, 0]; % red 
Color7 = [0, 0.75, 0.75]; % cyan
Color5 = [0.75, 0.75, 0];
Color6 = [0.75, 0, 0.75]; % purple
Color4 = [0.25, 0.25, 0.25]; % black
Color8 = [0.8500, 0.3250, 0.0980]; % dark orange

Gb_factor = (1/1e9);
figure(1)
box on
plot(HI_Vec,Gb_factor*WSR_main_FP_AL(1,:),'-*','Color',Color1,'MarkerIndices',1:2:length(HI_Vec),'LineWidth',0.8)
hold on
plot(HI_Vec,Gb_factor*WSR_Main_ADMM_AL(1,:),'-*','Color',Color2,'MarkerIndices',1:2:length(HI_Vec),'LineWidth',0.8)
hold on
plot(HI_Vec,Gb_factor*WSR_EqPow_AL(1,:),'-*','Color',Color4,'MarkerIndices',1:2:length(HI_Vec),'LineWidth',0.8)
hold on
plot(HI_Vec,Gb_factor*WSR_main_FP_Rnd_Assoc_AL(1,:),'-*','Color',Color6,'MarkerIndices',1:2:length(HI_Vec),'LineWidth',0.8)
hold on
plot(HI_Vec,Gb_factor*WSR_Main_Sing_Conn_AL(1,:),'-s','Color',Color3,'MarkerIndices',1:2:length(HI_Vec),'LineWidth',0.8)


grid on

legend({'Algo3-Cent. PA','Algo3-Dist. PA',...
    'Eq. PA+Opt. UASA','Cent. PA+Rnd. UASA',...
    'Algo3-Cent. PA-SC'},'FontSize',7)
ylabel('System Sum Rate [Gbps]', 'FontSize', 10)
xlabel('Levels of hardware imperfections at BSs and users, k_t=k_r', 'FontSize', 10)



xlabel('Levels of hardware imperfections at BSs and users, k_t=k_r', 'FontSize', 10)


