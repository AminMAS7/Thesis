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
B = 6;
N_Orig = 12;

% Generating users and BSs locations:


[w,F_S_vec] = BW_Freq_Vec_Func(S,f_I,w_T,w_e,w_I,Guard_Band);
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
% Gam_L = S*B/N_Orig;
% if floor(Gam_L) == Gam_L
%     Gam_L = S*B/N_Orig;
% else
%     Gam_L = ceil(S*B/N_Orig) - 1;
% end
% Assoc_params.Gam_L = Gam_L;
SystemParams.Pmax = 1;
Assoc_params.Gam_L = 1;
MNT = 32*45;
BS_Vec = [2:1:25]; 
parfor mnt = 1 : MNT
    mnt
    seed = mnt+5;
    rng(seed)
    
    Ther_Noise = (Ther_Noise_PSD * w .* (F_S_vec).^2)/C;
    Ther_Noise = reshape(Ther_Noise,[1,S,1]);
    
%     [WSR_main_FP_Rnd_Assoc,Mean_Conn_Order_Rnd_Assoc,WSR_EqPow,Mean_Conn_Order_Eq_Pow,...
%     WSR_main_FP,Mean_Conn_Order_FP,WSR_Main_ADMM,Mean_Conn_Order_ADMM,WSR_main_FP_Sing_Conn,Mean_Conn_Order_FP_Sing_Conn,C3_Violation] =...
%     eta_Varies(eta_Vec,q_Vec,Assoc_params,B,S,N_Orig,Ther_Noise,OptParams,BiSecParam,...
%     SystemParams,seed);
    [WSR_main_FP_Rnd_Assoc,Mean_Conn_Order_Rnd_Assoc,WSR_EqPow,Mean_Conn_Order_Eq_Pow,...
    WSR_main_FP,Mean_Conn_Order_FP,WSR_Main_ADMM,Mean_Conn_Order_ADMM,WSR_main_FP_Sing_Conn,Mean_Conn_Order_FP_Sing_Conn,C3_Violation] = ...
    InnerFunc_BS_Varies(BS_Vec,Assoc_params,S,N_Orig,Ther_Noise,OptParams,BiSecParam,...
        SystemParams,seed);
    
    WSR_main_FP_Rnd_Assoc_MNT(mnt,:,:) = WSR_main_FP_Rnd_Assoc;
    Mean_Conn_Order_Rnd_Assoc_MNT(mnt,:,:) = Mean_Conn_Order_Rnd_Assoc;
    %Users_Rnd(i,j,:) = mean(Users_Rate_FP_Rnd_Assoc,1);
    WSR_EqPow_MNT(mnt,:,:) = WSR_EqPow;
    Mean_Conn_Order_Eq_Pow_MNT(mnt,:,:) = Mean_Conn_Order_Eq_Pow;
    %Users_EqPow(i,j,:) = mean(User_Rate_EqPow,1);
    WSR_main_FP_MNT(mnt,:,:) = WSR_main_FP;
    Mean_Conn_Order_FP_MNT(mnt,:,:) = Mean_Conn_Order_FP;
    %Users_FP(i,j,:) = mean(Users_Rate_FP,1);
    WSR_Main_ADMM_MNT(mnt,:,:) = WSR_Main_ADMM;
    Mean_Conn_Order_ADMM_MNT(mnt,:,:) = Mean_Conn_Order_ADMM;
    %Users_ADMM(i,j,:) = mean(Users_Rate_ADMM,1);
    %clear Users_Rate_FP_Rnd_Assoc User_Rate_EqPow Users_Rate_FP Users_Rate_ADMM
    WSR_Main_Sing_Conn_MNT(mnt,:,:) = WSR_main_FP_Sing_Conn;
    Mean_Conn_Order_Sing_Conn_MNT(mnt,:,:) = Mean_Conn_Order_FP_Sing_Conn;
    C3_Violation_MNT(mnt) = C3_Violation;
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
C3_Violation_AL = squeeze(mean(C3_Violation_MNT,1));




load('BS_varies-Aug29-CCV2-N_12-q_002-Rmax_30-GamL_1-IdlRndBS-MNT_1440.mat','WSR_main_FP_Rnd_Assoc_AL','Mean_Conn_Order_Rnd_Assoc_AL','WSR_EqPow_AL','Mean_Conn_Order_Eq_Pow_AL'...
    ,'WSR_main_FP_AL','Mean_Conn_Order_FP_AL','WSR_Main_ADMM_AL','Mean_Conn_Order_ADMM_AL','WSR_Main_Sing_Conn_AL','Mean_Conn_Order_Sing_Conn_AL')

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
plot(BS_Vec,Gb_factor*WSR_main_FP_AL(1,:),'-*','Color',Color1,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
hold on
%plot(BS_Vec,Gb_factor*WSR_main_FP_AL(2,:),'-*','Color',Color1,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
%hold on
plot(BS_Vec,Gb_factor*WSR_Main_ADMM_AL(1,:),'-*','Color',Color2,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
hold on
%plot(BS_Vec,Gb_factor*WSR_Main_ADMM_AL(2,:),'-*','Color',Color2,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
%hold on
plot(BS_Vec,Gb_factor*WSR_EqPow_AL(1,:),'-*','Color',Color4,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
hold on
%plot(BS_Vec,Gb_factor*WSR_EqPow_AL(2,:),'-*','Color',Color4,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
% hold on
plot(BS_Vec,Gb_factor*WSR_main_FP_Rnd_Assoc_AL(1,:),'-*','Color',Color6,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
hold on
% plot(BS_Vec,Gb_factor*WSR_main_FP_Rnd_Assoc_AL(:,2),'-*','Color',Color3,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
%hold on
plot(BS_Vec,Gb_factor*WSR_Main_Sing_Conn_AL(1,:),'-s','Color',Color3,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
%hold on
%plot(BS_Vec,Gb_factor*WSR_Main_Sing_Conn_AL(2,:),'-*','Color',Color3,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
grid on
xticks([2:2:25])
yticks([100:200:1500])
xlim([2,24])
legend({'Proposed-Algo1','Proposed-Algo2',...
    'Eq-Power','Algo1-Rnd-Assoc',...
    'Proposed-Algo1-SC'},'FontSize',7)
ylabel('System Sum Rate [Gbps]', 'FontSize', 10)
xlabel('Number of BSs, B', 'FontSize', 10)



figure(2)
box on
plot(BS_Vec,Mean_Conn_Order_FP_AL(1,:),'-*','Color',Color1,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
hold on
%plot(BS_Vec,Mean_Conn_Order_FP_AL(2,:),'-*','Color',Color1,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
%hold on
plot(BS_Vec,Mean_Conn_Order_ADMM_AL(1,:),'-*','Color',Color2,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
hold on
%plot(BS_Vec,Mean_Conn_Order_ADMM_AL(2,:),'-*','Color',Color2,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
%hold on
plot(BS_Vec,Mean_Conn_Order_Eq_Pow_AL(1,:),'-*','Color',Color4,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
hold on
%plot(BS_Vec,Mean_Conn_Order_Eq_Pow_AL(2,:),'-*','Color',Color4,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
% hold on
plot(BS_Vec,Mean_Conn_Order_Rnd_Assoc_AL(1,:),'-*','Color',Color6,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
hold on
% plot(BS_Vec,Mean_Conn_Order_Rnd_Assoc_AL(:,2),'-*','Color',Color3,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
%hold on
plot(BS_Vec,Mean_Conn_Order_Sing_Conn_AL(1,:),'-s','Color',Color3,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
%hold on
%plot(BS_Vec,Mean_Conn_Order_Sing_Conn_AL(2,:),'-*','Color',Color3,'MarkerIndices',1:2:length(BS_Vec),'LineWidth',0.8)
grid on
legend({'Proposed-Algo1','Proposed-Algo2',...
    'Eq-Power','Algo1-Rnd-Assoc',...
    'Proposed-Algo1-SC'},'FontSize',7)

ylabel('Average order of multi-connectivity', 'FontSize', 10)
xlabel('Number of BSs, B', 'FontSize', 10)
xticks([2:2:25])
%yticks([100:200:1500])
xlim([2,24])



function [Assoc_Out] = Association(Assoc_params,B,S,N,BL_mat,SE_mat)
Idle_BS_flg = Assoc_params.Idle_BS_flg;
Gam_L = Assoc_params.Gam_L;
Gam_U = Assoc_params.Gam_U;
% S = SystemParams.S;
% B = SystemParams.B;
% N = SystemParams.N;

Block_mat_1 = BL_mat;
Block_mat_1(Block_mat_1==0)=-100;
Block_mat = reshape(Block_mat_1,[B,1,N]);
Block_mat = repmat(Block_mat,[1,S,1]);
Block_mat = ones(B,S,N);
% Constraints matrices:
C = kron(ones(B,1)',eye(N*S,N*S));
L1 = kron(eye(S,S),ones(N,1)');
D = kron(eye(B), L1);
E = repmat(eye(N,N),[1,B*S]);
F = -E;
k1 = ones(N*S,1);
k2 = ones(B*S,1);
k3 = Gam_U*ones(N,1);
k4 = -Gam_L*ones(N,1);

Weight = Block_mat.*SE_mat; 
W_vec = reshape(permute(Weight,[3 2 1]), [B*S*N, 1]);
lb = zeros(size(W_vec)); 
ub = ones(size(W_vec));
if Idle_BS_flg == 0
    K_Cont_Ieq = [k1;k3;k4];
    T_Cont_Ieq = [C;E;F];
    options = optimoptions('linprog', 'Algorithm', 'dual-simplex',Display='none');
    [x, ~, ~, ~] = linprog(-W_vec, T_Cont_Ieq, K_Cont_Ieq, D, k2, lb, ub, [], options);
elseif Idle_BS_flg == 1
    K = [k1;k2;k3;k4];
    T = [C;D;E;F];
    options = optimoptions('linprog', 'Algorithm', 'dual-simplex',Display='none');
    [x, ~, ~, ~] = linprog(-W_vec, T, K, [], [], lb, ub, [], options);
end
if isempty(x)
    Assoc_Out = NaN;
else
    Assoc_Out = permute(reshape(x,[N,S,B]), [3,2,1]);
end

end



function [SE_mat] = SE_mat_Func(Power,B,S,N,Ther_Noise_mat,SystemParams,h,h_bar,h_tilde)
% S = SystemParams.S;
% B = SystemParams.B;
% N = SystemParams.N;
%Ther_Noise_mat = SystemParams.Ther_Noise_mat;
q = SystemParams.InvAntennaDirec;
% Channels:

P = Power;
H_bar = reshape(h_bar,[B,1,S,N]);
G_bar = repmat(H_bar,[1,B,1,1]);
E = G_bar .* reshape(P,[B,1,S,1]);
r_mask = reshape(eye(B,B),[B,B,1,1]);
i_mask = ones(B,B,S,N) - r_mask;
I = squeeze(sum(E.*i_mask,1));
Cum_Interf_mat = q * I;
% The total sum rate function:
Num_Des_Power = reshape(Power,[B,S,1]).*h; % Desired power for all users-BS-frequencies
Absorb_Noise = reshape(Power,[B,S,1]).*h_tilde; % Absorption noise for all users-BS-frequencies
Gamma_vec = (Num_Des_Power)./(Cum_Interf_mat + Absorb_Noise + Ther_Noise_mat);
SE_mat = log2(exp(1))*log(1+Gamma_vec); 
end


function [CDF, axis] = MyCDF(DataIn, NumberOfPoints)
 % Note: the input of this function may be a matrix or vector.
 % If the input is a matrix it first convert it into a vector.
 
Data = reshape(DataIn,[numel(DataIn),1]);
axisT = (linspace(min(Data),max(Data),NumberOfPoints)).';
NumberOfDataComponents = numel(DataIn);
     
count = zeros(NumberOfPoints,1);
for n = 1:NumberOfPoints
    count(n,1) = sum(Data <= axisT(n));
end
CDFT = count/NumberOfDataComponents;

axis = axisT(find(CDFT>0 & CDFT<1));  %Remove redundancies
CDF = CDFT(find(CDFT>0 & CDFT<1));  %Remove redundancies
end

