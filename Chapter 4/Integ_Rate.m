% Author: Mohammad Amin Saeidi
% Date: August 2023
%
% Usage and Redistribution Notice:
% - This code is provided for [educational/research/personal] purposes only.
% - Any use, modification, redistribution, or other use of this code
%   requires explicit permission from the author.
% - To request permission, please send an email to [amin96a@yorku.ca or m.amin.saeidi96@gmail.com].

function [Output_Int] = Integ_Rate(D_ue_RF_THzbs,Params,fadeRand_Integ)

% Loading parameters: 
GT = Params.RF_Tx_gain;
GR = Params.RF_Rx_gain;
fcRF = Params.fcRF;
PR = Params.P_RF;
GTT = Params.THz_Tx_gain;
GRR = Params.THz_Rx_gain;
fcTHz = Params.fcTHz;
PT = Params.P_THz;
Thermal_Noise_RF = Params.Thermal_Noise_RF;
Thermal_Noise_THz = Params.Thermal_Noise_THz;

Wr = Params.RF_Bandwidth;
Wt = Params.THz_Bandwidth;
kf = Params.kf; % Molecular absorption coefficient at f = 1 THz
prob = Params.MissAll_Probability;
% Constants:
c=3*10^8; % Velocity of light
alpha = 4.01; % Path-loss exponent
aT = 2;% Path-loss exponent for THz transmission (Based on the mobility-aware paper, it is set to be 2)
gammaI=(c)^2*GT*GR/(16*pi^2*fcRF^2); % Function of constants for RF
gammaII=(c)^2*GTT*GRR/(16*pi^2*fcTHz^2); % Function of constants for THz

Small_dis = find(D_ue_RF_THzbs < 1);
D_ue_RF_THzbs(Small_dis) = 1;

% Obtaining the list of all SINR from all RF/THz BSs at the user of interset
% SINR and rate of RF mode:
RPrAllu = gammaI.*fadeRand_Integ.*PR.*D_ue_RF_THzbs.^(-alpha);  % Received power from all RF/THz BS in RF mode at the user of interest
DenomVec_RF = ones(length(D_ue_RF_THzbs),1).*sum(RPrAllu) - RPrAllu;
SINR_RF_vec = RPrAllu./(DenomVec_RF + Thermal_Noise_RF);
RF_rate_vec = Wr.*log2(1+SINR_RF_vec);

% SINR and rate of THz mode:

TPrAllu = (gammaII.*PT*exp(-kf.*D_ue_RF_THzbs))./(D_ue_RF_THzbs.^aT); % Received power from all Rf/THz BS in THz mode at the user of interest
NoiseTPrAllu = (gammaII.*PT*(1-exp(-kf.*D_ue_RF_THzbs)))./(D_ue_RF_THzbs.^aT); % Molecular absorption noise from all RF/THz BSs in THz mode
TCumPr = sum(TPrAllu); % Sum of the received power from RF/THz BS in THz mode at the user of interest
TCumNr = sum(NoiseTPrAllu); % Cumulative molecular absorption noise from all RF/THz BSs in THz mode
DenomVec_THz = prob.*(TCumPr - TPrAllu + TCumNr - NoiseTPrAllu) + NoiseTPrAllu + Thermal_Noise_THz;
SINR_THz_vec = TPrAllu ./ DenomVec_THz; 
THz_rate_vec = Wt*log2(1+SINR_THz_vec);

%%%%%%%%%% BS mode selection (Association) %%%%%%%%%%%%%%%%%%%%%%%

%%%% Association metric: Data Rate (Upper bound)
RFassoc = 0;
THassoc = 0;
[max_rate,max_rate_idx] = max([RF_rate_vec;THz_rate_vec]);

if max_rate_idx <= length(D_ue_RF_THzbs)
    RFassoc = RFassoc + 1;
    max_SE_Ass_rate = max_rate/Wr;
elseif max_rate_idx > length(D_ue_RF_THzbs)
    THassoc = THassoc + 1;
    max_SE_Ass_rate = max_rate/Wt;
end

%%%% Association metric: Test
Rate_vec = [RF_rate_vec;THz_rate_vec];
RF_Assc = Wr * RPrAllu;
THz_Assc = Wt * TPrAllu .* exp(-kf*D_ue_RF_THzbs);
[~,max_Assc_idx] = max([RF_Assc;THz_Assc]);
max_rate_T1 = Rate_vec(max_Assc_idx);

RFassoc_T1 = 0;
THassoc_T1 = 0;
if max_Assc_idx <= length(D_ue_RF_THzbs)
    RFassoc_T1 = RFassoc_T1 + 1;
    max_SE_Ass_T1 = max_rate_T1/Wr;
elseif max_Assc_idx > length(D_ue_RF_THzbs)
    THassoc_T1 = THassoc_T1 + 1;
    max_SE_Ass_T1 = max_rate_T1/Wt;
end

%%%% Association metric: SINR
RFassoc_SINR = 0;
THassoc_SINR= 0;
[max_SINR,max_SINR_idx] = max([SINR_RF_vec;SINR_THz_vec]);

if max_SINR_idx <= length(D_ue_RF_THzbs)
    max_rate_SINR = Wr.*log2(1+max_SINR);
    RFassoc_SINR = RFassoc_SINR + 1;
    max_SE_Ass_SINR = max_rate_SINR/Wr;
elseif max_SINR_idx > length(D_ue_RF_THzbs)
    max_rate_SINR = Wt*log2(1+max_SINR);
    THassoc_SINR = THassoc_SINR + 1;
    max_SE_Ass_SINR = max_rate_SINR/Wt;
end

%%%% Association metric: Recieved Power
RFassoc_Pow = 0;
THassoc_Pow= 0;
[~,max_Pow_idx] = max([RPrAllu;TPrAllu]);
max_rate_Pow = Rate_vec(max_Pow_idx);

if max_Pow_idx <= length(D_ue_RF_THzbs)
    RFassoc_Pow  = RFassoc_Pow  + 1;
    max_SE_Ass_Power = max_rate_Pow/Wr;
elseif max_Pow_idx > length(D_ue_RF_THzbs)
    THassoc_Pow  = THassoc_Pow + 1;
    max_SE_Ass_Power = max_rate_Pow/Wt;
end

Output_Int = struct();

Output_Int.max_rate_Upper = max_rate;
Output_Int.RFassoc = RFassoc;
Output_Int.THassoc = THassoc;
Output_Int.max_SE_Ass_rate = max_SE_Ass_rate;

Output_Int.max_rate_T1 = max_rate_T1;
Output_Int.RFassoc_T1 = RFassoc_T1;
Output_Int.THassoc_T1 = THassoc_T1;
Output_Int.max_SE_Ass_T1 = max_SE_Ass_T1;

Output_Int.max_rate_SINR = max_rate_SINR;
Output_Int.RFassoc_SINR = RFassoc_SINR;
Output_Int.THassoc_SINR = THassoc_SINR;
Output_Int.max_SE_Ass_SINR = max_SE_Ass_SINR;

Output_Int.max_rate_Pow = max_rate_Pow;
Output_Int.RFassoc_Pow = RFassoc_Pow;
Output_Int.THassoc_Pow = THassoc_Pow;
Output_Int.max_SE_Ass_Power = max_SE_Ass_Power;


end
    

 