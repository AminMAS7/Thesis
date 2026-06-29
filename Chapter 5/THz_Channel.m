
function [PL_THz,PL_tilde_THz,arr_Resp_THz,H_orig,H_tilde_orig,Fron_Dist_THz] = THz_Channel(x,v,K,B_T,M,THz_User_Dist,AoD_THz_Users_Angle,Traj_data,General_Struct,NFC_flg)

c = General_Struct.c;
f_THz = General_Struct.f_THz;
Ther_Noise_THz_Orig = General_Struct.Ther_Noise_THz_Orig;
k_a = General_Struct.k_a;


tau_eps_GT = Traj_data.tau_eps_GT;
N_GT = Traj_data.N_GT;

%%%%%%%%%%% THz Antenna Gain
G_Tx = General_Struct.G_Tx_THz; % [dB]
G_Rx = General_Struct.G_Rx_THz; % [dB]


%%%%%%%%%%% THz path loss:
PL_THz = (c./(4*pi*f_THz*THz_User_Dist)) .* (exp(-k_a*THz_User_Dist)).^0.5 * (G_Tx * G_Rx)^0.5;
PL_tilde_THz = (c./(4*pi*f_THz*THz_User_Dist)) .* (1-exp(-k_a*THz_User_Dist)).^0.5 * (G_Tx * G_Rx)^0.5;

% Normalizing the path loss
Norm_fac = Ther_Noise_THz_Orig^(0.5);
PL_THz = PL_THz./Norm_fac;
PL_tilde_THz = PL_tilde_THz./Norm_fac;


% Effective array Resposne


Ant_idx = (0:M-1)';

Ant_idx = repmat(Ant_idx,[1,B_T,K,N_GT]);

Ant_idx = 0.5 * (2*Ant_idx - M +1);


Wavelen_THz = c/f_THz;
d = x * Wavelen_THz; % Antenna spacing

Fron_Dist_THz = 2*(M*d)^2 / Wavelen_THz;

THz_User_Distance = repmat(reshape(THz_User_Dist,[1,B_T,K,N_GT]),[M,1,1,1]);
AoD_angle = repmat(reshape(AoD_THz_Users_Angle,[1,B_T,K,N_GT]),[M,1,1,1]);




if NFC_flg == 1

% NFC distances ( distances from each antenna elements)
r_m = (THz_User_Distance.^2 + Ant_idx.^2*d^2 - 2*THz_User_Distance.*Ant_idx.*d.*cos(AoD_angle)).^0.5;

v_m = ((THz_User_Distance.*cos(AoD_angle) - Ant_idx.*d) ./ r_m) * v;


Phase_NFC_Dist = r_m - THz_User_Distance;
Phase_NFC_Doppler = tau_eps_GT*v_m;

arr_Resp_THz = exp((-1j*2*pi/Wavelen_THz) * (Phase_NFC_Dist  - Phase_NFC_Doppler));

else

Phase_FFC_Dist =  - d*Ant_idx.*cos(AoD_angle);

Phase_FFC_Doppler = tau_eps_GT*v*cos(AoD_angle);

arr_Resp_THz = exp((-1j*2*pi/Wavelen_THz) * (Phase_FFC_Dist - Phase_FFC_Doppler) );

end

H_orig = reshape(PL_THz*Norm_fac,[1,B_T,K,N_GT]) .* arr_Resp_THz;
H_tilde_orig = reshape(PL_tilde_THz*Norm_fac,[1,B_T,K,N_GT]) .* arr_Resp_THz;

end

