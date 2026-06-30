function [PL_RF,arr_Resp_RF,G_orig,Fron_Dist_RF] = RF_Channel(x,v,K,B_R,M_hat,RF_User_Dist,AoD_RF_Users_Angle,RF_Rayleigh_mat,Traj_data,General_Struct,NFC_flg)


% RF:
c = General_Struct.c;
f_RF = General_Struct.f_RF;
Ther_Noise_RF_Orig = General_Struct.Ther_Noise_RF_Orig;

RF_PL_expon = General_Struct.RF_PL_expon;

Rician_fac = General_Struct.RicianFac;


tau_eps_GT = Traj_data.tau_eps_GT;
N_GT = Traj_data.N_GT;


% RF Antenna Gain
G_Tx_RF = General_Struct.G_Tx_RF; % [dB]
G_Rx_RF = General_Struct.G_Rx_RF; % [dB]

%%%%%%%%%%% RF path loss
PL_RF = (c./(4*pi*f_RF)) .* (G_Tx_RF * G_Rx_RF)^0.5 .* (RF_User_Dist).^(-RF_PL_expon/2);

Norm_fac = (Ther_Noise_RF_Orig)^(0.5);

PL_RF = (PL_RF)./Norm_fac;


Ant_idx = (0:M_hat-1)';
Ant_idx = repmat(Ant_idx,[1,B_R,K,N_GT]);

Ant_idx = 0.5 * (2*Ant_idx - M_hat +1);

Wavelen_RF = c/f_RF;
d = x * Wavelen_RF; % Antenna spacing

Fron_Dist_RF = 2*(M_hat*d)^2 / Wavelen_RF;

RF_User_Distance = repmat(reshape(RF_User_Dist,[1,B_R,K,N_GT]),[M_hat,1,1,1]);
AoD_angle = repmat(reshape(AoD_RF_Users_Angle,[1,B_R,K,N_GT]),[M_hat,1,1,1]);


if NFC_flg == 1
    % NFC distances ( distances from each antenna elements)
    r_m = (RF_User_Distance.^2 + Ant_idx.^2*d^2 - 2*RF_User_Distance.*Ant_idx.*d.*cos(AoD_angle)).^0.5;
    v_m = ((RF_User_Distance.*cos(AoD_angle) - Ant_idx.*d) ./ r_m) * v;

    Phase_NFC_Dist = r_m - RF_User_Distance;
    Phase_NFC_Doppler = tau_eps_GT*v_m;

    arr_Resp_RF = exp((-1j*2*pi/Wavelen_RF) * (Phase_NFC_Dist - Phase_NFC_Doppler) );
else
    Phase_FFC_Dist =  - d*Ant_idx.*cos(AoD_angle);
    
    Phase_FFC_Doppler = tau_eps_GT*v*cos(AoD_angle);
    
    arr_Resp_RF = exp((-1j*2*pi/Wavelen_RF) * (Phase_FFC_Dist - Phase_FFC_Doppler) );
end


Ch_Resp_RF = sqrt(Rician_fac/(Rician_fac+1)).* arr_Resp_RF + sqrt(1/(Rician_fac+1)).* RF_Rayleigh_mat;
G_orig = reshape(PL_RF*Norm_fac,[1,B_R,K,N_GT]) .* Ch_Resp_RF;


end
