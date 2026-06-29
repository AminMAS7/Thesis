function [Eff_G_mat,Eff_Ch_Resp_RF] = Eff_Channel_RF_V4(Q_Analog,G_mat,PL_RF,B_R,M_hat,K,General_Struct,Traj_data)

Ther_Noise_RF_Orig = General_Struct.Ther_Noise_RF_Orig;

Full_Conn_RF = ~General_Struct.Part_Conn_RF;
N_GT = Traj_data.N_GT;


Channel_phasor = exp(1j * angle(G_mat));

Channel_phasor_per = permute(Channel_phasor,[1,3,2,4]);
Q_Analog_per = permute(Q_Analog,[1,3,2,4]);


Eff_Arr_Resp_per = pagemtimes(pagectranspose(Q_Analog_per),Channel_phasor_per);
Eff_Ch_Resp_RF = permute(Eff_Arr_Resp_per,[1,3,2,4]);


Eff_G_mat = reshape(PL_RF,[1,B_R,K,N_GT]) .* Eff_Ch_Resp_RF;


end