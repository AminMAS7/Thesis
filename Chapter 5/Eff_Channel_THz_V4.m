function [Eff_H_mat,Eff_H_tilde_mat,Eff_Arr_Resp_THz] = Eff_Channel_THz_V4(F_Analog,H_mat,PL_THz,PL_tilde_THz,B_T,K,M,General_Struct,Traj_data)

Ther_Noise_THz_Orig = General_Struct.Ther_Noise_THz_Orig;

Full_Conn_THz = ~General_Struct.Part_Conn_THz;
N_GT = Traj_data.N_GT;


Channel_phasor = exp(1j * angle(H_mat));

Channel_phasor_per = permute(Channel_phasor,[1,3,2,4]);

F_Analog_per = permute(F_Analog,[1,3,2,4]);


Eff_Arr_Resp_per = pagemtimes(pagectranspose(F_Analog_per),Channel_phasor_per);
Eff_Arr_Resp_THz = permute(Eff_Arr_Resp_per,[1,3,2,4]);


Eff_H_mat = reshape(PL_THz,[1,B_T,K,N_GT]) .* Eff_Arr_Resp_THz;

Eff_H_tilde_mat = reshape(PL_tilde_THz,[1,B_T,K,N_GT]) .* Eff_Arr_Resp_THz;


end