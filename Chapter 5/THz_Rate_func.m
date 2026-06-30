function [SumRate_HO_THz,SumRate_THz,Rate_HO_THz,Rate_THz,THz_SINR] = THz_Rate_func(General_Struct,Traj_data,H_mat,H_tilde_mat,W_mat,eta_THz,BL_mat,Ther_Noise_THz,THz_BW,THz_Hyb,alpha_n_1,alpha_n)


[M, B_T, K] = size(H_mat);



Ther_Noise_THz_Orig = General_Struct.Ther_Noise_THz_Orig;
Norm_fac_THz = Ther_Noise_THz_Orig^(0.5);

H_mat = H_mat * Norm_fac_THz;
H_tilde_mat = H_tilde_mat * Norm_fac_THz;
Ther_Noise_THz = Ther_Noise_THz_Orig;

if THz_Hyb == 1
    M = K;
end



BL_mat_resh = repmat(reshape(BL_mat,[1,B_T,K]),[M,1,1]);
H_mat = H_mat .* BL_mat_resh;
H_tilde_mat = H_tilde_mat .* BL_mat_resh;


alpha_n = alpha_n > 0.5;
alpha_n_1 = alpha_n_1 > 0.5;

W_mat = repmat(reshape(alpha_n,[1,B_T,K]),[M,1,1]) .* W_mat;


Num = squeeze(abs(sum(sum(H_mat.*W_mat))).^2)';
[I_THz] = Interf_calc(H_mat,W_mat);
THz_SINR = Num./(I_THz + THz_Molec_Noise(H_tilde_mat,W_mat) + Ther_Noise_THz);


Rate_THz = log2(exp(1))* THz_BW * log(1+THz_SINR);

SumRate_THz = sum(Rate_THz);

THz_HO_Cost = (1 - eta_THz * sum((1-alpha_n_1).*alpha_n));

Rate_HO_THz =  max( THz_HO_Cost .* Rate_THz, 0); 

SumRate_HO_THz = sum(Rate_HO_THz);



end

function [Interf] = Interf_calc(H,W)

[M, B_T, K] = size(H);

H_sq = reshape(H,[M,B_T,1,K]);
H_sq = repmat(H_sq,[1,1,K,1]);
W = reshape(W,[M,B_T,K,1]);
i_K_mask = ones(M,B_T,K,K) - reshape(eye(K,K),[1,1,K,K]);
E_1 = H_sq .* repmat(W,[1,1,1,K]) .* i_K_mask;
E_2 = squeeze(sum(E_1));
E_3 = squeeze(sum(E_2));
E = squeeze(pow_abs(E_3,2));
Interf = sum(E);


end


function [Molec_Noise] = THz_Molec_Noise(H,W)
[M, B_T, K] = size(H);

H_sq = reshape(H,[M,B_T,1,K]);
H_sq = repmat(H_sq,[1,1,K,1]);
W = reshape(W,[M,B_T,K,1]);
i_K_mask = ones(M,B_T,K,K);
E_1 = H_sq .* repmat(W,[1,1,1,K]) .* i_K_mask;
E2 = squeeze(sum(E_1));
E_3 = squeeze(sum(E2));
E = squeeze(pow_abs(E_3,2));
Molec_Noise = squeeze(sum(E));

end