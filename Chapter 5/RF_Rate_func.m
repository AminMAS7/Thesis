function [SumRate_HO_RF,SumRate_RF,Rate_HO_RF,Rate_RF,RF_SINR] = RF_Rate_func(General_Struct,Traj_data,G_mat,U_mat,eta_RF,Ther_Noise_RF,RF_BW,RF_Hyb,beta_n_1,beta_n)

[M_hat, B_R, K] = size(G_mat);



Ther_Noise_RF_Orig = General_Struct.Ther_Noise_RF_Orig;
Norm_fac_RF = Ther_Noise_RF_Orig^(0.5);

G_mat = G_mat * Norm_fac_RF;
Ther_Noise_RF = Ther_Noise_RF_Orig;

if RF_Hyb == 1
    M_hat = K;
end

beta_n = beta_n > 0.5;
beta_n_1 = beta_n_1 > 0.5;

U_mat = repmat(reshape(beta_n,[1,B_R,K]),[M_hat,1,1]) .* U_mat;

Num = squeeze(abs(sum(sum(G_mat.*U_mat))).^2)';
[I_RF] = Interf_calc(G_mat,U_mat);
RF_SINR = Num./(I_RF + Ther_Noise_RF);

Overhead = 1 - sum(beta_n,1) * 0;

Rate_RF = log2(exp(1)) * RF_BW* log(1+RF_SINR);

SumRate_RF = sum(Rate_RF);

RF_HO_Cost = (1 - eta_RF * sum((1-beta_n_1).*beta_n));

Rate_HO_RF =  max( RF_HO_Cost .* Rate_RF,0); 

SumRate_HO_RF = sum(Rate_HO_RF);




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

