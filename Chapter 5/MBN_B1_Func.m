% MBN Function

function [SumRate_HO_MBN,Rate_HO_MBN,SumRate_HO_THz,SumRate_HO_RF,SumRate_MBN,Rate_MBN,Infeas_flg,Tmax,THz_SINR,RF_SINR,alpha_n,beta_n,W_mat,U_mat] = ...
    MBN_B1_Func(General_Struct,Traj_data,eps_conv,Pmax_THz,H_mat,H_tilde_mat,B_T,K,M,Rth_vec,THz_BW,Clus_Siz_THz,THz_User_Dist,RF_User_Dist,...
    Pmax_RF,G_mat,M_hat,B_R,RF_BW_Norm,Clus_Siz_RF,Clus_Siz_MBN,RF_BW,BL_mat,eta_THz,eta_RF,alpha_n_1,beta_n_1)


Infeas_flg = 0;

THz_Hyb = General_Struct.THz_Hyb;
RF_Hyb = General_Struct.RF_Hyb;


Ther_Noise_THz_Orig = General_Struct.Ther_Noise_THz_Orig;
Norm_fac_THz = Ther_Noise_THz_Orig^(0.5);
Ther_Noise_THz = 1; % Normalized
Ther_Noise_RF_Orig = General_Struct.Ther_Noise_RF_Orig;
Norm_fac_RF = Ther_Noise_RF_Orig^(0.5);
Ther_Noise_RF = 1; % Normalized



Infeas_flg = 0;
if THz_Hyb == 1
    M = K;
end

if RF_Hyb == 1
    M_hat = K;
end


BL_mat = BL_mat *1;
BL_mat_resh = repmat(reshape(BL_mat,[1,B_T,K]),[M,1,1]);
H_mat = H_mat .* BL_mat_resh;
H_tilde_mat = H_tilde_mat .* BL_mat_resh;

THz_User_Dist_Eff = (~BL_mat*1000) + THz_User_Dist;

if Clus_Siz_MBN > 1
    [Amat_THz] = mink(THz_User_Dist_Eff,Clus_Siz_THz,1);
    alpha_n_Heur = ismember(THz_User_Dist_Eff,Amat_THz);
    
    [Amat_RF] = mink(RF_User_Dist,Clus_Siz_RF,1);
    beta_n_Heur = ismember(RF_User_Dist,Amat_RF);
elseif Clus_Siz_MBN == 1
    All_Dist = [THz_User_Dist_Eff;RF_User_Dist];
    Amat_all = mink(All_Dist,Clus_Siz_MBN,1);
    All_Assoc = ismember(All_Dist,Amat_all);
    alpha_n_Heur = All_Assoc(1:B_T,:);
    beta_n_Heur = All_Assoc(B_T+1:end,:);
end




alpha_n = alpha_n_Heur;
beta_n = beta_n_Heur;


T = 30;
OptEps = eps_conv;
Obj = zeros(1,T);


MRT_Init_We = General_Struct.MRT_Init_We;
W_mat_t = MRT_Init_We * Norm_fac_THz .* conj(H_mat);
U_mat_t = MRT_Init_We * Norm_fac_RF .* conj(G_mat);



% Sum Rate Maximization Problem:
for t = 2 : T

    % THz auxiliary variables:
    Num_Aux_THz = squeeze(sum(sum(H_mat.*W_mat_t)))';
    A_t = Num_Aux_THz ./ (Interf_calc(H_mat,W_mat_t) + THz_Molec_Noise(H_tilde_mat,W_mat_t) + Ther_Noise_THz);


    % RF auxiliary variables:
    Num_Aux_RF = squeeze(sum(sum(G_mat.*U_mat_t)))';
    B_t = Num_Aux_RF ./ (Interf_calc(G_mat,U_mat_t) + Ther_Noise_RF);
    t;
    cvx_begin
    cvx_quiet(true)
    % cvx_precision high

    variable W_mat(M,B_T,K) complex
    variable p_THz_mat(B_T,K) nonnegative
    variable U_mat(M_hat,B_R,K) complex
    variable p_RF_mat(B_R,K) nonnegative

    [F_QoS_THz] = Obj_THz(W_mat,H_mat,H_tilde_mat,A_t,Ther_Noise_THz);
    
    [F_QoS_RF] = Obj_RF(U_mat,G_mat,B_t,Ther_Noise_RF,RF_BW_Norm);


     maximize sum( log(F_QoS_THz) + log(F_QoS_RF) )

    subject to



            F_QoS_THz + F_QoS_RF >= Rth_vec/THz_BW;
        

        %  THz:
            S = alpha_n - p_THz_mat;
            S = reshape(S,[1,B_T,K]);
            S_p = [2*W_mat;S];
            squeeze(norms(S_p,2,1)) <= alpha_n + p_THz_mat;
            sum(p_THz_mat') <= Pmax_THz;
            p_THz_mat <= 1 .* alpha_n * Pmax_THz;


        % RF:
            Y = beta_n - p_RF_mat;
            Y = reshape(Y,[1,B_R,K]);
            Y_p = [2*U_mat;Y];
            squeeze(norms(Y_p,2,1)) <= beta_n + p_RF_mat;
            sum(p_RF_mat') <= Pmax_RF;
            p_RF_mat <= 1 .* beta_n * Pmax_RF;

            

    cvx_end


    if (cvx_status(1) == 'I') || (cvx_status(1) == 'F') || (cvx_status(1) == 'E')
        Infeas_flg = 1;
        SumRate_MBN =0;
        Rate_MBN = zeros(1,K);
        SumRate_THz_MBN = 0;
        Rate_THz_User_MBN = zeros(1,K);
        THz_SINR = zeros(1,K);
        SumRate_RF_MBN = 0;
        Rate_RF_User_MBN = zeros(1,K);
        RF_SINR = zeros(1,K);
        Binary_Pen = 1000;  
        W_mat = zeros(M,B_T,K);
        U_mat = zeros(M_hat,B_R,K);
        SumRate_HO_MBN = 0;
        SumRate_HO_THz = 0;
        SumRate_HO_RF = 0;
        Rate_HO_MBN = zeros(1,K);
        alpha_n = ones(B_T,K);
        beta_n = ones(B_R,K);
        Tmax = t;

        break
    end

    Obj(t) = cvx_optval;
    % % t

    if abs(Obj(t)-Obj(t-1)) <= eps_conv || (Obj(t)-Obj(t-1) < 0) || t == T
        W_mat = W_mat_t;
        U_mat = U_mat_t;
        Tmax = t;
       break
    end

    W_mat_t = W_mat;

    U_mat_t = U_mat;

end


if Infeas_flg == 0

    [SumRate_HO_THz,SumRate_THz,Rate_HO_THz,Rate_THz,THz_SINR] = THz_Rate_func(General_Struct,Traj_data,H_mat,H_tilde_mat,W_mat_t,eta_THz,BL_mat,Ther_Noise_THz,THz_BW,THz_Hyb,alpha_n_1,alpha_n);

    [SumRate_HO_RF,SumRate_RF,Rate_HO_RF,Rate_RF,RF_SINR] = RF_Rate_func(General_Struct,Traj_data,G_mat,U_mat_t,eta_RF,Ther_Noise_RF,RF_BW,RF_Hyb,beta_n_1,beta_n);


    SumRate_HO_MBN = SumRate_HO_THz + SumRate_HO_RF;
    Rate_HO_MBN = Rate_HO_THz + Rate_HO_RF;

    if any(Rate_HO_MBN < Rth_vec)
        SumRate_HO_MBN = 0;
    end


    SumRate_MBN = SumRate_THz + SumRate_RF;
    Rate_MBN = Rate_THz + Rate_RF;


end



end



function [F_QoS_THz] = Obj_THz(W_mat,H_mat,H_tilde_mat,A_t,Ther_Noise_THz)

THz_SINR = 2*real(conj(A_t).* squeeze(sum(sum(H_mat.*W_mat)))') ...
    - abs(A_t).^2 .*(Interf_calc(H_mat,W_mat) + THz_Molec_Noise(H_tilde_mat,W_mat) + Ther_Noise_THz);


F_QoS_THz = log2(exp(1))* (log(1+THz_SINR));

end



function [F_QoS_RF] = Obj_RF(U_mat,G_mat,B_t,Ther_Noise_RF,RF_BW_Norm)

RF_SINR = 2*real(conj(B_t).* squeeze(sum(sum(G_mat.*U_mat)))') - abs(B_t).^2 .*(Interf_calc(G_mat,U_mat) + Ther_Noise_RF);

F_QoS_RF = RF_BW_Norm * log2(exp(1))* (log(1+RF_SINR));

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


