% MBN Function

function [SumRate_HO_THzOn,Rate_HO_THzOn,SumRate_THzOn,Rate_THzOn,Infeas_flg,Tmax,THz_SINR,alpha_n,W_mat] = ...
    THzOn_B1_Func(General_Struct,Traj_data,eps_conv,Pmax_THz,H_mat,H_tilde_mat,B_T,K,M,Rth_vec,THz_BW,Clus_Siz_MBN,THz_User_Dist,BL_mat,eta_THz,alpha_n_1)

Infeas_flg = 0;

THz_Hyb = General_Struct.THz_Hyb;


if THz_Hyb == 1
    M = K;
end



BL_mat = BL_mat *1;
BL_mat_resh = repmat(reshape(BL_mat,[1,B_T,K]),[M,1,1]);
H_mat = H_mat .* BL_mat_resh;
H_tilde_mat = H_tilde_mat .* BL_mat_resh;


    
THz_User_Dist_Eff = (~BL_mat*1000) + THz_User_Dist;
[Amat_THz] = mink(THz_User_Dist_Eff,Clus_Siz_MBN,1);
alpha_n_Heur = ismember(THz_User_Dist_Eff,Amat_THz);


alpha_n = alpha_n_Heur;


T = 30;
OptEps = 1e-2;
Obj = zeros(1,T);

Ther_Noise_THz_Orig = General_Struct.Ther_Noise_THz_Orig;
Norm_fac_THz = Ther_Noise_THz_Orig^(0.5);
Ther_Noise_THz = 1; % Normalized


MRT_Init_We = General_Struct.MRT_Init_We;
W_mat_t = MRT_Init_We * Norm_fac_THz .* conj(H_mat);


S_HO = General_Struct.S_HO;


Obj = zeros(1,T);


for t = 2 : T

    % THz auxiliary variables:
    Num_Aux_THz = squeeze(sum(sum(H_mat.*W_mat_t)))';
    A_t = Num_Aux_THz ./ (Interf_calc(H_mat,W_mat_t) + THz_Molec_Noise(H_tilde_mat,W_mat_t) + Ther_Noise_THz);

    t;
    cvx_begin
    cvx_quiet(true)
    % cvx_precision high

    variable W_mat(M,B_T,K) complex
    variable p_THz_mat(B_T,K) nonnegative


    [F_QoS_THz] = Obj_THz(W_mat,H_mat,H_tilde_mat,A_t,Ther_Noise_THz);
    

    maximize sum( log(F_QoS_THz) )
  

    subject to
    F_QoS_THz >= 5e-2;

        F_QoS_THz >= Rth_vec/THz_BW;

        %  THz:
            S = alpha_n - p_THz_mat;
            S = reshape(S,[1,B_T,K]);
            S_p = [2*W_mat;S];
            squeeze(norms(S_p,2,1)) <= alpha_n + p_THz_mat;
            sum(p_THz_mat') <= Pmax_THz;
            p_THz_mat <= 1 .* alpha_n * Pmax_THz;

    cvx_end



    if (cvx_status(1) == 'I') || (cvx_status(1) == 'F') || (cvx_status(1) == 'E')
        Infeas_flg = 1;
        SumRate_THzOn =0;
        Rate_THzOn = zeros(1,K);
        THz_SINR = zeros(1,K);
        Binary_Pen = 1000;  
        W_mat = zeros(M,B_T,K);
        SumRate_HO_THzOn = 0;
        SumRate_HO_THzOn = 0;
        Rate_HO_THzOn = zeros(1,K);
        alpha_n = ones(B_T,K);
        Tmax = t;

        break
    end

    Obj(t) = cvx_optval;

    if abs(Obj(t)-Obj(t-1)) <= eps_conv || (Obj(t)-Obj(t-1) < 0) || t == T
        W_mat = W_mat_t;
        Tmax = t;
       break
    end

    W_mat_t = W_mat;


end


if Infeas_flg == 0

    [SumRate_HO_THz,SumRate_THz,Rate_HO_THz,Rate_THz,THz_SINR] = THz_Rate_func(General_Struct,Traj_data,H_mat,H_tilde_mat,W_mat_t,eta_THz,BL_mat,Ther_Noise_THz,THz_BW,THz_Hyb,alpha_n_1,alpha_n);


    SumRate_HO_THzOn = SumRate_HO_THz;
    Rate_HO_THzOn = Rate_HO_THz ;

    Rate_HO_THzOn(Rate_HO_THzOn<Rth_vec) = 0;
    SumRate_HO_THzOn = sum(Rate_HO_THzOn);

    SumRate_THzOn = SumRate_THz;
    Rate_THzOn = Rate_THz;


end



end






function [F_QoS_THz] = Obj_THz(W_mat,H_mat,H_tilde_mat,A_t,Ther_Noise_THz)

THz_SINR = 2*real(conj(A_t).* squeeze(sum(sum(H_mat.*W_mat)))') ...
    - abs(A_t).^2 .*(Interf_calc(H_mat,W_mat) + THz_Molec_Noise(H_tilde_mat,W_mat) + Ther_Noise_THz);


F_QoS_THz = log2(exp(1))* (log(1+THz_SINR));

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



