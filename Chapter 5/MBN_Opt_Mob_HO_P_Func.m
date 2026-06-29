% MBN Function

function [SumRate_HO_MBN,Rate_HO_MBN,SumRate_HO_THz,SumRate_HO_RF,SumRate_MBN,Rate_MBN,Infeas_flg,Tmax,THz_SINR,RF_SINR,alpha_n,beta_n,W_mat,U_mat] = ...
    MBN_Opt_Mob_HO_P_Func(General_Struct,Traj_data,eps_conv,Pmax_THz,H_mat,H_tilde_mat,B_T,K,M,Rth_vec,THz_BW,Clus_Siz_THz,HO_Pen,UA_flg,...
    Pmax_RF,G_mat,M_hat,B_R,RF_BW_Norm,Clus_Siz_RF,Clus_Siz_MBN,RF_BW,Max_Min_Init,BL_mat,eta_THz,eta_RF,alpha_n_1,beta_n_1)

Infeas_flg = 0;

THz_Hyb = General_Struct.THz_Hyb;
RF_Hyb = General_Struct.RF_Hyb;


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


delta_1 = 1e3;
delta_2 = 1e3;


alpha_n_t = 0.5* ones(B_T,K);
beta_n_t = 0.5* ones(B_R,K);


OptEps = 1e-2;

if UA_flg == 1
    T = 50;
else
    T = General_Struct.T_NUA;
end
Obj = zeros(1,T);


Ther_Noise_THz_Orig = General_Struct.Ther_Noise_THz_Orig;
Norm_fac_THz = Ther_Noise_THz_Orig^(0.5);
Ther_Noise_THz = 1; % Normalized
Ther_Noise_RF_Orig = General_Struct.Ther_Noise_RF_Orig;
Norm_fac_RF = Ther_Noise_RF_Orig^(0.5);
Ther_Noise_RF = 1; % Normalized

MRT_Init_We = General_Struct.MRT_Init_We;
W_mat_t = MRT_Init_We * Norm_fac_THz .* conj(H_mat);
U_mat_t = MRT_Init_We * Norm_fac_RF .* conj(G_mat);


S_HO = General_Struct.S_HO;

MaxMin_flg = General_Struct.MaxMin_flg;

if Max_Min_Init == 1 && UA_flg == 1

for t = 2 : 10

    % THz auxiliary variables:
    Num_Aux_THz = squeeze(sum(sum(H_mat.*W_mat_t)))';
    A_t = Num_Aux_THz ./ (Interf_calc(H_mat,W_mat_t) + THz_Molec_Noise(H_tilde_mat,W_mat_t) + Ther_Noise_THz);


    % RF auxiliary variables:
    Num_Aux_RF = squeeze(sum(sum(G_mat.*U_mat_t)))';
    B_t = Num_Aux_RF ./ (Interf_calc(G_mat,U_mat_t) + Ther_Noise_RF);

    cvx_begin
    cvx_quiet(true)
    % cvx_precision high

    variable W_mat(M,B_T,K) complex
    variable U_mat(M_hat,B_R,K) complex
    if MaxMin_flg == 1
        variable Gamma_THz nonnegative
        variable Gamma_RF nonnegative
    end

    [F_QoS_THz] = Obj_THz(W_mat,H_mat,H_tilde_mat,A_t,Ther_Noise_THz);
    
    [F_QoS_RF] = Obj_RF(U_mat,G_mat,B_t,Ther_Noise_RF,RF_BW_Norm);

    if MaxMin_flg == 1
        maximize Gamma_THz + Gamma_RF
    else
        maximize sum(log(F_QoS_THz) + log(F_QoS_RF))
    end


    subject to
    
    if MaxMin_flg == 1
        F_QoS_THz >= Gamma_THz;
        F_QoS_RF >= Gamma_RF;
    end

        %  THz:
            W_mat_per = permute(W_mat,[1,3,2]);
            sum(sum(pow_abs(W_mat_per,2))) <= Pmax_THz;

        % RF:
           U_mat_per = permute(U_mat,[1,3,2]);
           sum(sum(pow_abs(U_mat_per,2))) <= Pmax_RF;

    cvx_end

    if (cvx_status(1) == 'I') || (cvx_status(1) == 'F') || (cvx_status(1) == 'E')

        W_mat_t = 1 * Norm_fac_THz .* conj(H_mat);
        U_mat_t = 1 * Norm_fac_RF .* conj(G_mat);
        break     
    end

    Obj(t) = cvx_optval;
    if abs(Obj(t)-Obj(t-1)) <= OptEps
        W_mat_t = W_mat;
        U_mat_t = U_mat;
       break
    end

    W_mat_t = W_mat;
    U_mat_t = U_mat;
end

end


%Sum Rate Maximization Problem:
Obj = zeros(1,T);

if UA_flg == 0
    alpha_n = alpha_n_1;
    beta_n = beta_n_1;
end



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

    if UA_flg == 1
        variable W_mat(M,B_T,K) complex
        variable alpha_n(B_T,K) nonnegative
        variable p_THz_mat(B_T,K) nonnegative
        variable U_mat(M_hat,B_R,K) complex
        variable beta_n(B_R,K) nonnegative
        variable p_RF_mat(B_R,K) nonnegative
    else
        variable W_mat(M,B_T,K) complex
        variable p_THz_mat(B_T,K) nonnegative
        variable U_mat(M_hat,B_R,K) complex
        variable p_RF_mat(B_R,K) nonnegative
    end

    [F_QoS_THz] = Obj_THz(W_mat,H_mat,H_tilde_mat,A_t,Ther_Noise_THz);
    
    [F_QoS_RF] = Obj_RF(U_mat,G_mat,B_t,Ther_Noise_RF,RF_BW_Norm);



    if UA_flg == 1
        maximize sum( log(F_QoS_THz) + log(F_QoS_RF) )...
            - HO_Pen * ( sum(sum((1-alpha_n_1).*alpha_n)) + sum(sum((1-beta_n_1).*beta_n)) )...
            - delta_1*(sum(sum(alpha_n.*(1-2*alpha_n_t) + alpha_n_t.^2))) - delta_2*(sum(sum(beta_n.*(1-2*beta_n_t) + beta_n_t.^2)))
    else
        maximize sum( log(F_QoS_THz) + log(F_QoS_RF) )
    end

    subject to



    if UA_flg == 1

        F_QoS_THz + F_QoS_RF >= Rth_vec/THz_BW;
       
        sum(alpha_n_1 .* alpha_n) + sum(beta_n_1 .* beta_n) >= S_HO;

        0 <= sum(alpha_n) <= Clus_Siz_THz;
        0 <= alpha_n <= 1;

        alpha_n <= BL_mat;

        0 <= sum(beta_n) <= Clus_Siz_RF;
        0 <= beta_n <= 1;

        sum(alpha_n) + sum(beta_n) <= Clus_Siz_MBN;


    else
        F_QoS_THz + F_QoS_RF >= Rth_vec/THz_BW;

    end


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
        alpha_n = alpha_n_1;
        beta_n = beta_n_1;
        Tmax = t;

        break
    end
    Binary_Pen(t) = (sum(sum(alpha_n.*(1-2*alpha_n_t) + alpha_n_t.^2))) + sum(sum(beta_n.*(1-2*beta_n_t) + beta_n_t.^2));


    Obj(t) = cvx_optval;


    if UA_flg == 1

        if abs(Obj(t)-Obj(t-1)) <= eps_conv || (Obj(t)-Obj(t-1) < 0 && Binary_Pen(t)< 1e-2) || t == T
            W_mat = W_mat_t;
            alpha_n = alpha_n_t;
            U_mat = U_mat_t;
            beta_n = beta_n_t;
            Tmax = t;
           break
        end
    else
        if abs(Obj(t)-Obj(t-1)) <= eps_conv || (Obj(t)-Obj(t-1)) < 0 || t == T
            W_mat = W_mat_t;
            U_mat = U_mat_t;
            Tmax = t;
           break
        end
    end


    if UA_flg == 1

        W_mat_t = W_mat;
        alpha_n_t = alpha_n;
    
        U_mat_t = U_mat;
        beta_n_t = beta_n;
    else
        W_mat_t = W_mat;
        U_mat_t = U_mat;
    end


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



