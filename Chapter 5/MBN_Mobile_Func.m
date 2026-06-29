% MBN Function

function [F_HO_AveCh_err,F_HO_AveRate_err,F_HO_AveCh_GT,F_HO_AveRate_GT,NumHO_MBN,W_MAT,U_MAT,alpha_n_mat,beta_n_mat,Tmax] = MBN_Mobile_Func(Traj_data,K,B_R,M_hat,Channel_info,...
        H_min,General_Struct,B_T,M,...
    eps_conv,Pmax_THz,Rth_vec,THz_BW,...
    Clus_Siz_THz,Clus_Siz_RF,Clus_Siz_MBN,THz_Hyb,Func_Type,eta_BL,...
    Pmax_RF,RF_BW_Norm,RF_Hyb,RF_BW,Max_Min_Init,eta_THz,eta_RF,HO_Pen_1,HO_Pen_2,Ch_err,Users_Init_x,THz_BS_Loc,RF_BS_Loc,seed)

v = Traj_data.v;

N = Traj_data.N;

tau_eps_GT = Traj_data.tau_eps_GT;
N_GT = Traj_data.N_GT;
N_RF_Coh = Traj_data.N_RF_Coh;

N_Hori = Traj_data.N_Hori;

Ther_Noise_THz_Orig = General_Struct.Ther_Noise_THz_Orig;
Norm_fac_THz = Ther_Noise_THz_Orig^(0.5);
Ther_Noise_THz = 1; % Normalized
Ther_Noise_RF_Orig = General_Struct.Ther_Noise_RF_Orig;
Norm_fac_RF = Ther_Noise_RF_Orig^(0.5);
Ther_Noise_RF = 1; % Normalized



rng(seed)

RF_Rayleigh_MAT_TEMP = 1/sqrt(2) * (randn(M_hat,B_R,K,ceil(N_GT/N_RF_Coh),N) + 1j*randn(M_hat,B_R,K,ceil(N_GT/N_RF_Coh),N));
RF_Ray_idx = repelem(1:ceil(N_GT/N_RF_Coh), N_RF_Coh);
RF_Rayleigh_mat_TEMP = RF_Rayleigh_MAT_TEMP(:,:,:,RF_Ray_idx,:);


THz_Ch_err =  sqrt(Ch_err^2/2) * (randn(M,B_T,K,N) + 1j*randn(M,B_T,K,N));

RF_Ch_err =  sqrt(Ch_err^2/2) * (randn(M_hat,B_R,K,N) + 1j*randn(M_hat,B_R,K,N));

THz_User_Dist_err = sqrt(Ch_err/2) * randn(B_T,K);

RF_User_Dist_err = sqrt(Ch_err/2) * randn(B_R,K);


alpha_n_mat = zeros(B_T,K,N+1);
beta_n_mat = zeros(B_R,K,N+1);

W_MAT = zeros(K,B_T,K,N);
U_MAT = zeros(K,B_R,K,N);

for n = 1 : N
    UA_flg = 0;

    if mod(n,N_Hori) == 0 || n == 1
        UA_flg = 1;
    end


    THz_Ch_err_n = THz_Ch_err(:,:,:,n);
    RF_Ch_err_n = RF_Ch_err(:,:,:,n);

    RF_Rayleigh_mat_GT = RF_Rayleigh_mat_TEMP(:,:,:,:,n);

    [Users_x_mat_GT,Users_Init_y] = Users_pos(K,THz_BS_Loc,N_GT,tau_eps_GT,v,Users_Init_x,Channel_info);
    Users_Init_x = Users_x_mat_GT(:,end);

    [THz_User_Dist_GT,AoD_THz_Users_Angle_GT] = BS_Users_Dis_ang(Users_x_mat_GT,Users_Init_y,THz_BS_Loc,H_min);

    Blockage_Prob_MBN = exp(-THz_User_Dist_GT(:,:,1)*eta_BL);
    rng(seed*n)
    BL_mat = rand(size(THz_User_Dist_GT(:,:,1))) < Blockage_Prob_MBN;

    [PL_THz_GT,PL_tilde_THz_GT,~,H_GT,H_tilde_GT,~] = THz_Channel(v,K,B_T,M,THz_User_Dist_GT,AoD_THz_Users_Angle_GT,Traj_data,General_Struct);
    

    H_err =  H_GT -  (vecnorm(H_GT,2,1)) .* THz_Ch_err_n;


    F_Analog = exp(1j*angle(H_err))/sqrt(M);
    [Eff_H_mat_err_GT,Eff_H_tilde_mat_err_GT,~] = Eff_Channel_THz_V4(F_Analog,H_GT,PL_THz_GT,PL_tilde_THz_GT,B_T,K,M,General_Struct,Traj_data,0);
    [Eff_H_mat_err,Eff_H_tilde_mat_err,~] = Eff_Channel_THz_V4(F_Analog,H_err,PL_THz_GT,PL_tilde_THz_GT,B_T,K,M,General_Struct,Traj_data,0);
    
    

    Ave_Eff_H_mat_err = Eff_H_mat_err(:,:,:,1);
    Ave_Eff_H_mat_tilde_err = Eff_H_tilde_mat_err(:,:,:,1);


    H_mat = Ave_Eff_H_mat_err;
    H_tilde_mat = Ave_Eff_H_mat_tilde_err;

    THz_User_Dist_B1 = THz_User_Dist_GT(:,:,1) + THz_User_Dist_err;

    [RF_User_Dist_GT,AoD_RF_Users_Angle_GT] = BS_Users_Dis_ang(Users_x_mat_GT,Users_Init_y,RF_BS_Loc,H_min);
    
    [PL_RF_GT,~,G_GT,~] = RF_Channel(v,K,B_R,M_hat,RF_User_Dist_GT,AoD_RF_Users_Angle_GT,RF_Rayleigh_mat_GT,Traj_data,General_Struct);

    G_err = G_GT -  (vecnorm(G_GT,2,1)) .* RF_Ch_err_n;


    Q_Analog = exp(1j*angle(G_err))/sqrt(M_hat);
    [Eff_G_mat_err_GT,~] = Eff_Channel_RF_V4(Q_Analog,G_GT,PL_RF_GT,B_R,M_hat,K,General_Struct,Traj_data,0);
    [Eff_G_mat_err,~] = Eff_Channel_RF_V4(Q_Analog,G_err,PL_RF_GT,B_R,M_hat,K,General_Struct,Traj_data,0);

    Ave_Eff_G_mat_err = Eff_G_mat_err(:,:,:,1);

    G_mat = Ave_Eff_G_mat_err;

    RF_User_Dist_B1 = RF_User_Dist_GT(:,:,1) + RF_User_Dist_err;

    if n == 1


        idx = 1;

        alpha_n_1 = ones(B_T,K);
        beta_n_1 = ones(B_R,K);

        alpha_n_mat(:,:,idx) = alpha_n_1;
        beta_n_mat(:,:,idx) = beta_n_1;

    end

    idx = idx + 1;


    switch Func_Type

        case 'WithCost'
            [SumRate_HO_MBN(n),Rate_HO_MBN(:,n),SumRate_HO_THz(n),SumRate_HO_RF(n),SumRate_MBN(n),Rate_MBN(:,n),~,Tmax(n),THz_SINR(:,n),RF_SINR(:,n),alpha_n,beta_n,W_mat_opt,U_mat_opt] = ...
            MBN_Opt_Mob_Func(General_Struct,Traj_data,eps_conv,Pmax_THz,H_mat,H_tilde_mat,B_T,K,M,Rth_vec,THz_BW,Clus_Siz_THz,UA_flg,...
            Pmax_RF,G_mat,M_hat,B_R,RF_BW_Norm,Clus_Siz_RF,Clus_Siz_MBN,RF_BW,Max_Min_Init,BL_mat,eta_THz,eta_RF,alpha_n_1,beta_n_1);
        case 'Pen-1'
            HO_Pen = HO_Pen_1;
            [SumRate_HO_MBN(n),Rate_HO_MBN(:,n),SumRate_HO_THz(n),SumRate_HO_RF(n),SumRate_MBN(n),Rate_MBN(:,n),~,Tmax(n),THz_SINR(:,n),RF_SINR(:,n),alpha_n,beta_n,W_mat_opt,U_mat_opt] = ...
            MBN_Opt_Mob_HO_P_Func(General_Struct,Traj_data,eps_conv,Pmax_THz,H_mat,H_tilde_mat,B_T,K,M,Rth_vec,THz_BW,Clus_Siz_THz,HO_Pen,UA_flg,...
            Pmax_RF,G_mat,M_hat,B_R,RF_BW_Norm,Clus_Siz_RF,Clus_Siz_MBN,RF_BW,Max_Min_Init,BL_mat,eta_THz,eta_RF,alpha_n_1,beta_n_1);
       case 'Pen-2'
            HO_Pen = HO_Pen_2;
            [SumRate_HO_MBN(n),Rate_HO_MBN(:,n),SumRate_HO_THz(n),SumRate_HO_RF(n),SumRate_MBN(n),Rate_MBN(:,n),~,Tmax(n),THz_SINR(:,n),RF_SINR(:,n),alpha_n,beta_n,W_mat_opt,U_mat_opt] = ...
            MBN_Opt_Mob_HO_P_Func(General_Struct,Traj_data,eps_conv,Pmax_THz,H_mat,H_tilde_mat,B_T,K,M,Rth_vec,THz_BW,Clus_Siz_THz,HO_Pen,UA_flg,...
            Pmax_RF,G_mat,M_hat,B_R,RF_BW_Norm,Clus_Siz_RF,Clus_Siz_MBN,RF_BW,Max_Min_Init,BL_mat,eta_THz,eta_RF,alpha_n_1,beta_n_1);
        case 'NoCost'
            [SumRate_HO_MBN(n),Rate_HO_MBN(:,n),SumRate_HO_THz(n),SumRate_HO_RF(n),SumRate_MBN(n),Rate_MBN(:,n),~,Tmax(n),THz_SINR(:,n),RF_SINR(:,n),alpha_n,beta_n,W_mat_opt,U_mat_opt] = ...
            MBN_Opt_Stat_Func(General_Struct,Traj_data,eps_conv,Pmax_THz,H_mat,H_tilde_mat,B_T,K,M,Rth_vec,THz_BW,Clus_Siz_THz,UA_flg,...
            Pmax_RF,G_mat,M_hat,B_R,RF_BW_Norm,Clus_Siz_RF,Clus_Siz_MBN,RF_BW,Max_Min_Init,BL_mat,eta_THz,eta_RF,alpha_n_1,beta_n_1);
        case 'MBN-B1'
            [SumRate_HO_MBN(n),Rate_HO_MBN(:,n),SumRate_HO_THz(n),SumRate_HO_RF(n),SumRate_MBN(n),Rate_MBN(:,n),~,Tmax(n),THz_SINR(:,n),RF_SINR(:,n),alpha_n,beta_n,W_mat_opt,U_mat_opt] = ...
            MBN_B1_Func(General_Struct,Traj_data,eps_conv,Pmax_THz,H_mat,H_tilde_mat,B_T,K,M,Rth_vec,THz_BW,Clus_Siz_THz,THz_User_Dist_B1,RF_User_Dist_B1,...
                Pmax_RF,G_mat,M_hat,B_R,RF_BW_Norm,Clus_Siz_RF,Clus_Siz_MBN,RF_BW,BL_mat,eta_THz,eta_RF,alpha_n_1,beta_n_1);
    end

    [AveSumRate_HO_THz_err(n),AveSumRate_THz_err(n),AveRate_HO_THz_err(:,n),AveRate_THz_err(:,n),AveTHz_SINR_err(:,:,n)] = THz_AveRate_func(General_Struct,Traj_data,Eff_H_mat_err,Eff_H_tilde_mat_err,W_mat_opt,eta_THz,BL_mat,Ther_Noise_THz,THz_BW,THz_Hyb,alpha_n_1,alpha_n);

    [AveSumRate_HO_RF_err(n),AveSumRate_RF_err(n),AveRate_HO_RF_err(:,n),AveRate_RF_err(:,n),AveRF_SINR_err(:,:,n)] = RF_AveRate_func(General_Struct,Traj_data,Eff_G_mat_err,U_mat_opt,eta_RF,Ther_Noise_RF,RF_BW,RF_Hyb,beta_n_1,beta_n);

    AveRate_HO_MBN_err_err(:,n) = AveRate_HO_THz_err(:,n) + AveRate_HO_RF_err(:,n);


    Rate_AvRate_HO_MBN_err = (AveRate_HO_MBN_err_err(:,n)).';

    AveSumRate_HO_MBN_err_err(n) = sum(Rate_AvRate_HO_MBN_err);

    if any(Rate_AvRate_HO_MBN_err < Rth_vec)
        AveSumRate_HO_MBN_err_err(n) = 0;
    end


    AveSumRate_MBN_err_err(n) = AveSumRate_THz_err(n) + AveSumRate_RF_err(n);
    AveRate_MBN_err_err(:,n) = AveRate_THz_err(:,n) + AveRate_RF_err(:,n);

    
    [SumRate_HO_THz_AveCh_err_GT(n),SumRate_THz_AveCh_err_GT(n),Rate_HO_THz_AveCh_err_GT(:,n),Rate_THz_AveCh_err_GT(:,n),THz_SINR_AveCh_err_GT(:,n)]  = THz_Rate_func(General_Struct,Traj_data,Eff_H_mat_err_GT(:,:,:,1),Eff_H_tilde_mat_err_GT(:,:,:,1),W_mat_opt,eta_THz,BL_mat,Ther_Noise_THz,THz_BW,THz_Hyb,alpha_n_1,alpha_n);
    [SumRate_HO_RF_AveCh_err_GT(n),SumRate_RF_AveCh_err_GT(n),Rate_HO_RF_AveCh_err_GT(:,n),Rate_RF_AveCh_err_GT(:,n),RF_SINR_AveCh_err_GT(:,n)]  = RF_Rate_func(General_Struct,Traj_data,Eff_G_mat_err_GT(:,:,:,1),U_mat_opt,eta_RF,Ther_Noise_RF,RF_BW,RF_Hyb,beta_n_1,beta_n);

    SumRate_HO_MBN_AveCh_err_GT(n) = SumRate_HO_THz_AveCh_err_GT(n) + SumRate_HO_RF_AveCh_err_GT(n);
    Rate_HO_MBN_AveCh_err_GT(:,n) = Rate_HO_THz_AveCh_err_GT(:,n) + Rate_HO_RF_AveCh_err_GT(:,n);

    Rate_HO_MBN_AveCh_GT = (Rate_HO_MBN_AveCh_err_GT(:,n)).';

    SumRate_HO_MBN_AveCh_GT(n) = sum(Rate_HO_MBN_AveCh_GT);
    
    if any(Rate_HO_MBN_AveCh_GT < Rth_vec)
        SumRate_HO_MBN_AveCh_GT(n) = 0;
    end


    SumRate_MBN_AveCh_err_GT(n) = SumRate_THz_AveCh_err_GT(n) + SumRate_RF_AveCh_err_GT(n);
    Rate_MBN_AveCh_err_GT(:,n) = Rate_THz_AveCh_err_GT(:,n) + Rate_RF_AveCh_err_GT(:,n);


    [AveSumRate_HO_THz_err_GT(n),AveSumRate_THz_err_GT(n),AveRate_HO_THz_err_GT(:,n),AveRate_THz_err_GT(:,n),AveTHz_SINR_err_GT(:,:,n)] = THz_AveRate_func(General_Struct,Traj_data,Eff_H_mat_err_GT,Eff_H_tilde_mat_err_GT,W_mat_opt,eta_THz,BL_mat,Ther_Noise_THz,THz_BW,THz_Hyb,alpha_n_1,alpha_n);

    [AveSumRate_HO_RF_err_GT(n),AveSumRate_RF_err_GT(n),AveRate_HO_RF_err_GT(:,n),AveRate_RF_err_GT(:,n),AveRF_SINR_err_GT(:,:,n)] = RF_AveRate_func(General_Struct,Traj_data,Eff_G_mat_err_GT,U_mat_opt,eta_RF,Ther_Noise_RF,RF_BW,RF_Hyb,beta_n_1,beta_n);

    AveSumRate_HO_MBN_err_GT(n) = AveSumRate_HO_THz_err_GT(n) + AveSumRate_HO_RF_err_GT(n);
    AveRate_HO_MBN_err_GT(:,n) = AveRate_HO_THz_err_GT(:,n) + AveRate_HO_RF_err_GT(:,n);

    AveSumRate_MBN_err_GT(n) = AveSumRate_THz_err_GT(n) + AveSumRate_RF_err_GT(n);
    AveRate_MBN_err_GT(:,n) = AveRate_THz_err_GT(:,n) + AveRate_RF_err_GT(:,n);


    THz_HO_Count_per_User(:,n) = (sum((1-alpha_n_1).*alpha_n));
    RF_HO_Count_per_User(:,n) = (sum((1-beta_n_1).*beta_n));

    HO_Count_per_User(:,n) = (sum((1-alpha_n_1).*alpha_n)) + (sum((1-beta_n_1).*beta_n));


    W_MAT(:,:,:,n) = W_mat_opt;
    U_MAT(:,:,:,n) = U_mat_opt;

    alpha_n_1 = alpha_n > 0.5;
    beta_n_1 = beta_n > 0.5;
    alpha_n_mat(:,:,idx) = alpha_n > 0.5;
    beta_n_mat(:,:,idx) = beta_n > 0.5;

end





alpha_n_mat = alpha_n_mat > 0.5;
beta_n_mat = beta_n_mat > 0.5;

Assos_n_mat = [alpha_n_mat;beta_n_mat];

j = 0;
for i = 2:N+1
    j = j + 1;
    NumHO_MBN(j) = sum(sum((1-Assos_n_mat(:,:,i-1)).*Assos_n_mat(:,:,i)));
end
NumHO_MBN = sum(NumHO_MBN);


F_HO_AveCh_err = mean(SumRate_HO_MBN);


F_HO_AveRate_err = mean(AveSumRate_HO_MBN_err_err);


F_HO_AveCh_GT = mean(SumRate_HO_MBN_AveCh_GT);


F_HO_AveRate_GT = mean(AveSumRate_HO_MBN_err_GT);


end



% function [alpha_n_Heur,beta_n_Heur] = Init_Assoc(BL_mat,THz_User_Dist,Clus_Siz_THz,RF_User_Dist,Clus_Siz_RF)
% 
% THz_User_Dist_Eff = (~BL_mat*1000) + THz_User_Dist;
% 
% [Amat_THz] = mink(THz_User_Dist_Eff,Clus_Siz_THz,1);
% alpha_n_Heur = ismember(THz_User_Dist_Eff,Amat_THz);
% 
% [Amat_RF] = mink(RF_User_Dist,Clus_Siz_RF,1);
% beta_n_Heur = ismember(RF_User_Dist,Amat_RF);
% 
% 
% end



