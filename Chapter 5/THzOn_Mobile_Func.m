% MBN Function

function [F_HO_AveCh_err,F_HO_AveRate_err,F_HO_AveCh_GT,F_HO_AveRate_GT,NumHO_THzOn,W_MAT,alpha_n_mat,Tmax] = THzOn_Mobile_Func(Traj_data,K,B_R,Channel_info,...
        H_min,General_Struct,B_T,M,...
    eps_conv,Pmax_THz,Rth_vec,THz_BW,...
    Clus_Siz_MBN,THz_Hyb,Func_Type,eta_BL,...
    Max_Min_Init,eta_THz,HO_Pen_1,HO_Pen_2,Ch_err,Users_Init_x,THz_BS_Loc,seed)

B_T = B_T + B_R;

v = Traj_data.v;

N = Traj_data.N;
tau_eps = Traj_data.tau_eps;

tau_eps_GT = Traj_data.tau_eps_GT;
N_GT = Traj_data.N_GT;

N_Hori = Traj_data.N_Hori;


Ther_Noise_THz_Orig = General_Struct.Ther_Noise_THz_Orig;
Norm_fac_THz = Ther_Noise_THz_Orig^(0.5);
Ther_Noise_THz = 1; % Normalized




rng(seed)

THz_Ch_err =  sqrt(Ch_err^2/2) * (randn(M,B_T,K,N) + 1j*randn(M,B_T,K,N));

THz_User_Dist_err = sqrt(Ch_err/2) * randn(B_T,K);


alpha_n_mat = zeros(B_T,K,N+1);

W_MAT = zeros(K,B_T,K,N);

for n = 1 : N
    n
    UA_flg = 0;

    if mod(n,N_Hori) == 0 || n == 1
        UA_flg = 1;
    end

    THz_Ch_err_n = THz_Ch_err(:,:,:,n);


    [Users_x_mat_GT,Users_Init_y] = Users_pos(K,THz_BS_Loc,N_GT,tau_eps_GT,v,Users_Init_x,Channel_info);
    Users_Init_x = Users_x_mat_GT(:,end);

    [THz_User_Dist_GT,AoD_THz_Users_Angle_GT] = BS_Users_Dis_ang(Users_x_mat_GT,Users_Init_y,THz_BS_Loc,H_min);

    Blockage_Prob_MBN = exp(-THz_User_Dist_GT(:,:,1)*eta_BL);
    rng(seed*n)
    BL_mat = rand(size(THz_User_Dist_GT(:,:,1))) < Blockage_Prob_MBN;

    [PL_THz_GT,PL_tilde_THz_GT,~,H_GT,H_tilde_GT,~] = THz_Channel(v,K,B_T,M,THz_User_Dist_GT,AoD_THz_Users_Angle_GT,Traj_data,General_Struct);
    

    H_err = 1 * H_GT - 1 * (vecnorm(H_GT,2,1)) .* THz_Ch_err_n;
    H_tilde_err = 1 * H_tilde_GT - 1 * (vecnorm(H_tilde_GT,2,1)) .* THz_Ch_err_n;


    F_Analog = exp(1j*angle(H_err))/sqrt(M);
    [Eff_H_mat_err_GT,Eff_H_tilde_mat_err_GT,~] = Eff_Channel_THz_V4(F_Analog,H_GT,PL_THz_GT,PL_tilde_THz_GT,B_T,K,M,General_Struct,Traj_data,0);
    [Eff_H_mat_err,Eff_H_tilde_mat_err,~] = Eff_Channel_THz_V4(F_Analog,H_err,PL_THz_GT,PL_tilde_THz_GT,B_T,K,M,General_Struct,Traj_data,0);
    

    Ave_Eff_H_mat_err = Eff_H_mat_err(:,:,:,1);
    Ave_Eff_H_mat_tilde_err = Eff_H_tilde_mat_err(:,:,:,1);


    H_mat = Ave_Eff_H_mat_err;
    H_tilde_mat = Ave_Eff_H_mat_tilde_err;

    THz_User_Dist_B1 = THz_User_Dist_GT(:,:,1) + THz_User_Dist_err;

    if n == 1


        idx = 1;

        alpha_n_1 = ones(B_T,K);

        alpha_n_mat(:,:,idx) = alpha_n_1;

    end

    idx = idx + 1;


    switch Func_Type

        case 'THzOn-WithCost'
            [SumRate_HO_THzOn(n),Rate_HO_THzOn(:,n),SumRate_THzOn(n),Rate_THzOn(:,n),Infeas_flg,Tmax(n),THz_SINR,alpha_n,W_mat_opt] = ...
                THzOn_Opt_Mob_Func(General_Struct,Traj_data,eps_conv,Pmax_THz,H_mat,H_tilde_mat,B_T,K,M,Rth_vec,THz_BW,...
                Clus_Siz_MBN,UA_flg,Max_Min_Init,BL_mat,eta_THz,alpha_n_1);
        case 'THzOn-Pen-1'
            HO_Pen = HO_Pen_1;
            [SumRate_HO_THzOn(n),Rate_HO_THzOn(:,n),SumRate_THzOn(n),Rate_THzOn(:,n),Infeas_flg,Tmax(n),THz_SINR,alpha_n,W_mat_opt] = ...
                THzOn_Opt_Mob_HO_P_Func(General_Struct,Traj_data,eps_conv,Pmax_THz,H_mat,H_tilde_mat,B_T,K,M,Rth_vec,THz_BW,...
                Clus_Siz_MBN,HO_Pen,UA_flg,Max_Min_Init,BL_mat,eta_THz,alpha_n_1);
       case 'THzOn-Pen-2'
            HO_Pen = HO_Pen_2;
            [SumRate_HO_THzOn(n),Rate_HO_THzOn(:,n),SumRate_THzOn(n),Rate_THzOn(:,n),Infeas_flg,Tmax(n),THz_SINR,alpha_n,W_mat_opt] = ...
                THzOn_Opt_Mob_HO_P_Func(General_Struct,Traj_data,eps_conv,Pmax_THz,H_mat,H_tilde_mat,B_T,K,M,Rth_vec,THz_BW,...
                Clus_Siz_MBN,HO_Pen,UA_flg,Max_Min_Init,BL_mat,eta_THz,alpha_n_1);
        case 'THzOn-NoCost'
            [SumRate_HO_THzOn(n),Rate_HO_THzOn(:,n),SumRate_THzOn(n),Rate_THzOn(:,n),Infeas_flg,Tmax(n),THz_SINR,alpha_n,W_mat_opt] = ...
            THzOn_Opt_Stat_Func(General_Struct,Traj_data,eps_conv,Pmax_THz,H_mat,H_tilde_mat,B_T,K,M,Rth_vec,THz_BW,...
            Clus_Siz_MBN,UA_flg,Max_Min_Init,BL_mat,eta_THz,alpha_n_1);
        case 'THzOn-B1'
            [SumRate_HO_THzOn(n),Rate_HO_THzOn(:,n),SumRate_THzOn(n),Rate_THzOn(:,n),Infeas_flg,Tmax(n),THz_SINR,alpha_n,W_mat_opt] = ...
            THzOn_B1_Func(General_Struct,Traj_data,eps_conv,Pmax_THz,H_mat,H_tilde_mat,B_T,K,M,Rth_vec,THz_BW,Clus_Siz_MBN,THz_User_Dist_B1,BL_mat,eta_THz,alpha_n_1);
    end

    [AveSumRate_HO_THz_err(n),AveSumRate_THz_err(n),AveRate_HO_THz_err(:,n),AveRate_THz_err(:,n),AveTHz_SINR_err(:,:,n)] = THz_AveRate_func(General_Struct,Traj_data,Eff_H_mat_err,Eff_H_tilde_mat_err,W_mat_opt,eta_THz,BL_mat,Ther_Noise_THz,THz_BW,THz_Hyb,alpha_n_1,alpha_n);

    AveRate_HO_THzOn_err_err(:,n) = AveRate_HO_THz_err(:,n);

    Rate_AvRate_HO_THzOn_err = (AveRate_HO_THzOn_err_err(:,n)).';

    AveSumRate_HO_THzOn_err_err(n) = sum(Rate_AvRate_HO_THzOn_err);

    if any(Rate_AvRate_HO_THzOn_err < Rth_vec)
        AveSumRate_HO_THzOn_err_err(n) = 0;
    end



    AveSumRate_THzOn_err_err(n) = AveSumRate_THz_err(n);
    AveRate_THzOn_err_err(:,n) = AveRate_THz_err(:,n);

    
    [SumRate_HO_THz_AveCh_err_GT(n),SumRate_THz_AveCh_err_GT(n),Rate_HO_THz_AveCh_err_GT(:,n),Rate_THz_AveCh_err_GT(:,n),THz_SINR_AveCh_err_GT(:,n)]  = THz_Rate_func(General_Struct,Traj_data,Eff_H_mat_err_GT(:,:,:,1),Eff_H_tilde_mat_err_GT(:,:,:,1),W_mat_opt,eta_THz,BL_mat,Ther_Noise_THz,THz_BW,THz_Hyb,alpha_n_1,alpha_n);

    SumRate_HO_THzOn_AveCh_err_GT(n) = SumRate_HO_THz_AveCh_err_GT(n);
    Rate_HO_THzOn_AveCh_err_GT(:,n) = Rate_HO_THz_AveCh_err_GT(:,n);

    Rate_HO_THzOn_AveCh_GT = (Rate_HO_THzOn_AveCh_err_GT(:,n)).';

    SumRate_HO_THzOn_AveCh_GT(n) = sum(Rate_HO_THzOn_AveCh_GT);

    if any(Rate_HO_THzOn_AveCh_GT < Rth_vec)
        SumRate_HO_THzOn_AveCh_GT(n) = 0;
    end



    SumRate_THzOn_AveCh_err_GT(n) = SumRate_THz_AveCh_err_GT(n);
    Rate_THzOn_AveCh_err_GT(:,n) = Rate_THz_AveCh_err_GT(:,n);


    [AveSumRate_HO_THz_err_GT(n),AveSumRate_THz_err_GT(n),AveRate_HO_THz_err_GT(:,n),AveRate_THz_err_GT(:,n),AveTHz_SINR_err_GT(:,:,n)] = THz_AveRate_func(General_Struct,Traj_data,Eff_H_mat_err_GT,Eff_H_tilde_mat_err_GT,W_mat_opt,eta_THz,BL_mat,Ther_Noise_THz,THz_BW,THz_Hyb,alpha_n_1,alpha_n);

    AveSumRate_HO_THzOn_err_GT(n) = AveSumRate_HO_THz_err_GT(n);
    AveRate_HO_THzOn_err_GT(:,n) = AveRate_HO_THz_err_GT(:,n);

    AveSumRate_THzOn_err_GT(n) = AveSumRate_THz_err_GT(n);
    AveRate_THzOn_err_GT(:,n) = AveRate_THz_err_GT(:,n);



    HO_Count_per_User(:,n) = (sum((1-alpha_n_1).*alpha_n));


    W_MAT(:,:,:,n) = W_mat_opt;

    alpha_n_1 = alpha_n > 0.5;
    alpha_n_mat(:,:,idx) = alpha_n > 0.5;

end





alpha_n_mat = alpha_n_mat > 0.5;

Assos_n_mat = [alpha_n_mat];

j = 0;
for i = 2:N+1
    j = j + 1;
    NumHO_THzOn(j) = sum(sum((1-Assos_n_mat(:,:,i-1)).*Assos_n_mat(:,:,i)));
end
NumHO_THzOn = sum(NumHO_THzOn);


F_HO_AveCh_err = mean(SumRate_HO_THzOn);


F_HO_AveRate_err = mean(AveSumRate_HO_THzOn_err_err);


F_HO_AveCh_GT = mean(SumRate_HO_THzOn_AveCh_GT);


F_HO_AveRate_GT = mean(AveSumRate_HO_THzOn_err_GT);


end





