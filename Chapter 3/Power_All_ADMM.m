
function [Output_ADMM] = Power_All_ADMM(A,B,S,N,Ther_Noise_mat,h,h_bar,h_tilde,SystemParams,OptParams)



w = SystemParams.w;
Pmax = SystemParams.Pmax;
q = SystemParams.InvAntennaDirec;
k_r = SystemParams.k_r;
k_t = SystemParams.k_t;
% Channels:

Init_scale = OptParams.Init_scale;
T_Power_ADMM = OptParams.Outer_ADMM_Iterations;
eps_Power_ADMM = OptParams.Outer_ADMM_Accuracy;
T_ADMM = OptParams.Inner_ADMM_Iterations;
epsilon_AD = OptParams.Inner_ADMM_Accuracy;
zeta = OptParams.ADMM_Penalty_Fac;

% Initializing A,P,Gamma
Power_Init = Init_scale*(sqrt(Pmax/(S).*ones(B,S))); % Power initialization
P_bar = Power_Init;
for t = 1 : T_Power_ADMM

    Nu_ADMM = zeros(B,S);
    Z_ADMM = zeros(B,S);

    [Cum_Interf_mat] = q*Interference_calc(P_bar,h_bar,B,N,S);
    [Cum_Interf_HWI] = q*Interference_calc(P_bar,h,B,N,S);

    [Sum_SE_i_1(t),~,~] = WSR_Calc(A,h,h_tilde,P_bar,S,B,Ther_Noise_mat,Cum_Interf_mat,Cum_Interf_HWI,k_r,k_t);
    Power_squared = P_bar.^2;
    Num_Des_Power = reshape(Power_squared,[B,S,1]).*h; % Desired power for all users-BS-frequencies
    Absorb_Noise = reshape(Power_squared,[B,S,1]).*h_tilde; % Absorption noise for all users-BS-frequencies
    
    % Updating Gamma:
    Gamma_vec = (Num_Des_Power)./(Cum_Interf_mat + Absorb_Noise + (k_t^2 + k_r^2)*Num_Des_Power + k_r^2 *Cum_Interf_HWI + Ther_Noise_mat);
    
    % Updating Y for the fixed gamma obtained in the previous step:
    y_Num = sqrt(A.*(1+Gamma_vec).*reshape(Power_squared,[B,S,1]).*h);
    y_Dnom = Cum_Interf_mat + reshape(Power_squared,[B,S,1]).*(h_bar) + (k_t^2 + k_r^2)*Num_Des_Power + k_r^2 *Cum_Interf_HWI + Ther_Noise_mat;
    y_opt_mat =  y_Num./y_Dnom;
    [y_inter_cal_mat] = q * y_square(y_opt_mat,h_bar,B,N,S);
    [y_inter_cal_HWI] = q * y_square(y_opt_mat,h,B,N,S);

    for i = 1 : T_ADMM
        i;
        % Obtaining power
        Num_n_power = y_opt_mat.*sqrt(A.*(1+Gamma_vec).*h);
        Num_power = sum(Num_n_power,3) + 0.5*zeta*(Nu_ADMM + Z_ADMM);
        Denom_power = y_inter_cal_mat + sum((y_opt_mat.^2).*h_bar,3) + (k_t^2 + k_r^2)*sum((y_opt_mat.^2).*h,3) + k_r^2 * y_inter_cal_HWI + + 0.5*zeta;
        P_bar_ADMM = Num_power./Denom_power;
        sum(P_bar_ADMM.^2,2);
        % Obtaining the optimal dual variable:
        alpha_mat = P_bar_ADMM - Z_ADMM;
        b = sqrt(Pmax.*ones(B,1)./sum(alpha_mat.^2,2));
        Nu_1 = min(ones(B,1),b);
        Nu_ADMM = Nu_1 .* alpha_mat;
        
        % Updating Lagrange multiplier z:
        Z_ADMM = 1*(Nu_ADMM - P_bar_ADMM) + Z_ADMM;

        T_inner_Cove_ADMM(t) = i;

        
        sqrt(sum((P_bar_ADMM - Nu_ADMM).^2,"all")) / sqrt(sum((P_bar_ADMM).^2,"all"));
        if sqrt(sum((P_bar_ADMM - Nu_ADMM).^2,"all")) / sqrt(sum((P_bar_ADMM).^2,"all")) <= epsilon_AD
            P_bar = Nu_ADMM;
            break
        end

    end

    P_bar_opt = P_bar;
    P_bar_opt(isnan(P_bar_opt)) = 0;
    P_orig = P_bar.^2;
    POW_Total_Vec(:,t) = sum(P_orig,2);
    
    [Cum_Interf_mat] = q*Interference_calc(P_bar,h_bar,B,N,S);
    [Cum_Interf_HWI] = q*Interference_calc(P_bar,h,B,N,S);

    [Sum_SE_i(t),SE_mat,Sum_SE_user] = WSR_Calc(A,h,h_tilde,P_bar,S,B,Ther_Noise_mat,Cum_Interf_mat,Cum_Interf_HWI,k_r,k_t);

    Power_opt_converge = P_orig;
    WSR_f = w * Sum_SE_i(t);
    WSR_Vec = w * Sum_SE_i_1;
    WSR_mat_f = w * SE_mat;
    User_rate_f = w * Sum_SE_user;
    Cum_Interf_mat_f = Cum_Interf_mat;
    Total_Power = sum(P_orig,2);
    Converg_pow = t;

    Output_ADMM.Power_Opt_ADMM = Power_opt_converge;
    Output_ADMM.WSR_ADMM = WSR_f;
    Output_ADMM.Sum_SE_ADMM = Sum_SE_i(t);
    Output_ADMM.WSR_mat_ADMM = WSR_mat_f;
    Output_ADMM.SE_mat_ADMM = SE_mat;
    Output_ADMM.User_rate_ADMM = User_rate_f;
    Output_ADMM.Cum_Interf_mat_ADMM = Cum_Interf_mat_f;
    Output_ADMM.Total_Power_ADMM = Total_Power;
    Output_ADMM.Converg_FP_ADMM = Converg_pow;
    Output_ADMM.WSR_Vec_ADMM = WSR_Vec;
    Output_ADMM.Sum_SE_Vec_ADMM = Sum_SE_i_1;
    Output_ADMM.Converge_Inner_ADMM = T_inner_Cove_ADMM;
    Output_ADMM.Total_Power_Vec_ADMM = POW_Total_Vec;

    if abs(Sum_SE_i(t) - Sum_SE_i_1(t)) < eps_Power_ADMM || t == T_Power_ADMM
        break
    end
end
end


function [Sum_SE,SE_mat,Sum_SE_user] = WSR_Calc(A,h,h_tilde,Power,S,B,Ther_Noise_mat,Cum_Interf_mat,Cum_Interf_HWI,k_r,k_t)
Power_squared = Power.^2;
% The total sum rate function:
Num_Des_Power = reshape(Power_squared,[B,S,1]).*h; % Desired power for all users-BS-frequencies
Absorb_Noise = reshape(Power_squared,[B,S,1]).*h_tilde; % Absorption noise for all users-BS-frequencies
Gamma_vec = (Num_Des_Power)./(Cum_Interf_mat + Absorb_Noise + (k_t^2 + k_r^2)*Num_Des_Power + k_r^2 *Cum_Interf_HWI  + Ther_Noise_mat);
SE_mat = log2(exp(1)) * A.*log(1+Gamma_vec); 
Sum_SE_user = sum(SE_mat,[1,2]);
Sum_SE = sum(SE_mat,'all');
end

function [I] = Interference_calc(P_bar,h_bar,B,N,S)
P = (P_bar.^2);
H_bar = reshape(h_bar,[B,1,S,N]);
G_bar = repmat(H_bar,[1,B,1,1]);
E = G_bar .* reshape(P,[B,1,S,1]);
r_mask = reshape(eye(B,B),[B,B,1,1]);
i_mask = ones(B,B,S,N) - r_mask;
I = squeeze(sum(E.*i_mask,1));
end

function [I_Y_out] = y_square(y_opt_mat,h_bar,B,N,S)
y_opt_mat_pow = y_opt_mat.^2;
Y = reshape(y_opt_mat_pow,[B,1,S,N]);
Y_bar = repmat(Y,[1,B,1,1]);
H = reshape(h_bar,[1,B,S,N]);
H_bar = repmat(H,[B,1,1,1]);
E = Y_bar.*H_bar;
r_mask = reshape(eye(B,B),[B,B,1,1]);
i_mask = ones(B,B,S,N) - r_mask;
I_Y_bar = E .* i_mask;
I_Y = squeeze(sum(I_Y_bar,1));
I_Y_out = sum(I_Y,3);
end
