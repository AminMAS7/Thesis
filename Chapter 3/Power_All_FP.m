
function [Output_FP] = Power_All_FP(A,B,S,N,Ther_Noise_mat,SystemParams,h,h_bar,h_tilde,OptParams,BiSecParam)


k_r = SystemParams.k_r;
k_t = SystemParams.k_t;

w = SystemParams.w;
Pmax = SystemParams.Pmax;
Pmax_vec = Pmax*ones(1,B);
% Ther_Noise_mat = SystemParams.Ther_Noise_mat;
q = SystemParams.InvAntennaDirec;

% Channels:

Init_scale = OptParams.Init_scale;
T_Power = OptParams.FP_Iterations;
eps_Power = OptParams.eps_Power;


% Initializing A,P,Gamma
Power_Init = Init_scale*(sqrt(Pmax/(S).*ones(B,S))); % Power initialization
P_bar = Power_Init;
for i = 1 : T_Power
    i;
    [Cum_Interf_mat] = q*Interference_calc(P_bar,h_bar,B,N,S);
    [Cum_Interf_HWI] = q*Interference_calc(P_bar,h,B,N,S);

    [Sum_SE_i_1(i),~,~] = WSR_Calc(A,h,h_tilde,P_bar,S,B,Ther_Noise_mat,Cum_Interf_mat,Cum_Interf_HWI,k_r,k_t);
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
    % Updating Power for fixed Gamma and Y:

    % Bi-section search for Lagrange multiplier (mu_b) for each BS
    mu_b_Min = BiSecParam.mu_b_Min;
    mu_b_Max = BiSecParam.mu_b_Max;
    errTol = BiSecParam.errTol_Bisec;
    mu_b_vec = zeros(B,1);
    for b = 1 : B
        if Power_Calc_Bisec(A,y_opt_mat,S,Gamma_vec,h,h_bar,y_inter_cal_mat,y_inter_cal_HWI,k_t,k_r,0,Pmax_vec,b) <= 0
            mu_b_vec(b,1) = 0;
            continue
        end
        [mu_b_out] = BiSection_FP_Func(A,h,h_bar,S,Gamma_vec,y_opt_mat,b,Pmax_vec,mu_b_Min,mu_b_Max,errTol,y_inter_cal_mat,y_inter_cal_HWI,k_t,k_r);
        mu_b_vec(b,1) = mu_b_out;
    end
    mu_b_mat = repmat(mu_b_vec,[1,S]);

    Num_n_power = y_opt_mat.*sqrt(A.*(1+Gamma_vec).*h);
    Num_power = sum(Num_n_power,3);
    Denom_power = y_inter_cal_mat + sum((y_opt_mat.^2).*h_bar,3) + (k_t^2 + k_r^2)*sum((y_opt_mat.^2).*h,3) + k_r^2 * y_inter_cal_HWI + mu_b_mat;
    P_bar = Num_power./Denom_power;
    P_bar_opt = P_bar;
    P_bar_opt(isnan(P_bar_opt)) = 0;
    P_bar = P_bar_opt;
    P_orig = P_bar.^2;
    POW_Total(:,i) = sum(P_orig,2);
   [Cum_Interf_mat] = q*Interference_calc(P_bar,h_bar,B,N,S);
   [Cum_Interf_HWI] = q*Interference_calc(P_bar,h,B,N,S);

   [Sum_SE_i(i),SE_mat,Sum_SE_user] = WSR_Calc(A,h,h_tilde,P_bar,S,B,Ther_Noise_mat,Cum_Interf_mat,Cum_Interf_HWI,k_r,k_t);

    Power_opt_converge = P_orig;
    WSR_f = w * Sum_SE_i(i);
    WSR_Vec = w * Sum_SE_i_1;
    WSR_mat_f = w * SE_mat;
    User_rate_f = w * Sum_SE_user;
    Cum_Interf_mat_f = Cum_Interf_mat;
    Total_Power = sum(P_orig,2);
    Converg_pow = i;
    Output_FP = {};
    Output_FP.Power_FP = Power_opt_converge;
    Output_FP.Sum_SE_vec_FP = Sum_SE_i_1;
    Output_FP.WSR_Vec_FP = WSR_Vec;
    Output_FP.Sum_SE_FP = Sum_SE_i(i);
    Output_FP.Cum_Interf_mat_FP = Cum_Interf_mat_f;
    Output_FP.WSR_FP = WSR_f;
    Output_FP.SE_mat_FP = SE_mat;
    Output_FP.WSR_mat_FP = WSR_mat_f;
    Output_FP.User_rate_FP = User_rate_f;
    Output_FP.Total_Power_FP = Total_Power;
    Output_FP.Converg_FP = Converg_pow;

    if abs(Sum_SE_i(i) - Sum_SE_i_1(i)) < eps_Power || i == T_Power
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

function [mu_b_out] = BiSection_FP_Func(A,h,h_bar,S,Gamma_vec,y_opt_mat,b,Pmax_vec,mu_b_Min,mu_b_Max,errTol,y_inter_cal_mat,y_inter_cal_HWI,k_t,k_r)
    a1 = Power_Calc_Bisec(A,y_opt_mat,S,Gamma_vec,h,h_bar,y_inter_cal_mat,y_inter_cal_HWI,k_t,k_r,mu_b_Min,Pmax_vec,b);
    a2 = Power_Calc_Bisec(A,y_opt_mat,S,Gamma_vec,h,h_bar,y_inter_cal_mat,y_inter_cal_HWI,k_t,k_r,mu_b_Max,Pmax_vec,b);
    if a1 * a2 > 0
        disp('Wrong choice')
    else
        mu_b_opt = (mu_b_Min + mu_b_Max)/2;
        err = abs(Power_Calc_Bisec(A,y_opt_mat,S,Gamma_vec,h,h_bar,y_inter_cal_mat,y_inter_cal_HWI,k_t,k_r,mu_b_opt,Pmax_vec,b));
        count = 0;
        while err > errTol
            count = count + 1;
            if a1 * Power_Calc_Bisec(A,y_opt_mat,S,Gamma_vec,h,h_bar,y_inter_cal_mat,y_inter_cal_HWI,k_t,k_r,mu_b_opt,Pmax_vec,b) < 0
                mu_b_Max = mu_b_opt;
            else
                mu_b_Min = mu_b_opt;
            end
            mu_b_opt = (mu_b_Max + mu_b_Min)/2;
            err = abs(Power_Calc_Bisec(A,y_opt_mat,S,Gamma_vec,h,h_bar,y_inter_cal_mat,y_inter_cal_HWI,k_t,k_r,mu_b_opt,Pmax_vec,b));
            a1 = Power_Calc_Bisec(A,y_opt_mat,S,Gamma_vec,h,h_bar,y_inter_cal_mat,y_inter_cal_HWI,k_t,k_r,mu_b_Min,Pmax_vec,b);
            if count > 1e6
                disp('Error in bi-section search')
            break
            end
        end
    end
    mu_b_out = mu_b_opt;
end
    

function [Out] = Power_Calc_Bisec(A,y_opt_mat,S,Gamma_vec,h,h_bar,y_inter_cal_mat,y_inter_cal_HWI,k_t,k_r,mu_b,Pmax_vec,b)
Num_n = y_opt_mat.*sqrt(A.*(1+Gamma_vec).*h);
Num = sum(Num_n,3);
pow_Num_b = Num(b,:);
Denom_1 = y_inter_cal_mat + sum((y_opt_mat.^2).*h_bar,3) + (k_t^2 + k_r^2)*sum((y_opt_mat.^2).*h,3) + k_r^2 * y_inter_cal_HWI;
Denom_1_1 = Denom_1(b,:);
pow_Denom_b = (Denom_1_1 + mu_b);
Pow_square = (pow_Num_b ./ pow_Denom_b).^2;
Pow_square(isnan(Pow_square)) = 0;
Pow_out = sum(Pow_square,2);
Out = Pow_out - Pmax_vec(b);
end

