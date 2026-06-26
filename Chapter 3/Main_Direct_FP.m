function [Output_Main_FP] = Main_Direct_FP(A_init_Decoupled,B,S,N,Ther_Noise_mat,Assoc_params,SystemParams,h,h_bar,h_tilde,BL_mat,OptParams,BiSecParam)



T_main_FP = OptParams.T_main_Direct_FP;
A = A_init_Decoupled;

Sum_SE_main_FP = zeros(1,T_main_FP);
for t = 2 : T_main_FP
    if t == 2
        A_hat = A;
    else
        A_hat = A_2D;
    end

    [Output_PowerAll_FP] = Power_All_FP(A_hat,B,S,N,Ther_Noise_mat,SystemParams,h,h_bar,h_tilde,OptParams,BiSecParam);
    Power_Opt_FP = Output_PowerAll_FP.Power_FP;
    WSR_main_FP(t) = Output_PowerAll_FP.WSR_FP;
    Sum_SE_main_FP(t) = Output_PowerAll_FP.Sum_SE_FP;

    [SE_mat_2D] = SE_mat_Func(Power_Opt_FP,B,S,N,Ther_Noise_mat,SystemParams,h,h_bar,h_tilde);
    [A_2D] = Association(Assoc_params,B,S,N,BL_mat,SE_mat_2D);
    %A_2D(Block_mat==0)=0;

    Output_Main_FP.Assoc_mat_FP_main = A_2D;
    Output_Main_FP.Power_Opt_FP_main = Power_Opt_FP;
    Output_Main_FP.WSR_Vec_FP_main = WSR_main_FP;
    Output_Main_FP.Sum_SE_Vec_FP_main = Sum_SE_main_FP;
    Output_Main_FP.Sum_SE_FP_main = Sum_SE_main_FP(t);
    Output_Main_FP.WSR_FP_main = WSR_main_FP(t);
    Output_Main_FP.OutPut_Pow_All_Direct_FP = Output_PowerAll_FP;

    
    if abs(Sum_SE_main_FP(t) - Sum_SE_main_FP(t-1)) < OptParams.eps_main_FP
        break
    end


end

end

function [SE_mat] = SE_mat_Func(Power,B,S,N,Ther_Noise_mat,SystemParams,h,h_bar,h_tilde)

k_r = SystemParams.k_r;
k_t = SystemParams.k_t;
q = SystemParams.InvAntennaDirec;
% Channels:

P = Power;
H_bar = reshape(h_bar,[B,1,S,N]);
G_bar = repmat(H_bar,[1,B,1,1]);
E = G_bar .* reshape(P,[B,1,S,1]);
r_mask = reshape(eye(B,B),[B,B,1,1]);
i_mask = ones(B,B,S,N) - r_mask;
I = squeeze(sum(E.*i_mask,1));
Cum_Interf_mat = q * I;


H = reshape(h,[B,1,S,N]);
G = repmat(H,[1,B,1,1]);
E = G .* reshape(P,[B,1,S,1]);
r_mask = reshape(eye(B,B),[B,B,1,1]);
i_mask = ones(B,B,S,N) - r_mask;
I_HWI = squeeze(sum(E.*i_mask,1));
Cum_Interf_HWI = q * I_HWI;


% The total sum rate function:
Num_Des_Power = reshape(Power,[B,S,1]).*h; % Desired power for all users-BS-frequencies
Absorb_Noise = reshape(Power,[B,S,1]).*h_tilde; % Absorption noise for all users-BS-frequencies
Gamma_vec = (Num_Des_Power)./(Cum_Interf_mat + Absorb_Noise + (k_t^2 + k_r^2)*Num_Des_Power + k_r^2 *Cum_Interf_HWI  + Ther_Noise_mat);
SE_mat = log2(exp(1))*log(1+Gamma_vec); 
end

function [Assoc_Out] = Association(Assoc_params,B,S,N,BL_mat,SE_mat)
Idle_BS_flg = Assoc_params.Idle_BS_flg;
Gam_L = Assoc_params.Gam_L;
Gam_U = Assoc_params.Gam_U;


Block_mat_1 = ones(B,N);
Block_mat_1(BL_mat==0)= -100;
Block_mat = reshape(Block_mat_1,[B,1,N]);
Block_mat = repmat(Block_mat,[1,S,1]);
Block_mat = ones(B,S,N);
% Constraints matrices:
C = kron(ones(B,1)',eye(N*S,N*S));
L1 = kron(eye(S,S),ones(N,1)');
D = kron(eye(B), L1);
E = repmat(eye(N,N),[1,B*S]);
F = -E;
k1 = ones(N*S,1);
k2 = ones(B*S,1);
k3 = Gam_U*ones(N,1);
k4 = -Gam_L*ones(N,1);

Weight = Block_mat.*SE_mat; 
%Weight = SE_mat; 
W_vec = reshape(permute(Weight,[3 2 1]), [B*S*N, 1]);
lb = zeros(size(W_vec)); 
ub = ones(size(W_vec));
if Idle_BS_flg == 0
    K_Cont_Ieq = [k1;k3;k4];
    T_Cont_Ieq = [C;E;F];
    options = optimoptions('linprog', 'Algorithm', 'dual-simplex',Display='none');
    [x, ~, ~, ~] = linprog(-W_vec, T_Cont_Ieq, K_Cont_Ieq, D, k2, lb, ub, options);
elseif Idle_BS_flg == 1
    K = [k1;k2;k3;k4];
    T = [C;D;E;F];
    options = optimoptions('linprog', 'Algorithm', 'dual-simplex',Display='none');
    [x, ~, ~, ~] = linprog(-W_vec, T, K, [], [], lb, ub, options);
end

Assoc_Out = permute(reshape(x,[N,S,B]), [3,2,1]);

end