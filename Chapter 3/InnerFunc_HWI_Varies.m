function [WSR_main_FP_Rnd_Assoc,Mean_Conn_Order_Rnd_Assoc,WSR_EqPow,Mean_Conn_Order_Eq_Pow,...
    WSR_main_FP,Mean_Conn_Order_FP,WSR_Main_ADMM,Mean_Conn_Order_ADMM,WSR_main_FP_Sing_Conn,Mean_Conn_Order_FP_Sing_Conn] = ...
    InnerFunc_HWI_Varies(HI_Vec,Assoc_params,S,N_Orig,Ther_Noise,OptParams,BiSecParam,...
        SystemParams,B,seed)

alpha = SystemParams.alpha;
k_f = SystemParams.KF;

eta = 0.003;


for i = 1 : length(HI_Vec)
    i;
    k_HWI = HI_Vec(i);
    SystemParams.k_r = k_HWI;
    SystemParams.k_t = k_HWI;

    for j = 1 
        j;
        rng(seed)
        Fading_mat = exprnd(1,B,S,N_Orig);
        Ther_Noise_mat = repmat(Ther_Noise,[B,1,N_Orig]); % Thermal noise for all BS, users and frequencies

        R_max = 30; % The maximum radius of the users
        [D_user_BS_2d,~] = BS_User_Distances_mean(B,N_Orig,R_max,0,0,seed);
        
        % Channles: 
        %%%%% Blockage:

        Blockage_Prob = exp(-D_user_BS_2d*eta);
        BL_mat = rand(size(D_user_BS_2d)) < Blockage_Prob;
        Outage_Users_idx = find(sum(BL_mat,1)*S < Assoc_params.Gam_L);
        D_user_BS_2d(:,Outage_Users_idx)= [];
        BL_mat(:,Outage_Users_idx) = [];
        [~, N] = size(D_user_BS_2d);
        N;
        Block_mat_1 = reshape(BL_mat,[B,1,N]);
        Block_mat_1 = repmat(Block_mat_1,[1,S,1]);
        Num_Infeasible_Users = length(Outage_Users_idx);
        
        D_user_BS = reshape(D_user_BS_2d,[B,1,N]);
        % Channles:
        Molec_Loss = exp(-k_f.*repmat(D_user_BS,[1,S,1]));
        h = ((repmat(D_user_BS,[1,S,1])).^(alpha)).*Molec_Loss.*Fading_mat;
        h_bar = ((repmat(D_user_BS,[1,S,1])).^(alpha)).*Fading_mat;
        h_tilde = ((repmat(D_user_BS,[1,S,1])).^(alpha)).*(1-Molec_Loss).*Fading_mat;
        h = h .*Block_mat_1;
        h_bar = h_bar.*Block_mat_1;
        h_tilde = h_tilde.*Block_mat_1;

        % Initial association with equal power allocation:
        [SE_mat_Eq_Power] = SE_mat_Func(SystemParams.Pmax/S*ones(B,S),B,S,N,Ther_Noise_mat,SystemParams,h,h_bar,h_tilde);
        [A_init_Eq_Pow] = Association(Assoc_params,B,S,N,BL_mat,SE_mat_Eq_Power);
        Mean_Conn_Order_Eq_Pow(i,j) = mean(sum((squeeze(sum(A_init_Eq_Pow,2)) > 0),1));


        %%%% Equal Power Allocation:
        [SE_mat_EqPow_MRT] = SE_mat_Func(SystemParams.Pmax/S*ones(B,S),B,S,N,Ther_Noise_mat,SystemParams,h,h_bar,h_tilde);
        WSR_EqPow(i,j) = sum(A_init_Eq_Pow.*SE_mat_EqPow_MRT*SystemParams.w,'all');

        
        %%%%%%%% Random Association:
        [A_init_Rnd] = Association(Assoc_params,B,S,N,BL_mat,rand(B,S,N));
        Mean_Conn_Order_Rnd_Assoc(i,j) = mean(sum((squeeze(sum(A_init_Rnd,2)) > 0),1));


        [Output_RndAss] = Main_Direct_FP_Rnd_Assoc(A_init_Rnd,B,S,N,Ther_Noise_mat,SystemParams,h,h_bar,h_tilde,BL_mat,OptParams,BiSecParam);        
        WSR_main_FP_Rnd_Assoc(i,j) = Output_RndAss.WSR_FP_main;
        


        %%%%%%%%%%%%%%%%%%%%  Main FP
        [Output_Main_FP] = Main_Direct_FP(A_init_Eq_Pow,B,S,N,Ther_Noise_mat,Assoc_params,SystemParams,h,h_bar,h_tilde,BL_mat,OptParams,BiSecParam);

        Mean_Conn_Order_FP(i,j) = mean(sum((squeeze(sum(Output_Main_FP.Assoc_mat_FP_main,2)) > 0),1));   
        WSR_main_FP(i,j) = Output_Main_FP.WSR_FP_main;

        


        %%%%%%%%%%%%%%% Single Connectivity
        A_init_Sing_Conn = Assoc_Single_Conn(Assoc_params,B,S,N,BL_mat,SE_mat_Eq_Power);
        [Output_Main_FP_Sing_Conn] = Main_Direct_FP_Sing_Conn(A_init_Sing_Conn,B,S,N,Ther_Noise_mat,Assoc_params,SystemParams,h,h_bar,h_tilde,BL_mat,OptParams,BiSecParam);
        Mean_Conn_Order_FP_Sing_Conn(i,j) = mean(sum((squeeze(sum(Output_Main_FP_Sing_Conn.Assoc_mat_FP_main,2)) > 0),1));  

        WSR_main_FP_Sing_Conn(i,j) = Output_Main_FP_Sing_Conn.WSR_FP_main;


        %%%%%%%%%%%%% ADMM:


        [Output_Main_ADMM] = Main_ADMM_based_FP(A_init_Eq_Pow,B,S,N,Ther_Noise_mat,Assoc_params,SystemParams,h,h_bar,h_tilde,BL_mat,OptParams);
        Mean_Conn_Order_ADMM(i,j) = mean(sum((squeeze(sum(Output_Main_ADMM.Assoc_mat_ADMM_main,2)) > 0),1));  

        WSR_Main_ADMM(i,j) = Output_Main_ADMM.WSR_ADMM_main;

    end

end
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
if isempty(x)
    Assoc_Out = NaN;
else
    Assoc_Out = permute(reshape(x,[N,S,B]), [3,2,1]);
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


