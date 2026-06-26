function [A_Out_LP] = Assoc_Single_Conn(Assoc_params,B,S,N,BL_mat,SE_mat)

Idle_BS_flg = Assoc_params.Idle_BS_flg;
Gam_L = Assoc_params.Gam_L;
Gam_U = Assoc_params.Gam_U;


Block_mat_1 = ones(B,N);
Block_mat_1(BL_mat==0)= -100;
Block_mat = reshape(Block_mat_1,[B,1,N]);
Block_mat = repmat(Block_mat,[1,S,1]);
Block_mat = ones(B,S,N);
% Constraints matrices:

Q = repmat(eye(N,N),[S,1]);
Z = kron(eye(S),Q);
W = repmat(eye(N*S,N*S),[S,1]);
Z1 = kron(eye(B),Z);
W_mask = ones(B,B) - eye(B,B);
C = kron(W_mask,W) + Z1;

L1 = kron(eye(S,S),ones(N,1)');
D = kron(eye(B), L1);
E = repmat(eye(N,N),[1,B*S]);
F = -E;
k1 = ones(N*S*S*B,1);
k2 = ones(B*S,1);
k3 = Gam_U*ones(N,1);
k4 = -Gam_L*ones(N,1);


delta = 10;
T_Assoc = 20;

A_t_LP = 0.55*ones(B,S,N);

for t_LP = 1 : T_Assoc
    t_LP;
    Weight = Block_mat.*SE_mat - delta * (1-2*A_t_LP); 
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
        A_LP = NaN;
    else
        A_LP = permute(reshape(x,[N,S,B]), [3,2,1]);
    end

    Penalty_func_LP(t_LP) = sum(sum(sum(A_LP.*(1-2*A_t_LP) + A_t_LP.^2)));
    A_t_LP = A_LP;

    if Penalty_func_LP(t_LP) < 1e-3 || t_LP == 5
        A_Out_LP = A_LP > 0.5;
        break
    else
        A_t_LP = A_LP;
    end

end


end