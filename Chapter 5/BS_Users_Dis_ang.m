function [THz_User_Dist,AoD_THz_Users_Angle] = BS_Users_Dis_ang(Users_x_mat,Users_Init_y,THz_BS_Loc,H_min)


% pre‑sizes
[K,N]      = size(Users_x_mat);
B          = size(THz_BS_Loc,1);

% --- build 3D arrays ---

% 1×K×N array for user x‑coords
users_x3   = reshape(Users_x_mat, [1, K, N]);

% 1×K×N array for user y‑coords: replicate the 1×K vector along the 3rd dim
users_y3   = reshape(Users_Init_y, [1, K, 1]);
users_y3   = repmat(users_y3, [1, 1, N]);

% B×1×1 arrays for BS coords (will broadcast across K×N)
bs_x3      = reshape(THz_BS_Loc(:,1), [B, 1, 1]);
bs_y3      = reshape(THz_BS_Loc(:,2), [B, 1, 1]);

% --- deltas (B×K×N) ---
deltaX     = users_x3 - bs_x3;
deltaY     = users_y3 - bs_y3;

% --- distance (B×K×N) ---
THz_User_Dist = sqrt( H_min^2 + deltaX.^2 + deltaY.^2 );

% --- AoD angle in radians (B×K×N) ---

AoD_THz_Users_Angle = abs( atan2(deltaY, deltaX) );



end