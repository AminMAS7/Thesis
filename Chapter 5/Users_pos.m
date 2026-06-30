
function [Users_x_mat,Users_Init_y] = Users_pos(K,THz_BS_Loc,N,tau_eps,v,Users_in_X,Channel_info)

W_y_max = Channel_info.W_y_max;
W_min_THz = Channel_info.W_min_THz;



Users_x_Incre = repmat((0:N-1)*tau_eps*v,[K,1]);
Users_x_mat_Orig = repmat(Users_in_X,[1,N]) + Users_x_Incre;

% Users_x_mat = mod(Users_x_mat_Orig,max(THz_BS_Loc(:,1))+10);

Users_x_mat = Users_x_mat_Orig;

Users_Init_y = linspace(0+W_min_THz,W_y_max-W_min_THz,K);

end
