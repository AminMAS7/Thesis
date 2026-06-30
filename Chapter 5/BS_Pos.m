function [THz_BS_Loc,RF_BS_Loc,BS_Loc] = BS_Pos(B_T,B_R,Channel_info)

L_x_max = Channel_info.L_x_max;
W_y_max = Channel_info.W_y_max;

% Positions loading:

%%%%% THz BS: 
if mod(B_T,2) == 0
    B_T_1 = B_T/2;
    B_T_2 = B_T/2;
else
    B_T_1 = floor(B_T/2);
    B_T_2 = B_T - B_T_1;
end
    

X_vec_1 = 5:L_x_max+5;
X_vec_2 = 20:L_x_max-20;

if B_T > 1
    THz_BS_Loc_1 = [quantile(X_vec_1,B_T_1)',0*ones(B_T_1,1)];
    THz_BS_Loc_2 = [quantile(X_vec_2,B_T_2)',W_y_max*ones(B_T_2,1)];
    
    if B_T_1 == 1
        THz_BS_Loc_1 = [median(X_vec_1),0*ones(B_T_1,1)];
    end
    
    if B_T_2 == 1
        THz_BS_Loc_2 = [median(X_vec_2),W_y_max*ones(B_T_2,1)];
    end
    
    THz_BS_Loc_Orig = [THz_BS_Loc_1;THz_BS_Loc_2];
elseif B_T == 1
    THz_BS_Loc_Orig = [median(X_vec_2),W_y_max*ones(B_T_2,1)];
end




%%%% RF BS:

if mod(B_R,2) == 0
    B_R_1 = B_R/2;
    B_R_2 = B_R/2;
else
    B_R_1 = floor(B_R/2);
    B_R_2 = B_R - B_R_1;
end


if B_R > 1
    RF_BS_Loc_1 = [quantile(X_vec_1,B_R_1)',(0)*ones(B_R_1,1)];
    RF_BS_Loc_2 = [quantile(X_vec_2,B_R_2)',(W_y_max)*ones(B_R_2,1)];
    
    if B_R_1 == 1
        RF_BS_Loc_1 = [median(X_vec_1),(0)*ones(B_R_1,1)];
    end
    
    if B_R_2 == 1
        RF_BS_Loc_2 = [median(X_vec_2),(W_y_max)*ones(B_R_2,1)];
    end
    
    RF_BS_Loc_Orig = [RF_BS_Loc_1;RF_BS_Loc_2];
elseif B_R == 1
    RF_BS_Loc_Orig = [median(X_vec_2),(0)*ones(B_R_2,1)];
end


THz_BS_Loc(:,1) = THz_BS_Loc_Orig(:,1) - min(THz_BS_Loc_Orig(:,1)) + 5;
THz_BS_Loc(:,2) = THz_BS_Loc_Orig(:,2);
RF_BS_Loc(:,1) = RF_BS_Loc_Orig(:,1) - min(THz_BS_Loc_Orig(:,1)) + 5;
RF_BS_Loc(:,2) = RF_BS_Loc_Orig(:,2);


BS_Loc = [THz_BS_Loc;RF_BS_Loc];

end