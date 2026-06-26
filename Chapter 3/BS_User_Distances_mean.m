% The average of BS-User distances: 
function [D_user_BS,AoD_THz_Users_Angle] = BS_User_Distances_mean(B,N,R_max,xBS,yBS,seed)


% xu = linspace(-R_max,R_max,N);
% yu = linspace(-R_max,R_max,N);
% [XU_mat, XBS_mat] = meshgrid(xu,xBS);
% [YU_mat, YBS_mat] = meshgrid(yu,yBS);
% D_user_BS = sqrt((XU_mat-XBS_mat).^2 + (YU_mat-YBS_mat).^2); % D_user_BS(b,n) = d(b,n)

rng(seed)
Rmin = 4;
% Rmax = 75;
c= 1;
while true
    c = c + 1;
    rUsers = R_max*sqrt(rand(N,1));
    thetaUsers = 2*pi*rand(N,1);
    xu = rUsers.*cos(thetaUsers);
    yu = rUsers.*sin(thetaUsers);
%     [XU_mat, XBS_mat] = meshgrid(xu,xBS);
%     [YU_mat, YBS_mat] = meshgrid(yu,yBS);

    rBS = R_max*sqrt(rand(B,1));
    thetaBS = 2*pi*rand(B,1);
    xBS = rBS.*cos(thetaBS);
    yBS = rBS.*sin(thetaBS);
    [XU_mat, XBS_mat] = meshgrid(xu,xBS);
    [YU_mat, YBS_mat] = meshgrid(yu,yBS);
    D_user_BS = sqrt((XU_mat-XBS_mat).^2 + (YU_mat-YBS_mat).^2); % D_user_BS(b,n) = d(b,n)

    AoD_THz_Users_Angle = atan(abs(XU_mat-XBS_mat)./abs(YU_mat-YBS_mat));

    % scatter(xu,yu)
    % hold on
    % scatter(xBS,yBS)

    flg = min(D_user_BS,[],1);
%     flg_2 = max(D_user_BS,[],1);
    if all(flg > Rmin)
        break
%         if all(flg_2 < Rmax)
%             break
%         end
    end
%     Small_dis = find(D_user_BS < 1);
%     D_user_BS(Small_dis) = 1;
end

end

