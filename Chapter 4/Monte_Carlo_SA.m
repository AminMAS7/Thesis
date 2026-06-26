% Author: Mohammad Amin Saeidi
% Date: August 2023
%
% Usage and Redistribution Notice:
% - This code is provided for [educational/research/personal] purposes only.
% - Any use, modification, redistribution, or other use of this code
%   requires explicit permission from the author.
% - To request permission, please send an email to [amin96a@yorku.ca or m.amin.saeidi96@gmail.com].

function [max_rate_SA_mean,RFassoc_SA_mean,THassoc_SA_mean] = Monte_Carlo_SA(NR,NT,T,Params,Power_flg)
max_rate_SA = zeros(T,1);
RFassoc_SA = zeros(T,1);
THassoc_SA = zeros(T,1);
Xu = 0;
Yu = 0;
R_max = Params.R_max;

for i = 1 : T
    % Loading BSs positions for Stand-alone deployment
    % THz:
    rTHzbs = R_max*sqrt(rand(NT,1));
    thetaTb = 2*pi*rand(NT,1);
    Xb = rTHzbs.*cos(thetaTb);
    Yb = rTHzbs.*sin(thetaTb);
    % RF:
    rRFzbs = R_max*sqrt(rand(NR,1));
    thetaRb = 2*pi*rand(NR,1);
    Xrb = rRFzbs.*cos(thetaRb);
    Yrb = rRFzbs.*sin(thetaRb);
    % Distacnces from THz BSs
    [Xmp_Tmat, Xp_Tmat] = meshgrid(Xu,Xb);
    [Ymp_Tmat, Yp_Tmat] = meshgrid(Yu,Yb);
    D_ue_Tbs = sqrt((Xmp_Tmat-Xp_Tmat).^2 + (Ymp_Tmat-Yp_Tmat).^2);
    % Distacnces from RF BSs
    [Xmp_Rmat, Xp_Rmat] = meshgrid(Xu,Xrb);
    [Ymp_Rmat, Yp_Rmat] = meshgrid(Yu,Yrb);
    D_ue_Rbs = sqrt((Xmp_Rmat-Xp_Rmat).^2 + (Ymp_Rmat-Yp_Rmat).^2);
    % Fading channel for RF BSs
    fadeRand = exprnd(1,NR,1); % Channel fade RV based on number of users of interest (1) and number of RF BSs

    % Computing rates
    % Stand-alone:
    if NR == 0
        D_ue_Rbs = 1e100;
        fadeRand = 0;
    end
    if NT == 0
        D_ue_Tbs = 1e100;
    end
    
    [max_rate_SA(i),RFassoc_SA(i),THassoc_SA(i)] = Stand_alone_Rate(D_ue_Rbs,D_ue_Tbs,Params,fadeRand,Power_flg);
    
end

% Stand-alone Average rates:
max_rate_SA_mean = mean(max_rate_SA);
RFassoc_SA_mean = mean(RFassoc_SA);
THassoc_SA_mean = mean(THassoc_SA);

end