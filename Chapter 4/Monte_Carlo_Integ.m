% Author: Mohammad Amin Saeidi
% Date: August 2023
%
% Usage and Redistribution Notice:
% - This code is provided for [educational/research/personal] purposes only.
% - Any use, modification, redistribution, or other use of this code
%   requires explicit permission from the author.
% - To request permission, please send an email to [amin96a@yorku.ca or m.amin.saeidi96@gmail.com].

function [max_rate_Intg_mean,RFassoc_Intg_mean,THassoc_Intg_mean] = Monte_Carlo_Integ(NInteg,T,Params,Power_flg)
max_rate_Intg = zeros(T,1);
RFassoc_Intg = zeros(T,1);
THassoc_Intg = zeros(T,1);
Xu = 0;
Yu = 0;
R_max = Params.R_max;

parfor i = 1 : T
    % Loading BSs positions for Integrated deployment
    rRF_THzbs = R_max*sqrt(rand(NInteg,1));
    thetaRF_TH = 2*pi*rand(NInteg,1);
    XbHyb = rRF_THzbs.*cos(thetaRF_TH);
    YbHyb = rRF_THzbs.*sin(thetaRF_TH);

    % Distacnces from RF_THz BSs
    [Xmp_Tmat_Integ, Xp_Tmat_Integ] = meshgrid(Xu,XbHyb);
    [Ymp_Tmat_Integ, Yp_Tmat_Integ] = meshgrid(Yu,YbHyb);
    D_ue_RF_THzbs = sqrt((Xmp_Tmat_Integ-Xp_Tmat_Integ).^2 + (Ymp_Tmat_Integ-Yp_Tmat_Integ).^2);
    % Fading channels for RF mode of Hybrid BSs
    fadeRand_Integ = exprnd(1,NInteg,1);

    % Computing rates
    % Integrated:
    [max_rate_Intg(i),RFassoc_Intg(i),THassoc_Intg(i)] = Integ_Rate(D_ue_RF_THzbs,Params,fadeRand_Integ,Power_flg);

end

% Integrated Average rates:
max_rate_Intg_mean = mean(max_rate_Intg);
RFassoc_Intg_mean = mean(RFassoc_Intg);
THassoc_Intg_mean = mean(THassoc_Intg);


end