% Author: Mohammad Amin Saeidi
% Date: August 2023
%
% Usage and Redistribution Notice:
% - This code is provided for [educational/research/personal] purposes only.
% - Any use, modification, redistribution, or other use of this code
%   requires explicit permission from the author.
% - To request permission, please send an email to [amin96a@yorku.ca or m.amin.saeidi96@gmail.com].

% MonteCarlo Funtion

function [Rate_SA,Sum_rate_SA,Sum_SE_SA,SE_SA,RF_ASSOC_SA,THz_ASSOC_SA,Rate_Intg,Sum_rate_Intg,Sum_SE_Intg,SE_Intg,RF_ASSOC_Intg,THz_ASSOC_Intg] = MonteCarlo_func(N,T,Params,R_max)
Xu = 0;
Yu = 0;

max_rate_SA = ones(T,1);
RFassoc_SA = ones(T,1);
THassoc_SA = ones(T,1);
max_rate_Int = ones(T,1);
RFassoc_Int = ones(T,1);
THassoc_Int = ones(T,1);

parfor i = 1 : T
    % Loading BSs positions for Stand-alone deployment
    % THz:
    rTHzbs = R_max*sqrt(rand(N,1));
    thetaTb = 2*pi*rand(N,1);
    Xb = rTHzbs.*cos(thetaTb);
    Yb = rTHzbs.*sin(thetaTb);
    % RF:
    rRFzbs = R_max*sqrt(rand(N,1));
    thetaRb = 2*pi*rand(N,1);
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
    fadeRand = exprnd(1,N,1); % Channel fade RV based on number of users of interest (1) and number of RF BSs

    % Loading BSs positions for Integrated deployment
    rRF_THzbs = R_max*sqrt(rand(N,1));
    thetaRF_TH = 2*pi*rand(N,1);
    XbHyb = rRF_THzbs.*cos(thetaRF_TH);
    YbHyb = rRF_THzbs.*sin(thetaRF_TH);

    % Distacnces from RF_THz BSs
    [Xmp_Tmat_Integ, Xp_Tmat_Integ] = meshgrid(Xu,XbHyb);
    [Ymp_Tmat_Integ, Yp_Tmat_Integ] = meshgrid(Yu,YbHyb);
    D_ue_RF_THzbs = sqrt((Xmp_Tmat_Integ-Xp_Tmat_Integ).^2 + (Ymp_Tmat_Integ-Yp_Tmat_Integ).^2);
    % Fading channels for RF mode of Hybrid BSs
    fadeRand_Integ = exprnd(1,N,1);


    [Output_SA] = Stand_alone_Rate(D_ue_Rbs,D_ue_Tbs,Params,fadeRand);
    max_rate_SA(i) = Output_SA.max_rate_Upper;
    max_SE_SA(i) = Output_SA.max_SE_Ass_rate;
    RFassoc_SA(i) = Output_SA.RFassoc;
    THassoc_SA(i) = Output_SA.THassoc;
    % Integrated:
    [Output_Int] = Integ_Rate(D_ue_RF_THzbs,Params,fadeRand_Integ);
    max_rate_Int(i) = Output_Int.max_rate_Upper;
    max_SE_Int(i) = Output_Int.max_SE_Ass_rate;
    RFassoc_Int(i) = Output_Int.RFassoc;
    THassoc_Int(i) = Output_Int.THassoc;
end  
    % Stand-alone Average rates:
    Rate_SA = mean(max_rate_SA);
    Sum_rate_SA = 2 * N * mean(max_rate_SA);
    Sum_SE_SA = 2 * N * mean(max_SE_SA);
    SE_SA = mean(max_SE_SA);

    RF_ASSOC_SA = mean(RFassoc_SA);
    THz_ASSOC_SA = mean(THassoc_SA);

    % Integrated Average rates:
    Rate_Intg = mean(max_rate_Int);
    Sum_rate_Intg = 2 * N * mean(max_rate_Int);
    Sum_SE_Intg = 2 * N * mean(max_SE_Int);
    SE_Intg = mean(max_SE_Int);

    RF_ASSOC_Intg = mean(RFassoc_Int);
    THz_ASSOC_Intg = mean(THassoc_Int);
end