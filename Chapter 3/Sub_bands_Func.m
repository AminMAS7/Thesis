function [S_hat_L,S_hat_U,w_T] = Sub_bands_Func(w_I,w_e,f_I,f_e,Fractional_BW,betta,Guard_Band)
w_g = Guard_Band;
B_Th = Fractional_BW;
% Total bandwidth:
w_T = f_e - f_I; 
% Lower bound of S:
S_hat_L = ceil((2*(w_T-w_I+w_g-w_e)+B_Th*(-w_T+w_I-w_g+w_e))/(2*w_g+2*B_Th*(f_I+w_I)+B_Th*w_g));
% Upper bound of S:
S_hat_U = ceil((w_T-w_I+w_g-w_e)/((1+betta)*w_g))-1;
end