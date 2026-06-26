function [w,F_S_vec] = BW_Freq_Vec_Func(S,f_I,w_T,w_e,w_I,w_g)
w = (w_T - w_e - w_I - (S - 1)*w_g)/S;
F_S_vec = zeros(S,1);
for s = 1 : S
    F_S_vec(s) = f_I + w_I + (s-0.5)*w + (s-1)*w_g;
end
end