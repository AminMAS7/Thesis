function [w_I,w_e,f_I,f_e,K_f] = EdgeBands(Tolerance)
%Tolerance = 0.01;
BW = 0.76e12:1e9:0.915e12; % Range of THz bandwidth
x = BW./1e12;
b0 = 0.5;
b1 = 9.6221;
b2 = -8.1526;
b3 = 0.0139;
%y = b0*exp(-1./((b1*x+b2).^2))+b3;
f_I = BW(1)./1e12; % Starting edge band
f_e = BW(end)./1e12; % End edge band
% Solving for feasibility problem
w_I = optimvar('w_I');
w_e = optimvar('w_e');

prob = optimproblem;
cons1 = sqrt(((b0*exp(-1./((b1*(f_e-w_e)+b2).^2))+b3) - (b0*exp(-1./((b1*(f_I+w_I)+b2).^2))+b3))^2) <= Tolerance;
cons2 = w_e >= 0;
cons3 = w_I >= 0;
prob.Constraints.cons1 = cons1;
prob.Constraints.cons2 = cons2;
prob.Constraints.cons3 = cons3;

x0.w_I = 0.00001;
x0.w_e = 0.00001;
opts=optimoptions(@fmincon,'Display','off');
[sol,~,~,~] = solve(prob,x0,'Options',opts);



C25a = (b0*exp(-1./((b1*(f_e-sol.w_e)+b2).^2))+b3) - (b0*exp(-1./((b1*(f_I+sol.w_I)+b2).^2))+b3) <= Tolerance;
C25b = (b0*exp(-1./((b1*(f_I+sol.w_I)+b2).^2))+b3) - (b0*exp(-1./((b1*(f_e-sol.w_e)+b2).^2))+b3) <= Tolerance;


Mid_freq = ((f_e-sol.w_e)+(f_I+sol.w_I))/2;
K_f = b0*exp(-1./((b1*Mid_freq+b2).^2))+b3;
f_I = f_I *1e12;
f_e = f_e *1e12;
w_I = sol.w_I *1e12;
w_e = sol.w_e *1e12;

end