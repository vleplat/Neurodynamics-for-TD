function [alphak,prod] = step_CubReg_Newton(H,g)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Hinv = inv(H);

% g_norm = sqrt((g)'*(g/H));
prod = g/H;
g_norm = sqrt(trace(g'*(prod)));

Lf = max(eig(H));
alphak = (-1+sqrt(1+2*Lf*g_norm))/(Lf*g_norm);
end