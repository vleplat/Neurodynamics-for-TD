function [ A,B,C ] = ParUnfold( theta,R1,R2,R3,R )
%Unfolds theta to compute factors A B C

M   =     reshape(theta,R1+R2+R3,R);

A   =     M(1:R1,:);
B   =     M(1+R1:R1+R2,:);
C   =     M(1+R1+R2:R1+R2+R3,:);


end



