function [ theta ] = ParFold( A,B,C )
%Folds A B and C to compute parameters vector theta

[R1,R]    =     size(A);
[R2,R]    =     size(B);
[R3,R]    =     size(C);

theta     =     reshape([A;B;C],R*(R1+R2+R3),1);

end