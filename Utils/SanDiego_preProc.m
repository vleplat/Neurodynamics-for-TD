function [Xsub] = SanDiego_preproc(flag)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
load('SanDiego.mat'); 
X = A'; clear A; 
% preprocess - remove outliers  
minX = min(X);
t = find( minX < 0); 
X(:,t) = 0; 
nx = sum( X.^2 );
[a,b] = sort(-nx);
X(:,b(1:1000)) = 0; 

% select columns
%sub = [1 40 80 120 158]; 
%sub = [1 30 60 90 120 150];  
sub = FastSepNMF(X',5)
% Preprocess
Xsub = X(sub,:); 
end