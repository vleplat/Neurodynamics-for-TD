function [Y_hat]=     miss_data(Y)
%-------------------------------------------------------------------------%
% [Y_hat]=     miss_data(Y)
%
% Missing data interpolation for three-way tensors. Interpolation is done
% through the third way by replacing 'NaN' by the average of the previous
% and next valus.
%
% The inputs are:
% 
% - Y           : data block.
%
% The outputs are:
% 
% - Y_hat       : interpolated tensor.
%
% List of updates                 -     28/12/2014  -     Rodrigo Cabral
%                                       Creation of the file 
%-------------------------------------------------------------------------%

%-------------------------------Parameters--------------------------------%
% Tensor dimensions
dim =     size(Y);
%-------------------------------------------------------------------------%

%--------------------------Tensor interpolation---------------------------%
% Copying the original tensor
Y_hat     =     Y;
% Counting the number of missing data
n_nan     =     sum(sum(sum(isnan(Y_hat))));
% Interpolate while there are still missing data
while     n_nan >     0
% Subscripts of the missing data
[I,J,K]   =     ind2sub(dim,find(isnan(Y_hat)));
% Interpolation
for i     =     1:length(I)
    % Indexes in the first slice through the third mode
    if          K(i)  ==    1
    Y_hat(I(i),J(i),1)      =     Y_hat(I(i),J(i),2);
    % Indexes in the last slice through the third mode
    elseif      K(i)  ==    dim(3)
    Y_hat(I(i),J(i),dim(3)) =     Y_hat(I(i),J(i),dim(3)-1);
    % Correction for close missing data
    elseif      isnan(Y_hat(I(i),J(i),K(i)-1))
    Y_hat(I(i),J(i),K(i))   =     Y_hat(I(i),J(i),K(i)+1);
    % Correction for close missing data
    elseif      isnan(Y_hat(I(i),J(i),K(i)+1))
    Y_hat(I(i),J(i),K(i))   =     Y_hat(I(i),J(i),K(i)-1);
    % Interpolation
    else
    Y_hat(I(i),J(i),K(i))   =     (Y_hat(I(i),J(i),K(i)-1)+...
                                  Y_hat(I(i),J(i),K(i)+1))/2;
    end
end
% Number of remaining missing data after interpolation
n_nan     =     sum(sum(sum(isnan(Y_hat))));
end
%-------------------------------------------------------------------------%