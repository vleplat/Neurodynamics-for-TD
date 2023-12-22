load Cuprite_f970619t01p02_r02_sc03.a.rfl.mat
X(:,:,end-3:end)=[];
imageData = X;
load AvirisCuprite_spec.mat;
imageData(:,:,specBand_removed')=[];
X=imageData(2:2+120-1,2:2+120-1,:); % take the subimage 120x120 pixels
V=double(X);  % convert to double format
save('V.mat','V')