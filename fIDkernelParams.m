function [bestpredictions,ptype,plambda] = fIDkernelParams(predictions,types,lambdas)
% fIDkernelParams identifies the parameters in kernel ridge CPM that led to
% the best prediction
%   
% Input arguments:
%    'predictions': is a struct containing 3 fields (one for each
%       prediction type). Each field is an MxN matrix, where M is the number
%       of permutations and N is the number of lambdas used
%    'types': is a cell array indicating the types of kernels used
%    'lambdas': a vector of the lambdas used in each type of kernel
%       ridge CPM
%
% Output arguments:
%    'bestpredictions' will be a vector containing the rho values of the
%        type/lambda with the strongest predictions
%    'ptype' is the kernel type that led to the strongest prediction
%    'plambda' is the lambda that led to the strongest prediction

nlambda = length(lambdas);

fnames = fieldnames(predictions);

ix=0;
for f = 1:length(fnames)
    for l = 1:nlambda
        ix=ix+1;
        medians(1,ix) = eval(sprintf('median(predictions.%s(:,l))',fnames{f}));
    end
end

maxmed = find(medians == max(medians));

% only use the first lambda if the strongest prediction occurs at
% multiple lambdas/types
if size(maxmed,2) > 1
    maxmed = maxmed(1,1);
end

if maxmed <= 22
    kr_vec = eval(sprintf('predictions.%s(:,maxmed)',fnames{1}));
    ptype = types{1,1};
    plambda = lambdas(1,maxmed);
elseif maxmed > 22 && maxmed <= 44
    kr_vec = eval(sprintf('predictions.%s(:,maxmed-nlambda)',fnames{2}));
    ptype = types{1,2};
    plambda = lambdas(1,maxmed-nlambda);
elseif maxmed > 44
    kr_vec = eval(sprintf('predictions.%s(:,maxmed-(nlambda*2))',fnames{3}));
    ptype = types{1,3};
    plambda = lambdas(1,maxmed-(nlambda*2));
end

bestpredictions(:,1) = kr_vec(:,1);

end