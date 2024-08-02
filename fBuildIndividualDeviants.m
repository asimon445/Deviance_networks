function [Deviant_conns] = fBuildIndividualDeviants(X1,X2)
%
% fBuildIndividualDeviants will build a deviant network for each individual
% in a sample relative to a normative population. 
%
% Input arguments:
%  'X1': a 3-dimensional matrix of connectivity matrices from the patients
%  'X2': a 3-dimensional matrix of connectivity matrices from a normative
%      sample
%  
% Output arguments:
%  'Deviant_conns': This is a 3D matrix containing deviance connectomes for 
%      all participants 


norm_conn_avg = squeeze(mean(X2,3));
norm_conn_stdev = squeeze(std(X2,0,3));

D = diag(norm_conn_avg);
norm_avg_vec = [squareform((norm_conn_avg-diag(D)).')];
clear D

D = diag(norm_conn_stdev);
norm_std_vec = [squareform((norm_conn_stdev-diag(D)).')];
clear D

% figure out how many z-scores from the norm each patient's connectome is 
for s = 1:size(X1,3)

    conn = squeeze(X1(:,:,s));
    D = diag(conn);
    vec = [squareform((conn-diag(D)).')];

    for e = 1:length(vec)
        zvec(1,e) = (vec(e) - norm_avg_vec(e)) / norm_std_vec(e);
    end

    conn_no_diag = squareform(zvec(1,:));

    % add the diagonal elements back to the matrix
    conn_reconstructed = conn_no_diag + diag(D);

    % Since we removed the diagonal initially, we ensure symmetry by:
    Deviant_conns(:,:,s) = triu(conn_reconstructed) + triu(conn_reconstructed, 1)';

    clear conn D vec zvec
end

end