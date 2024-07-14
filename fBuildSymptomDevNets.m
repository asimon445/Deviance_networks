function [Deviance_network] = fBuildSymptomDevNets(X1,X2,Y,direction,target_nedges,nDevs)
% 
% fBuildSymptomDeviants will build a deviant network based on the symptom
% you are trying to model. For example, if you are trying to build a
% deviant network of the most depressed patients in a sample, then this
% will identify the patients with the worst symptoms, and find the networks
% that they all share that are abnormal relative to a normative sample (e.g., 
% HPC healthy aging) 
%
% Input arguments:
%  'X1': a 3-dimensional matrix of connectivity matrices from the patients
%  'X2': a 3-dimensional matrix of connectivity matrices from a normative
%      sample
%  'Y': the clinical feature of interest (e.g., depression symptom scores)
%  'Direction': specifies whether you want to build a deviant network based
%      on the subjects with the highest (1) or lowest (0) symptom scores. 
%      High = 1, low = 0.
%  'target_nedges': the target number of edges to include in the deviant
%      network mask
%  'nDevs': the number of patients used to construct the deviance networks
%  
% Output arguments:
%  'Deviant_network': This is a struct containing network masks for
%      positive and negative deviant edges that range from 500-2000 edges in
%      500 edge increments

z_thresh = 0.5:0.1:5;
N = size(X1,3);

[feats_ascending,ix_feats] = sort(Y);

% high symptoms = deviant network
if direction == 1
    n_high_symp_subs = length(find(feats_ascending == max(feats_ascending)));   % the n subjects with the worst symptom scores will be compared to the normative model
    
    % build a deviant network on a minimum of 5 subjects
    if n_high_symp_subs < nDevs
        m = 0;
        ix=1;
        while m==0
            ix=ix+1;
            maxn = maxk(feats_ascending,ix);
            n_high_symp_subs = length(find(feats_ascending >= maxn(end)));
    
            if n_high_symp_subs >= nDevs
                m=1;
            end
        end
        clear maxn ix m
    end
    
    highest_symp_ix = ix_feats(N-n_high_symp_subs+1:N,1);

% Low symptoms = deviant network
elseif direction == 0
    n_high_symp_subs = length(find(feats_ascending == min(feats_ascending)));   % the n subjects with the worst symptom scores will be compared to the normative model
    
    % build a deviant network on a minimum of 5 subjects
    if n_high_symp_subs <= 4
        m = 0;
        ix=1;
        while m==0
            ix=ix+1;
            maxn = maxk(feats_ascending,ix);
            n_high_symp_subs = length(find(feats_ascending >= maxn(end)));
    
            if n_high_symp_subs > 4
                m=1;
            end
        end
        clear maxn ix m
    end
    
    highest_symp_ix = ix_feats(1:n_high_symp_subs,1);
else
    error('You can only specifiy a 1 (high symptom levels) or a 0 (low symptoms levels) for the ''direction'' parameter');
end

norm_conn_avg = squeeze(mean(X2,3));
norm_conn_stdev = squeeze(std(X2,0,3));

for i = 1:n_high_symp_subs
    high_symp_conn(:,:,i) = X1(:,:,highest_symp_ix(i,1));
end

high_symp_conn_avg = squeeze(mean(high_symp_conn,3));
high_symp_conn_std = squeeze(std(high_symp_conn,0,3));

D = diag(high_symp_conn_avg);
high_symp_vec = [squareform((high_symp_conn_avg-diag(D)).')];
clear D

D = diag(high_symp_conn_std);
high_symp_std = [squareform((high_symp_conn_std-diag(D)).')];
clear D

D = diag(norm_conn_avg);
norm_avg_vec = [squareform((norm_conn_avg-diag(D)).')];
clear D

D = diag(norm_conn_stdev);
norm_std_vec = [squareform((norm_conn_stdev-diag(D)).')];
clear D

% find the right z-score for building the deviant network
for st = 1:length(z_thresh)
    deviants_pos = zeros(length(high_symp_vec),1);
    deviants_neg = zeros(length(high_symp_vec),1);

    for e = 1:length(high_symp_vec)
        if high_symp_vec(1,e) > (norm_avg_vec(1,e)+(norm_std_vec(1,e)*z_thresh(st)))
            deviants_pos(e,1) = 1;
        elseif high_symp_vec(1,e) < (norm_avg_vec(1,e)-(norm_std_vec(1,e)*z_thresh(st)))
            deviants_neg(e,1) = 1;
        end
    end

    deviants_pos = logical(deviants_pos);
    deviants_neg = logical(deviants_neg);

    npos(st) = sum(deviants_pos);
    nneg(st) = sum(deviants_neg);

    clear deviants_neg deviants_pos
end

alldevs = npos+nneg;

Deviance_network.num_subjects = n_high_symp_subs;
Deviance_network.highest_symp_ix = highest_symp_ix;

% Find numbers greater than target number of edeges
for te = 1:length(target_nedges)
    over = alldevs(alldevs > target_nedges(te));

    % Find the closest number
    closest_to_targ = find(min(over) == over);

    % build the deviant network
    deviants_pos = zeros(length(high_symp_vec),1);
    deviants_neg = zeros(length(high_symp_vec),1);

    deviants_pos_var = zeros(length(high_symp_vec),1);
    deviants_neg_var = zeros(length(high_symp_vec),1);

    for e = 1:length(high_symp_vec)
        if high_symp_vec(1,e) > (norm_avg_vec(1,e)+(norm_std_vec(1,e)*z_thresh(closest_to_targ)))
            deviants_pos(e,1) = 1;
            deviants_pos_var(e,1) = high_symp_std(1,e);
        elseif high_symp_vec(1,e) < (norm_avg_vec(1,e)-(norm_std_vec(1,e)*z_thresh(closest_to_targ)))
            deviants_neg(e,1) = 1;
            deviants_neg_var(e,1) = high_symp_std(1,e);
        end
    end

    n_dev_edges_temp = sum(deviants_pos) + sum(deviants_neg);
    Deviance_network.z_threshold(te) = z_thresh(closest_to_targ);

    clear over closest_to_targ

    vlen = length(deviants_pos);
    smat = (1 + sqrt(1 + 8*vlen)) / 2; % Size of the original matrix
    M_p = zeros(smat); % Initialize matrix with zeros
    M_n = zeros(smat);

    M_pv = zeros(smat); % Initialize matrix with zeros
    M_nv = zeros(smat);

    % Reconstruct matrix of the pos deviant network
    M_p(tril(true(smat), -1)) = deviants_pos;
    M_pv(tril(true(smat), -1)) = deviants_pos_var;

    % Transpose to get upper triangle
    M_p = M_p + M_p.'; % This ensures symmetry
    M_pv = M_pv + M_pv.'; % This ensures symmetry

    % Restore diagonal
    M_p(1:N+1:end) = 0;
    M_pv(1:N+1:end) = 0;

    % Reconstruct lower matrix of the pos deviant network
    M_n(tril(true(smat), -1)) = deviants_neg;
    M_nv(tril(true(smat), -1)) = deviants_neg_var;

    % Transpose to get upper triangle
    M_n = M_n + M_n.'; % This ensures symmetry
    M_nv = M_nv + M_nv.'; % This ensures symmetry

    % Restore diagonal
    M_n(1:N+1:end) = 0;
    M_nv(1:N+1:end) = 0;

    Deviance_network.pos_masks(:,:,te) = M_p;
    Deviance_network.neg_masks(:,:,te) = M_n;
    Deviance_network.all_masks(:,:,te) = M_p + M_n;

    Deviance_network.pos_net_variance(:,:,te) = M_pv;
    Deviance_network.neg_net_variance(:,:,te) = M_nv;
    Deviance_network.comb_net_variance(:,:,te) = M_pv + M_nv;

    Deviance_network.target_n_edges(te) = target_nedges(te);

    clear n_dev_edges_temp deviants_pos deviants_neg M_p M_n vlen smat
end

end