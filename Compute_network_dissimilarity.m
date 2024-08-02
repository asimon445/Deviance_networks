% This script will calculate dissimilarity metrics from Kraus, ..., Gratton
% et al (2021) for measure in the transdiagnositic dataset compared to the
% deviance networks. Then, it will correlate each subject's (dis)similarity
% to the deviance network to the symptom score.
clear;

path = '/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/Full_prediction_results_new/';
files = dir([path '*.mat']); 

load('/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/InputData/connmats_235.mat');
load('/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/InputData/data_full.mat');
load('/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/InputData/header_full.mat');

sim_all_high = zeros(2,length(header));
sim_pos_high = zeros(2,length(header));
sim_neg_high = zeros(2,length(header));

sim_all_low = zeros(2,length(header));
sim_pos_low = zeros(2,length(header));
sim_neg_low = zeros(2,length(header));

% loop through symptoms
for f = 1:length(files)

    load([files(f).folder '/' files(f).name]);

    numnets = length(Deviance_network.target_n_edges);   %index the number of different sized deviance networks that were computed

    % find the column in the data that contains this variable 
    ishigh = strfind(files(f).name,'highSymp');
    islow = strfind(files(f).name,'lowSymp');

    if ~isempty(ishigh)
        index = find(strcmp(files(f).name(1:end-20), header));
    elseif ~isempty(islow)
        index = find(strcmp(files(f).name(1:end-19), header));
    end

    [~, col] = ind2sub(size(header), index);
    clear index

    % Use the deviance networks mask to pull the deviance network
    % Loop through network size
    for n = 1:numnets

        % Turn masks into vectors
        posmask = [squareform(squeeze(Deviance_network.pos_masks(:,:,n)).')];
        negmask = [squareform(squeeze(Deviance_network.neg_masks(:,:,n)).')];
        allmask = [squareform(squeeze(Deviance_network.all_masks(:,:,n)).')];

        % loop through subjects
        ix=0;
        for s = 1:size(data,1)
            if ismember(s,Deviance_network.highest_symp_ix)
                ix=ix+1;
                
                % pull the connectivity from the deviance network masks for
                % the subjects that composed the networks
                thismat = squeeze(mats(:,:,s));
                for d = 1:268
                    thismat(d,d)=0;
                end
                thissub = squareform(thismat);
                
                ixp = 0;
                ixn = 0;
                ixa = 0;
                
                for e = 1:length(thissub)
                    if posmask(e)==1
                        ixp = ixp+1;
                        ixa = ixa+1;
                        devnet_pos(ix,ixp) = thissub(e);
                        devnet_all(ix,ixa) = thissub(e);
                    elseif negmask(e) == 1
                        ixn = ixn+1;
                        ixa = ixa+1;
                        devnet_neg(ix,ixn) = thissub(e);
                        devnet_all(ix,ixa) = thissub(e);
                    end
                end
                clear thissub thismat
            end
        end   % loop through subjects

        devnet_all = mean(devnet_all,1);
        devnet_pos = mean(devnet_pos,1);
        devnet_neg = mean(devnet_neg,1);

        % loop through subjects, and compute (dis)similarity to deviance
        % networks for each subject that did not make up the deviance
        % networks

        ix=0;
        for s = 1:size(data,1)
            if ~ismember(s,Deviance_network.highest_symp_ix)
                ix=ix+1;

                % pull the connectivity from the deviance network masks for
                % the subjects that composed the networks
                thismat = squeeze(mats(:,:,s));
                for d = 1:268
                    thismat(d,d)=0;
                end

                thissub = squareform(thismat);
                
                ixp = 0;
                ixn = 0;
                ixa = 0;
                
                for e = 1:length(thissub)
                    if posmask(e)==1
                        ixp = ixp+1;
                        ixa = ixa+1;
                        net_pos(ix,ixp) = thissub(e);
                        net_all(ix,ixa) = thissub(e);
                    elseif negmask(e) == 1
                        ixn = ixn+1;
                        ixa = ixa+1;
                        net_neg(ix,ixn) = thissub(e);
                        net_all(ix,ixa) = thissub(e);
                    end
                end
                clear thissub thismat

                % compute (dis)similarity
                r_all(ix) = corr(net_all(ix,:)',devnet_all(1,:)','rows','complete','type','Spearman');
                r_pos(ix) = corr(net_pos(ix,:)',devnet_pos(1,:)','rows','complete','type','Spearman');
                r_neg(ix) = corr(net_neg(ix,:)',devnet_neg(1,:)','rows','complete','type','Spearman');

                % pull behavior for this sub
                beh(ix) = cell2mat(data(ix,col));
            end
        end   % loop through subjects

        % correlate (dis)similarity with behavior
        if ~isempty(ishigh)
            sim_all_high(n,col) = corr(beh',r_all','rows','complete','type','Spearman');
            sim_pos_high(n,col) = corr(beh',r_pos','rows','complete','type','Spearman');
            sim_neg_high(n,col) = corr(beh',r_neg','rows','complete','type','Spearman');
        elseif ~isempty(islow)
            sim_all_low(n,col) = corr(beh',r_all','rows','complete','type','Spearman');
            sim_pos_low(n,col) = corr(beh',r_pos','rows','complete','type','Spearman');
            sim_neg_low(n,col) = corr(beh',r_neg','rows','complete','type','Spearman');
        end
        clear r_* beh devnet_* clear net_*
    end   % loop through network sizes
    clear col
end   % loop through symptoms



%% manual for now -- fix this -- it's for plotting

ix=0;
for i = 49:80
    ix=ix+1;
    temp{1,ix}=header{1,i};
end

figure;
bar(sim_all_low(1,49:80))
set(gca,'XTickLabel',temp,'TickLabelInterpreter','none');
ylabel('Correlation coefficient (rho)','FontSize',14)
title('Correlations between symtoms and deviance network (dis)similarity (1000 edge deviance network)','FontSize',14)
ylim([-0.2 0.2])
xticks(1:33)


figure;
bar(sim_all_low(2,49:80))
set(gca,'XTickLabel',temp,'TickLabelInterpreter','none');
ylabel('Correlation coefficient (rho)','FontSize',14)
title('Correlations between symtoms and deviance network (dis)similarity (2000 edge deviance network)','FontSize',14)
ylim([-0.2 0.2])
xticks(1:33)


