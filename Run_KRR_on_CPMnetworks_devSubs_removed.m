% This script will compute the CPM networks for each symptom in half of 
% the sample, and then use the CPM network as a mask for KRR.

clear;

%% load stuff
load('/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/InputData/connmats_235.mat');
load('/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/InputData/header_full.mat');
load('/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/InputData/data_full.mat');
load('/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/InputData/HCP_avg_conns.mat');  % load HCP data

%% set parameters
outdir = '/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/KRR_on_CPM_nets_preds/';
path = '/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/Deviance_prediction_results/HighScoresModeled/Symp/';
dev_files = dir([path '*_highSymp_devnet.mat']);

lambdas = [0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];   % DO NOT MODIFY THIS!
nperms = 10;
N = size(mats,3);
devpredictions = 85; %1:87;
kfolds = 10;   % number of folds used in defining CPM networks in half of the sample

types = {'corr','exp','gaus'};

idx=0;

%% loop through symptoms
for t = 1:length(devpredictions)

    israw = strfind(header{1,devpredictions(t)},'_raw');   % skip raw scores

    if isempty(israw)
        
        outfile = sprintf('%s/%s_CPM_network_KRRpreds_noDevSubs.mat',outdir,header{1,devpredictions(t)});
    
        if ~exist(outfile,'file')
            fprintf('%s\n',header{1,devpredictions(t)});
            symp_ix = devpredictions(t);
            
            for i = 1:size(data,1)
                if ~isempty(strfind(data{i,symp_ix},'na'))
                    data{i,symp_ix} = NaN;
                elseif isempty(data{i,symp_ix})
                    data{i,symp_ix}=NaN;
                end
                feat_of_interest(i,1) = cell2mat(data(i,symp_ix));
            end

            % find the corresponding deviance file to remove the subjects
            % that made up the deviance network from this prediction
            for i = 1:length(dev_files)
                match = strfind(dev_files(i).name,header{1,devpredictions(t)});

                if ~isempty(match)
                    load([dev_files(i).folder '/' dev_files(i).name]);
                    break
                end
            end
            
            % structure data for deviants predictions
            ixx=0;
            for i = 1:N
                if sum(i ~= Deviance_network.highest_symp_ix) == Deviance_network.num_subjects && ~isnan(feat_of_interest(i,1))
                    ixx=ixx+1;
                    x(:,:,ixx) = mats(:,:,i);
                    y(ixx,1) = feat_of_interest(i,1);
                end
            end

            nsubs = length(y);

            % loop through permutations
            for np = 1:nperms

                if np == 1
                    fprintf('Running KRR on CPM networks for %s \n',header{1,devpredictions(t)});
                end

                % Partition the data into 2 halves
                indices = cvpartition(nsubs,'k',2);

                % Store indexes of the train/test partitions for each
                % permutation
                train_indx(:,np) = indices.training(1);
                test_indx(:,np) = indices.test(1);

                train_x = x(:,:,train_indx(:,np));
                train_y = y(train_indx(:,np),1);   

                CPM_network = main(train_x,train_y,kfolds,'CPM network');

                all_pos_edges(:,1:kfolds) = CPM_network.all_pos_edges;
                all_neg_edges(:,1:kfolds) = CPM_network.all_neg_edges;

                pos_edges_allfolds = sum(all_pos_edges,2);
                neg_edges_allfolds = sum(all_neg_edges,2);

                % identify the right threshold for making the CPM network ~1000 edges
                ix=0;
                for te = 0.1:0.1:1
                    ix=ix+1;

                    foldsThresh(ix) = te;
                    pos_edge_ix = find(pos_edges_allfolds > kfolds*te);
                    neg_edge_ix = find(neg_edges_allfolds > kfolds*te);
                    edges = [pos_edge_ix;neg_edge_ix];

                    numEdges(ix) = length(edges);

                    clear pos_edge_ix neg_edge_ix edges
                end

                [~, ix_thresh] = min(abs(numEdges - 1000));
                fold_thresh_used(np) = foldsThresh(ix_thresh);

                clear ix_thresh CPM_network

                pos_edge_ix = find(pos_edges_allfolds > kfolds*fold_thresh_used(np));
                neg_edge_ix = find(neg_edges_allfolds > kfolds*fold_thresh_used(np));
                edges = [pos_edge_ix;neg_edge_ix];

                % create CPM defined network masks for KRR
                CPM_netmask = zeros(length(pos_edges_allfolds),1);
                CPM_netmask(edges) = 1;
                CPM_netmask = logical(CPM_netmask);
                CPM_network(:,np) = CPM_netmask;

                test_x = x(:,:,test_indx(:,np));
                test_y = y(test_indx(:,np),1);   

                KRR = main_forceEdges_krCPM(test_x,test_y,kfolds,CPM_netmask,CPM_netmask,'results',1);

                preds.kr_rho_corr(np,:) = KRR.r_rank_corr;
                preds.kr_rho_exp(np,:) = KRR.r_rank_Exponential;
                preds.kr_rho_gaus(np,:) = KRR.r_rank_Gaussian;

                fprintf('%d ',np);

                clear train_x train_y indices KRR test_x test_y pos_edge_ix neg_edge_ix edges
            end
            fprintf('\n');

            [a,b,c] = fIDkernelParams(preds,types,lambdas);
            CPMnet_KRRpreds.results = a;
            CPMnet_KRRpreds.type = b;
            CPMnet_KRRpreds.lambda = c;
            CPMnet_KRRpreds.CPM_network_selection_threshold = fold_thresh_used;
            CPMnet_KRRpreds.CPM_netmask = CPM_network;

            save(outfile,'CPMnet_KRRpreds');

            clear a b c fold_thresh_used nsubs CPMnet_KRRpreds preds CPM_netmask outfile symp_ix feat_of_interest x y CPM_network

        end
    end
end

          
