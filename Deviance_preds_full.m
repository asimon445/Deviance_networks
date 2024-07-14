% This script will average the connectomes of the subjects with the highest
% symptom scores to construct a deviant network for each symptom, and run
% CPM using that binarized network as a mask
%
% The user should first specify the number of target edges to include in 
% the deviant network. The script will then step through a series of
% z-score cutoffs to identify the appropriate cutoff for building a deviant
% network containing the target number of edges. Finally, this script will
% then run traditional and kernel ridge CPM on the deviant network and on
% the whole brain (as well as a random network for comparison). 

clear;

%% load stuff
load('/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/InputData/connmats_235.mat');
load('/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/InputData/header_full.mat');
load('/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/InputData/data_full.mat');
load('/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/InputData/HCP_avg_conns.mat');  % load HCP data

%% set parameters
z_thresh = 0.5:0.1:5;
target_nedges = [1000 2000];     % step through z_thresh until we find the target number of pos+neg deviance edges

lambdas = [0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];   % DO NOT MODIFY THIS!
nperms = 10;
N = size(mats,3);
devpredictions = 1:142;
nDevs = 5;  % this chooses the number of high symptom subjects to use to make the deviance networks
direction = 1;   % use subjects with high symptom values to do deviance network predictions

types = {'corr','exp','gaus'};

idx=0;

%% loop through symptoms
for t = 1:length(devpredictions)

    fprintf('%s\n',header{1,devpredictions(t)});
    symp_ix = devpredictions(t);
    feat_of_interest = cell2mat(data(:,symp_ix));
    
    Deviance_network = fBuildSymptomDevNets(mats,HCP_conns,feat_of_interest,direction,target_nedges,nDevs);
    
    % Find numbers greater than target number of edeges
    for te = 1:length(target_nedges)
        
        % structure data for deviants CPM
        ix=0;
        for i = 1:N
            if sum(i ~= Deviance_network.highest_symp_ix) == Deviance_network.num_subjects && ~isnan(feat_of_interest(i,1))
                ix=ix+1;
                x(:,:,ix) = mats(:,:,i);
                y(ix,1) = feat_of_interest(i,1);
            end
        end

        for np = 1:nperms

            if np == 1
                fprintf('Iteration ')
            end
            fprintf('%d ', np)

            %% run predictions on deviance networks
            % run krCPM 
            D = diag(squeeze(Deviance_network.pos_masks(:,:,te)));
            posmask = [squareform(squeeze(Deviance_network.pos_masks(:,:,te)).')];
            clear D

            D = diag(squeeze(Deviance_network.neg_masks(:,:,te)));
            negmask = [squareform(squeeze(Deviance_network.neg_masks(:,:,te)).')];
            clear D

            posmask = logical(posmask);
            negmask = logical(negmask);

            deviants_kr = main_forceEdges_krCPM(x,y,10,posmask,negmask,'results',1);

            devnetpreds.kr_rho_corr(np,:) = deviants_kr.r_rank_corr;
            devnetpreds.kr_rho_exp(np,:) = deviants_kr.r_rank_Exponential;
            devnetpreds.kr_rho_gaus(np,:) = deviants_kr.r_rank_Gaussian;

            % run CPM using deviance network masks
            deviants_CPM = main_forceEdges_vanillaCPM(x,y,10,posmask,negmask,'results',1);
            Deviance_network.Deviance_preds_CPM(np,te) = deviants_CPM.r_rank;
            
            clear posmask negmask deviants_CPM deviants_kr

            %% run whole brain predictions
            % re-format x and y to include subjects that made up the
            % deviance networks
            clear x y

            ix=0;
            for i = 1:N
                if ~isnan(feat_of_interest(i,1))
                    ix=ix+1;
                    x(:,:,ix) = mats(:,:,i);
                    y(ix,1) = feat_of_interest(i,1);
                end
            end

            % run krCPM
            pos_bin_vec = true(35778,1);
            neg_bin_vec = true(35778,1);

            krCPMpreds = main_forceEdges_krCPM(x,y,10,pos_bin_vec,neg_bin_vec,'results',1);

            krCPM.rho_corr(np,:) = krCPMpreds.r_rank_corr;
            krCPM.rho_exp(np,:) = krCPMpreds.r_rank_Exponential;
            krCPM.rho_gaus(np,:) = krCPMpreds.r_rank_Gaussian;

            % run regular CPM
            CPMpreds = main(x,y,10,'results');
            Deviance_network.Wholebrain_preds_CPM(np,te) = CPMpreds.r_rank;

            clear krCPMpreds CPMpreds pos_bin_vec neg_bin_vec
            
        end

        fprintf('\n');

        % identify lambda and kernel type that predicts best
        % get the strongest predictions from krCPM on the deviance network
        [a,b,c] = fIDkernelParams(devnetpreds,types,lambdas);
        Deviance_network.Deviance_preds_krCPM(:,te) = a;
        Deviance_network.Dev_krCPM_type{te} = b;
        Deviance_network.Dev_krCPM_lambda(te) = c;

        clear devnetpreds a b c
       
        % get the strongest predictions from krCPM on the whole brain
        [a,b,c] = fIDkernelParams(krCPM,types,lambdas);
        Deviance_network.Wholebrain_preds_krCPM(:,te) = a;
        Deviance_network.Wholebrain_krCPM_type{te} = b;
        Deviance_network.Wholebrain_krCPM_lambda(te) = c;

        clear krCPM a b c

        fprintf('\n')

        clear x y

    end

    clear feat_of_interest n_high_symp_subs CPM_r

    % save the predictions
    outfile = sprintf('/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/Full_prediction_results/%s_highSymp_devnet.mat',...
        header{1,devpredictions(t)});
    save(outfile,'Deviance_network');
    clear Deviance_network outfile
end

