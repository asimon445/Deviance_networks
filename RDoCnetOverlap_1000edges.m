% This script identifies the symptoms that we can predict well using the
% deviance network modeling, and identifies the contribution of cognitive
% networks to maintaining symptomatology by removing edges from the 
% cognitive networks that are present in the deviance networks and
% re-running predictions. If deviance predictions are worse without the
% construct network edges present, then that suggests that that cognitive 
% construct is important for maintaining the symptoms.
%
% Then, for each symptom that we can predict, identify the cognitive
% measures that it is correlated with. Since we can predict cognitive
% measures using the construct networks, we will then remove the edges from
% the deviance networks that are embedded in the construct networks and run
% those predictions again. If the prediction strength drops in the same 
% construct networks that we identified are important in psych symptoms, 
% then that suggests that the relationship between symptoms and cognition 
% is supported by the deviance edges that are embedded in those construct
% networks.


clear;

path = '/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/Deviance_prediction_results/HighScoresModeled/Symp/';
outdir_symps = '/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/lesion_analysis/symptom_predictions_RDoC_lesions_1000edges/';
outdir_cog = '/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/lesion_analysis/cognitive_predictions_deviance_lesions_1000edges/';

load('/Users/ajsimon/Documents/Data/Constable_lab/Network_compensation/Input_data/RDoC_construct_networks.mat');
load('/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/InputData/connmats_235.mat');
load('/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/InputData/header_full.mat');
load('/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/InputData/data_full.mat');

files = dir([path '*_highSymp_devnet.mat']);

cutoff = 0.18;   % Don't look at symptoms with deviance predictions lower than this coefficient. If it's really weak, it won't be interesting.
nperms = 100;
N = size(mats,3);
types = {'corr','exp','gaus'};
lambdas = [0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];   % DO NOT MODIFY THIS!

% turn cog construct networks into vectors
for r = 1:268
    RDoC_masks_forCPM(r,r,:) = 0;
end

for n = 1:size(RDoC_masks_forCPM,3)
    RDoCnet(:,n) = [squareform(squeeze(RDoC_masks_forCPM(:,:,n)).')]';
end
clear RDoC_masks_forCPM

%% iterate through symptoms that have been modeled using the deviance network approach
ix=0;
for f = 24  %1:length(files)
    
    israw = strfind(files(f).name,'_raw');   % skip raw scores

    if isempty(israw)
        
        load([files(f).folder '/' files(f).name]);

        symp_ix = 0;
        for h = 1:length(header)
            try
                thisone = strcmp(files(f).name(1:length(header{1,h})),header{1,h});
                if thisone == 1
                    symp_ix = h;
                    break
                end
            end
        end

        try
            for i = 1:size(data,1)
                if ~isempty(strfind(data{i,symp_ix},'na'))
                    data{i,symp_ix} = NaN;
                elseif isempty(data{i,symp_ix})
                    data{i,symp_ix}=NaN;
                end
                feat_of_interest(i,1) = cell2mat(data(i,symp_ix));
            end

        catch
            error('Failed to find the index of the symptom this deviance network is modeling');
        end

        predStren_temp = Deviance_network.Deviance_preds_krCPM(:,1);   %vector of prediction strength values from 2000 edge deviance krCPM
    
        % Don't look at symptoms with really weak deviance predictions
        if median(predStren_temp) > cutoff
            ix=ix+1;

            predStren = predStren_temp;
            clear predStren_temp

            % convert the deviance network mask into a vector
            devnet = [squareform(squeeze(Deviance_network.all_masks(:,:,1)).')]';
        
            % structure data for deviants predictions
            ixx=0;
            for i = 1:N
                if sum(i ~= Deviance_network.highest_symp_ix) == Deviance_network.num_subjects && ~isnan(feat_of_interest(i,1))
                    ixx=ixx+1;
                    x(:,:,ixx) = mats(:,:,i);
                    y(ixx,1) = feat_of_interest(i,1);
                end
            end

            RDoCnet_lesioned = logical(RDoCnet);
            anyRDoCnet = sum(RDoCnet,2);

            % replace edges shared between networks with a 1
            for e = 1:length(anyRDoCnet)
                if anyRDoCnet(e) > 1
                    anyRDoCnet(e) = 1;
                end
            end

            % Find positions where devnet contains a 1, but anyRDoCnet contains a 0
            exclusiveDevnet_ix = find(devnet == 1 & anyRDoCnet == 0);
            exclusiveDevnet = zeros(length(devnet),1);
            exclusiveDevnet(exclusiveDevnet_ix) = 1;
            n_exclusiveDevnet = length(exclusiveDevnet_ix);

            for n = 1:size(RDoCnet,2)
                idx=0;

                shared_ix = find(devnet == 1 & RDoCnet(:,n) == 1);
                shared = zeros(length(devnet),1);
                shared(shared_ix) = 1;

                % remove the contribution from this construct network in
                % the deviance network prediction
                devnet_lesioned = logical(devnet);
                devnet_lesioned(shared_ix) = false;

                % lesion out shared edges from the construct networks to
                % use later
                for e = 1:size(RDoCnet,1)
                    if RDoCnet(e,n) == 1 && devnet(e) == 1
                        idx=idx+1;
                        RDoCnet_lesioned(e,n) = false;
                    end
                end

                % count number of edges that overlap between cognitive
                % construct networks and deviance networks
                percShared_edges(ix,n) = idx / sum(RDoCnet(:,n));
                clear idx

                % count number of shared edges between the deviance network
                % and the construct network -- this is how many we'll
                % lesion out from the deviance network with all shared
                % construct network edges removed
                ixshared = find(RDoCnet(:,n) == 1 & devnet == 1);
                nShared_edges(1,n) = length(ixshared);

                % run kr prediction on the deviance network with the RDoC
                % network edges lesioned out
                for np = 1:nperms
                    if np == 1
                        fprintf('Predicting %s from deviance networks with %s edges lesioned out \n',header{1,h},labels{n});
                        fprintf('%s network lesions: ',labels{n});
                    end
                    fprintf('%d ', np);

                    % build null model
                    randEdges = exclusiveDevnet_ix(randperm(n_exclusiveDevnet, nShared_edges(1,n)));
                    devnet_null_lesioned = logical(devnet);

                    % make sure that these edges ARE NOT in the RDoC networks
                    for rn = 1:length(randEdges)
                        devnet_null_lesioned(randEdges(rn)) = false;
                    end

                    test = main_forceEdges_krCPM(x,y,10,devnet_lesioned,devnet_lesioned,'results',1);
                    null = main_forceEdges_krCPM(x,y,10,devnet_null_lesioned,devnet_null_lesioned,'results',1); 

                    testpreds.kr_rho_corr(np,:) = test.r_rank_corr;
                    testpreds.kr_rho_exp(np,:) = test.r_rank_Exponential;
                    testpreds.kr_rho_gaus(np,:) = test.r_rank_Gaussian;

                    nullpreds.kr_rho_corr(np,:) = null.r_rank_corr;
                    nullpreds.kr_rho_exp(np,:) = null.r_rank_Exponential;
                    nullpreds.kr_rho_gaus(np,:) = null.r_rank_Gaussian;

                    clear test null devnet_null_lesioned randEdges
                end
                fprintf('\n');
                
                [a,b,c] = fIDkernelParams(testpreds,types,lambdas);
                Dev_RDoClesions_preds(:,n) = a;
                Dev_RDoClesions_predType{n} = b;
                Dev_RDoClesions_predLambda(n) = c;

                clear a b c

                [a,b,c] = fIDkernelParams(nullpreds,types,lambdas);
                Dev_RDoClesions_nulls(:,n) = a;
                Dev_RDoClesions_nulls_predType{n} = b;
                Dev_RDoClesions_nulls_predLambda(n) = c;

                shared_edges(:,n) = shared;

                clear a b c testpreds nullpreds devnet_lesioned rand_devix randomN devnetSize dev_ix shared shared_ix ixshared 
            end 
            clear RDoC_masks_forCPM x exclusiveDevnet_ix exclusiveDevnet y n_exclusiveDevnet

            % Save prediction results
            Deviance_lesions.lesion_results = Dev_RDoClesions_preds;
            Deviance_lesions.lesion_lambda = Dev_RDoClesions_predLambda;
            Deviance_lesions.lesion_type = Dev_RDoClesions_predLambda;
            Deviance_lesions.null_results = Dev_RDoClesions_nulls;
            Deviance_lesions.null_lambda = Dev_RDoClesions_nulls_predLambda;
            Deviance_lesions.null_type = Dev_RDoClesions_nulls_predLambda;

            Deviance_lesions.deviance_cogconstruct_shared_edges = shared_edges;
            Deviance_lesions.num_shared_edges = nShared_edges;
            Deviance_lesions.netlabels = labels;

            clear Dev_RDoClesions_* shared_edges 

            outfile = sprintf('%s/%s_deviance_preds_with_construct_nets_lesioned.mat',outdir_symps,header{1,h});
            save(outfile,'Deviance_lesions');

            clear Deviance_lesions

            %% identify which cog measures this symptom is correlated with
            idx=0;
            for cg = 88:142
                israw = strfind(header{1,cg},'_raw');   % skip raw scores
                if isempty(israw)

                    idx=idx+1;

                    % loop through subs, make sure its not a deviant sub
                    % and make sure it's not a nan. If it's a nan, remove
                    % the same position in 'y'
                    for i = 1:size(data,1)
                        if ~isempty(strfind(data{i,cg},'na'))
                            data{i,cg} = NaN;
                        elseif isempty(data{i,cg})
                            data{i,cg}=NaN;
                        end

                        cogofinterest(i,1) = cell2mat(data(i,cg));
                    end

                    cogmeasures(:,idx) = cogofinterest;

                    ixx = 0;
                    for i = 1:size(data,1)
                        if sum(i ~= Deviance_network.highest_symp_ix) == Deviance_network.num_subjects && ~isnan(cogofinterest(i,1)) && ~isnan(feat_of_interest(i,1))      
                            ixx=ixx+1;
                            symp(ixx,1) = feat_of_interest(i,1);
                            cogmeasure(ixx,1) = cogofinterest(i,1);
                        end
                    end

                    [symp_cog_corr_rho(1,idx),symp_cog_corr_p(1,idx)] = corr(symp,cogmeasure,'rows','complete','type','Spearman');
                    cogmeas_lbl_all{1,idx} = header{1,cg};

                    clear symp cogmeasure cogofinterest

                end   % if isempty(raw)
            end   % for cg = 88:142

            sig_corrs = find(symp_cog_corr_p(1,:) < ( 0.05 / 55 ));   % 55 is the number of cog scores (excluding raw ones)
            beep;
   
            % if there are not any significant symptom-cognition
            % correlations, then skip this one
            if ~isempty(sig_corrs)
                relevant_corrs_rho = symp_cog_corr_rho(1,sig_corrs);
                relevant_corrs_p = symp_cog_corr_p(1,sig_corrs);

                for lb = 1:length(sig_corrs)
                    relevant_cogmeas_lbl{1,lb} = cogmeas_lbl_all{1,sig_corrs(lb)};
                end

                % Run predictions on each RDoC cognitive construct network 
                % with shared deviance edges lesioned out
                for cg = 1:length(sig_corrs)

                    fprintf('%s and %s are correlated \n',header{1,h},relevant_cogmeas_lbl{1,cg});

                    % loop through subs, make sure its not a deviant sub
                    % and make sure it's not a nan. If it's a nan, remove
                    % the same position in 'y'
                    cogofinterest = cogmeasures(:,sig_corrs(cg));

                    ixx = 0;
                    for i = 1:size(data,1)
                        if sum(i ~= Deviance_network.highest_symp_ix) == Deviance_network.num_subjects && ~isnan(cogofinterest(i,1)) && ~isnan(feat_of_interest(i,1))      
                            ixx=ixx+1;
                            symp(ixx,1) = feat_of_interest(i,1);
                            cogmeasure(ixx,1) = cogofinterest(i,1);
                            x(:,:,ixx) = mats(:,:,i);
                        end
                    end

                    for n = 1:length(labels)
                        if n == 1
                            fprintf('Predicting %s from cognitive construct networks with %s deviance edges lesioned out \n',relevant_cogmeas_lbl{1,cg},header{1,h})
                        end

                        % Find positions where devnet contains a 0, but this RDoC network contains a 1
                        exclusiveCognet_ix = find(devnet == 0 & RDoCnet(:,n) == 1);
                        exclusiveCognet = zeros(length(RDoCnet(:,n)),1);
                        exclusiveCognet(exclusiveCognet_ix) = 1;
                        n_exclusiveCognet = length(exclusiveCognet_ix);

                        shared_ix = find(devnet == 1 & RDoCnet(:,n) == 1);
                        shared = zeros(length(RDoCnet(:,n)),1);
                        shared(shared_ix) = 1;

                        cognet_lesioned = RDoCnet_lesioned(:,n);
                        cognet = logical(RDoCnet(:,n));

                        for np = 1:nperms
                            if np == 1
                                fprintf('%s network predictions: ',labels{n});
                            end
                            fprintf('%d ', np);

                            % build null model
                            randEdges = exclusiveCognet_ix(randperm(n_exclusiveCognet, nShared_edges(1,n)));
                            cognet_null_lesioned = logical(RDoCnet(:,n));

                            % make sure that these edges ARE NOT in the RDoC networks
                            for rn = 1:length(randEdges)
                                cognet_null_lesioned(randEdges(rn)) = false;
                            end

                            cog = main_forceEdges_krCPM(x,cogmeasure,10,cognet,cognet,'results',1);
                            test = main_forceEdges_krCPM(x,cogmeasure,10,cognet_lesioned,cognet_lesioned,'results',1);
                            null = main_forceEdges_krCPM(x,cogmeasure,10,cognet_null_lesioned,cognet_null_lesioned,'results',1);

                            cogpreds.kr_rho_corr(np,:) = cog.r_rank_corr;
                            cogpreds.kr_rho_exp(np,:) = cog.r_rank_Exponential;
                            cogpreds.kr_rho_gaus(np,:) = cog.r_rank_Gaussian;

                            testpreds.kr_rho_corr(np,:) = test.r_rank_corr;
                            testpreds.kr_rho_exp(np,:) = test.r_rank_Exponential;
                            testpreds.kr_rho_gaus(np,:) = test.r_rank_Gaussian;

                            nullpreds.kr_rho_corr(np,:) = null.r_rank_corr;
                            nullpreds.kr_rho_exp(np,:) = null.r_rank_Exponential;
                            nullpreds.kr_rho_gaus(np,:) = null.r_rank_Gaussian;

                            clear test null cognet_null_lesioned randEdges
                        end
                        fprintf('\n');

                        [a,b,c] = fIDkernelParams(cogpreds,types,lambdas);
                        Cognet_lesions.full_results(:,n,cg) = a;
                        Cognet_lesions.full_type{n,cg} = b;
                        Cognet_lesions.full_lambda(n,cg) = c;

                        clear a b c

                        [a,b,c] = fIDkernelParams(testpreds,types,lambdas);
                        Cognet_lesions.lesion_results(:,n,cg) = a;
                        Cognet_lesions.lesion_type{n,cg} = b;
                        Cognet_lesions.lesion_lambda(n,cg) = c;

                        clear a b c

                        [a,b,c] = fIDkernelParams(nullpreds,types,lambdas);
                        Cognet_lesions.null_results(:,n,cg) = a;
                        Cognet_lesions.null_type{n,cg} = b;
                        Cognet_lesions.null_lambda(n,cg) = c;

                        clear a b c cognet cognet_lesioned exclusiveCognet_ix exclusiveCognet n_exclusiveCognet shared_ix shared nullpreds testpreds
                    end

                    clear symp cogmeasure cogofinterest x 
                end    % for cg = 1:length(sig_corrs)
                clear cogmeasures

                Cognet_lesions.netlabels = labels;
                Cognet_lesions.symp_cog_corrcoefs = relevant_corrs_rho;
                Cognet_lesions.symp_cog_corrs_p = relevant_corrs_p;

                % save the name of the cognitive measure for each index in the
                % output
                Cognet_lesions.cogmeasure_name = relevant_cogmeas_lbl;

                outfile = sprintf('%s/%s_deviance_lesioned_cog_meas_preds_from_constructNets.mat',outdir_cog,header{1,h});
                save(outfile,'Cognet_lesions');

                clear Cognet_lesions Dev_RDoClesions_* cogmeas_lbl relevant_corrs_* 
            else
                fprintf('%s is not correlated with any cognitive measures.\n',header{1,h});
            end
            clear sig_corrs symp_cog_corr_* anyRDoCnet cogmeas_lbl_all nShared_edges RDoCnet_lesioned relevant_cogmeas_lbl 
        end   % if ~isempty(symp_cog_corr_p)
    end
    clear Deviance_network feat_of_interest y devnet percShared_edges predStren
end


