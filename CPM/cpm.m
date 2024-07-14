classdef cpm < predictory
    methods
        function this = cpm(group,options)
            this = this@predictory(group,options);
        end
        function run(this)
            rng(this.seed);
            all_edges = this.group.all_edges;
            all_edges = permute(all_edges, [1, 3, 2]);
            all_edges = reshape(all_edges, [],this.num_sub_total);
            
            indices = cvpartition(this.num_sub_total,'k',this.k);
            for i_fold = 1 : this.k
               % fprintf('%dth fold\n', i_fold);
                test.indx = indices.test(i_fold);
                train.indx = indices.training(i_fold);
                test.x = all_edges(:,test.indx);
                train.x = all_edges(:,train.indx);
                test.y = this.phenotype.all_behav(indices.test(i_fold),:);
                train.y = this.phenotype.all_behav(indices.training(i_fold),:);
                
                % first step univariate edge selection
                if size(this.control,1)== size(this.phenotype.all_behav, 1)
                    train.control = this.control(indices.training(i_fold),:);
                    [edge_corr, edge_p] = partialcorr(train.x', train.y, train.control);
                else
                    [edge_corr, edge_p] = corr(train.x', train.y,'rows','complete');
                end
                edges_pos = (edge_p < this.thresh) & (edge_corr > 0);
                edges_neg = (edge_p < this.thresh) & (edge_corr < 0);
                
                this.all_pos_edges(:,i_fold) = edges_pos;
                this.all_neg_edges(:,i_fold) = edges_neg;

                this.edges_pos(:,i_fold) = edge_corr(:,1);
                this.edges_neg(:,i_fold) = edge_corr(:,1);
                
                % build model on TRAIN subs
                train_sum = (nansum(train.x(edges_pos, :), 1) - nansum(train.x(edges_neg, :), 1))';
                fit_train = polyfit(train_sum, train.y(:,1), 1);
                %fit_train = polyfit(train_sum, train.y, 1);
                
                % run model on TEST sub
                test_sum = sum(test.x(edges_pos, :), 1) - sum(test.x(edges_neg,:), 1);
                
                this.Y(test.indx) = (test_sum*fit_train(1)+fit_train(2))';

                % build and test a model for each edge that makes the p-cutoff
                % then save the intercept and slope for each edge -- if an
                % edge doesnt make the p-cutoff, then store that slope/int
                % as a NaN for this iteration


                

%                 % store model info (slope and intercept) and training data
%                 % for each iteration -- only works for LOOCV (need to
%                 edit as of 1/10/23)
%                 this.modelSlope(i_fold,1) = fit_train(1);
%                 this.modelInt(i_fold,1) = fit_train(2);
%                 this.trains{i_fold} = train_sum;
%                 this.testsum(i_fold,1) = test_sum;
            end
        end
        function this = evaluate(this)
            [this.r_pearson, this.p_pearson] = corr(this.Y, this.phenotype.all_behav,'rows','complete');
            [this.r_rank, this.p_rank] = corr(this.Y, this.phenotype.all_behav, 'type', 'spearman','rows','complete');
            this.mse = sum((this.Y - this.phenotype.all_behav).^2) / this.num_sub_total;
            this.q_s = 1 - this.mse / var(this.phenotype.all_behav, 1,'omitnan');
%             fprintf('q_s=%f\n',this.q_s);
%             fprintf('spearman=%f\n',this.r_rank);
%             fprintf('spearman pval=%.30f\n',this.p_rank);

        end
    end
end
