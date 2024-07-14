classdef cpm_forceEdges < predictory
    methods
        function this = cpm_forceEdges(group,options)
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

%                 % first step univariate edge selection
%                 if size(this.control,1)== size(this.phenotype.all_behav, 1)
%                     train.control = this.control(indices.training(i_fold),:);
%                     [edge_corr, edge_p] = partialcorr(train.x', train.y, train.control);
%                 else
%                     [edge_corr, edge_p] = corr(train.x', train.y);
%                 end
%                 edges_pos = (edge_p < this.thresh) & (edge_corr > 0);
%                 edges_neg = (edge_p < this.thresh) & (edge_corr < 0);

                pos = find(this.edges_pos == 1);
                neg = find(this.edges_neg == 1);

                comb = sort([pos,neg]);

                edges_pos = zeros(size(all_edges,1),1);
                edges_neg = zeros(size(all_edges,1),1);
                
                for i = 1:length(pos)
                    posix(i,1) = find(pos(i) == comb);
                end

                for i = 1:length(neg)
                    negix(i,1) = find(neg(i) == comb);
                end

                edges_pos(posix,1) = 1;
                edges_neg(negix,1) = 1;

                edges_pos = logical(edges_pos);
                edges_neg = logical(edges_neg);

                % build model on TRAIN subs
                this.all_pos_edges(:,i_fold) = edges_pos;
                this.all_pos_edges = logical(this.all_pos_edges);
                this.all_neg_edges(:,i_fold) = edges_neg;
                this.all_neg_edges = logical(this.all_neg_edges);
                train_sum = (sum(train.x(this.all_pos_edges(:,i_fold), :), 1) - sum(train.x(this.all_neg_edges(:,i_fold), :), 1))';

                % adding the pos and neg edges here because they
                % weren't thresholded based on correlations
                % train_sum = (sum(train.x(this.all_pos_edges(:,i_fold), :), 1) + sum(train.x(this.all_neg_edges(:,i_fold), :), 1))';
                fit_train = polyfit(train_sum, train.y(:,1), 1);
                %fit_train = polyfit(train_sum, train.y, 1);

                % run model on TEST sub
                test_sum = sum(test.x(this.all_pos_edges(:,i_fold), :), 1) - sum(test.x(this.all_neg_edges(:,i_fold),:), 1);

                this.Y(test.indx) = (test_sum*fit_train(1)+fit_train(2))';
            end
        end
        function this = evaluate(this)
            [this.r_pearson, this.p_pearson] = corr(this.Y, this.phenotype.all_behav);
            [this.r_rank, this.p_rank] = corr(this.Y, this.phenotype.all_behav, 'type', 'spearman');
            this.mse = sum((this.Y - this.phenotype.all_behav).^2) / this.num_sub_total;
            this.q_s = 1 - this.mse / var(this.phenotype.all_behav, 1);
        end
    end
end
