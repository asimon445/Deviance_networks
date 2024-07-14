classdef kernel_rcpm_forceEdges < predictory
    methods
        function this = kernel_rcpm_forceEdges(group,options)
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

                   % break it
                    error('This function is not setup to handle multiple controls'); 
                else
                     % put 'kernel_prediction' stuff here
                     lambda = [0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];
                     scale = 0.2;
                     type = {'corr','Exponential','Gaussian'};

                     for l = 1:length(lambda)
                         for t = 1:length(type)

                            [pred, estimate] = kernel_prediction(train.x, test.x, train.y, lambda(l), type{t}, scale);
                    
                            eval(sprintf('this.Y_%s(test.indx,l) = pred;',type{t}));
                            clear pred estimate
                         end
                     end
                end
                
            end
        end
        function this = evaluate(this)
            [this.r_pearson, this.p_pearson] = corr(this.Y, this.phenotype.all_behav);
            type = {'corr','Exponential','Gaussian'};
            lambda = [0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];

            for l = 1:length(lambda)
                for t = 1:length(type)
                    eval(sprintf('[this.r_rank_%s, this.p_rank_%s] = corr(this.Y_%s, this.phenotype.all_behav, ''type'', ''spearman'');',type{t},type{t},type{t}));
                end
            end

            this.mse = sum((this.Y - this.phenotype.all_behav).^2) / this.num_sub_total;
            this.q_s = 1 - this.mse / var(this.phenotype.all_behav, 1);
%             fprintf('q_s=%f\n',this.q_s);
%             fprintf('spearman=%f\n',this.r_rank);
%             fprintf('spearman pval=%.30f\n',this.p_rank);

        end
    end
end
