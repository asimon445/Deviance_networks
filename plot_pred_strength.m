% This script will plot prediction strength for all measures from all types
% of predictions (deviance/whole brain & kernel ridge/linear) in descending
% order

clear;

path = '/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/Full_prediction_results_new/';

net = {'highSymp_devnet','lowSymp_devnet'};
type = {'Wholebrain_preds_krCPM','Wholebrain_preds_CPM','Deviance_preds_krCPM','Deviance_preds_CPM'};

for n=1:length(net)

    files = dir([path sprintf('*%s.mat',net{n})]);
    
    ishigh = strfind(net{n},'highSymp');
    
    for f = 1:length(files)
        load([files(f).folder '/' files(f).name]);

        if ~isempty(ishigh)
            measure{f,1} = files(f).name(1:end-20);
        else
            measure{f,1} = files(f).name(1:end-19);
        end

        for t = 1:length(type)
            isdev = strfind(type{t},'Deviance');
            if ~isempty(isdev)

                for e = 1:length(Deviance_network.target_n_edges)

                    eval(sprintf('%s_%d_edges(:,f) = Deviance_network.%s(:,e);',type{t},Deviance_network.target_n_edges(e),type{t}));
                    eval(sprintf('%s_%d_edges_median(1,f) = median(%s_%d_edges(:,f));',...
                        type{t},Deviance_network.target_n_edges(e),type{t},Deviance_network.target_n_edges(e)));
                end
            else
                eval(sprintf('%s(:,f) = Deviance_network.%s(:,1);',type{t},type{t}));
                    eval(sprintf('%s_median(1,f) = median(%s(:,f));',...
                        type{t},type{t}));
            end
        end
    end

    for t = 1:length(type)
        isdev = strfind(type{t},'Deviance');
        if ~isempty(isdev)
            for e = 1:length(Deviance_network.target_n_edges)
                % sort in descending order
                [sorted,ix] = eval(sprintf('sort(%s_%d_edges_median(1,:),''ascend'')',type{t},Deviance_network.target_n_edges(e)));

                for i = 1:length(ix)
                    data(:,i) = eval(sprintf('%s_%d_edges(:,ix(i))',type{t},Deviance_network.target_n_edges(e)));
                    header{1,i} = measure{ix(i),1};
                end

                % plot
                figure;
                boxplot(data,'Orientation','horizontal','symbol', '');
                h=findobj('LineStyle','--');
                set(h, 'LineStyle','-');
                xlim([-0.2 0.45]);
                set(gca,'YTickLabel',header,'TickLabelInterpreter','none');
                eval(sprintf('title(''%s_%d_edges'',''interpreter'',''none'')',type{t},Deviance_network.target_n_edges(e)));
                xlabel('Prediction strength (rho)');

                clear data header h sorted ix
            end
        else   % whole brain

            [sorted,ix] = eval(sprintf('sort(%s_median)',type{t}));

            for i = 1:length(ix)
                data(:,i) = eval(sprintf('%s(:,ix(i))',type{t}));
                header{1,i} = measure{ix(i),1};
            end

            % plot
            figure;
            boxplot(data,'Orientation','horizontal','symbol', '');
            h=findobj('LineStyle','--');
            set(h, 'LineStyle','-');
            xlim([-0.2 0.4]);
            set(gca,'YTickLabel',header,'TickLabelInterpreter','none');
            eval(sprintf('title(''%s'',''interpreter'',''none'')',type{t}));
            xlabel('Prediction strength (rho)');

            clear data header h sorted ix
        end
    end
    clear files measure
end


