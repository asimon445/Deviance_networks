% This script will plot prediction strength for all measures from all types
% of predictions (deviance/whole brain & kernel ridge/linear) in descending
% order

clear;

path = '/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/Deviance_prediction_results/LowScoresModeled/Cog/';

outpath = '/Users/ajsimon/Documents/Data/Constable_lab/Deviance_networks/Deviance_transdiagnostic/Deviance_prediction_results/';
outfile = 'Symptom_highend_deviance_predictions.mat';

type = {'Deviance_preds_krCPM','Wholebrain_preds_krCPM'}; %{'Wholebrain_preds_krCPM','Wholebrain_preds_CPM','Deviance_preds_krCPM','Deviance_preds_CPM'};

files = dir([path '*.mat']);

ishigh = strfind(path,'HighScores');

ref = 0;   % use 'ref = 1' to sort the prediction strength based on deviance networks, use 'ref = 0' to sort based on whole brain predictions
thresh = 0.15;  % only plot coefficients greater than this

delta = linspace(0.4,0.7,2);
width = 0.2;
colors = {[0.0 0.65 0.65] [0.75 0 0]};  

%% Pull the prediction strength data
idx = 0;
for f = 1:length(files)

    israw = strfind(files(f).name,'_raw');

    if isempty(israw)
        idx=idx+1;

        load([files(f).folder '/' files(f).name]);

        if ~isempty(ishigh)
            measure{idx,1} = files(f).name(1:end-20);
        else
            measure{idx,1} = files(f).name(1:end-19);
        end

        for t = 1:length(type)
            isdev = strfind(type{t},'Deviance');
            if ~isempty(isdev)

                for e = 1:length(Deviance_network.target_n_edges)

                    eval(sprintf('%s_%d_edges(:,idx) = Deviance_network.%s(:,e);',type{t},Deviance_network.target_n_edges(e),type{t}));
                    eval(sprintf('%s_%d_edges_median(1,idx) = median(%s_%d_edges(:,idx));',...
                        type{t},Deviance_network.target_n_edges(e),type{t},Deviance_network.target_n_edges(e)));
                end
            else
                eval(sprintf('%s(:,idx) = Deviance_network.%s(:,1);',type{t},type{t}));
                eval(sprintf('%s_median(1,idx) = median(%s(:,idx));',...
                    type{t},type{t}));
            end
        end
    end
end

%% Plot results
for e = 1:length(Deviance_network.target_n_edges)
    fig = figure;
    for t = 1:length(type)
        isdev = strfind(type{t},'Deviance');

        if ref == 1
            [sorted,ix] = eval(sprintf('sort(Deviance_preds_krCPM_%d_edges_median(1,:),''ascend'')',Deviance_network.target_n_edges(e)));
            toplot = find(sorted>thresh);
            ix_use = ix(toplot);
        elseif ref == 0
            [sorted,ix] = sort(Wholebrain_preds_krCPM_median);
            toplot = find(sorted>thresh);
            ix_use = ix(toplot);
        else
            error('You can only set ref = 0 or 1. Use a different ref code')
        end

        for i = 1:length(ix_use)
            header{1,i} = measure{ix_use(i),1};
        end

        if ~isempty(isdev)
            for i = 1:length(ix_use)
                data_dev(:,i) = eval(sprintf('%s_%d_edges(:,ix_use(i))',type{t},Deviance_network.target_n_edges(e)));
            end
            data = data_dev;
        else   % whole brain
            for i = 1:length(ix_use)
                data_wb(:,i) = eval(sprintf('%s(:,ix_use(i))',type{t}));
            end
            data = data_wb;
        end

        % plot
        boxplot(data,'Orientation','horizontal','symbol', '','Color',colors{t},'position',...
            (1:size(data,2))+delta(t),'widths',width);
        h=findobj('LineStyle','--');
        set(h, 'LineStyle','-');
        %xlim([-0.08 0.38]);
        xlim([-0.08 0.58]);
        set(gca,'YTickLabel',header,'TickLabelInterpreter','none');
        
        hold on;

        bh{t} = plot(nan,'Color',colors{t});

        h = findobj(gca, 'Tag', 'Box');
        for j = 1:size(data,2)
            patch(get(h(j), 'XData'), get(h(j), 'YData'), colors{t}, 'FaceAlpha', .9, 'EdgeColor', 'k', 'LineWidth', 0.5); % Set the face color, transparency, and edge color
        end
        clear h

        h=findobj('LineStyle','--');
        set(h, 'LineStyle','-','LineWidth',1,'k');
        clear h

        clear data
    end

    eval(sprintf('title(''%d edge deviance network predictions'',''interpreter'',''none'')',Deviance_network.target_n_edges(e)));
    xlabel('Prediction strength (rho)');

    h=gca;
    h.XAxis.TickLength = [0 0];
    h.YAxis.TickLength = [0 0];

    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(lines, 'Color', 'k','linewidth',3);

    lgnd{1,1} = 'Deviance network predictions';
    lgnd{1,2} = 'Whole brain predictions';

    [~, hobj, ~, ~] = legend([bh{:}],lgnd,'location','southeastoutside','FontSize',12);
    hl = findobj(hobj,'type','line');
    set(hl,'LineWidth',6);

    eval(sprintf('header_%d = header;',Deviance_network.target_n_edges(e)));
    clear sorted ix header data h lgnd h1 bh data_dev data_wb
end

clearvars -except Deviance_preds* header*
save([outpath outfile])

