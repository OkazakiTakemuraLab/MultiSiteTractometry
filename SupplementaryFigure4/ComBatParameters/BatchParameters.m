function BatchParameters

% This script aims to replicate a panel in Supplementary Figure 4 of the following article:
% Taguma, D., Yokoi, I., Kinjo, T., Tsuchida, S., Miyata, T., Matsuda, T., Lerma-Usabiaga, G. & Takemura, H. (2026) {title}. {journal info}
%
% Daiki Taguma, NIPS SCBM/SOKENDAI

%% Settings

baseDirs = {'../../TractometryData/Protocol1', ...
            '../../TractometryData/Protocol2'};

subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05'};

sessions = {'ses-VerioA1','ses-VerioA2', ...
            'ses-VerioB1','ses-VerioB2', ...
            'ses-PrismaFitA1','ses-PrismaFitA2'};

props = {'md.csv','fa.csv','icvf.csv', ...
         'odi.csv','isov.csv','t1w_over_b0.csv'};

propNames = {'MD','FA','ICVF','ODI','IsoV','T1w/b0'};

addpath(genpath('../../ExternalTools'))

nSub = numel(subjects);
nSes = numel(sessions);
nProt = numel(baseDirs);

% Protocol and session indices
[sesIdx,protIdx] = ndgrid(1:nSes,1:nProt);
sesIdx = sesIdx(:);
protIdx = protIdx(:);

nCond = numel(sesIdx);

% Scanner identity (remove repetition number)
scanners = regexprep(sessions,'\d+$','');

parameters = cell(2,2,numel(props)); % gamma/delta x hat/star


%% ComBat

for p = 1:numel(props)

    data = nan(nSub,nCond,20);

    % Read data for each condition
    for c = 1:nCond

        % T1w/b0 excludes PrismaFitA2
        if p == 6 && sesIdx(c) == 6
            continue
        end

        for s = 1:nSub

            filename = fullfile( ...
                baseDirs{protIdx(c)}, ...
                subjects{s}, ...
                sessions{sesIdx(c)}, ...
                props{p});

            A = readmatrix(filename);

            data(s,c,:) = mean(A(1:100,1:20),1);

        end
    end

    % Compare all pairs of different conditions
    for i = 1:nCond

        % T1w/b0 excludes PrismaFitA2
        if p == 6 && sesIdx(i) == 6
            continue
        end

        for j = i+1:nCond

            % T1w/b0 excludes PrismaFitA2
            if p == 6 && sesIdx(j) == 6
                continue
            end

            % Skip repeated scans of the same condition
            if protIdx(i) == protIdx(j) && ...
                    strcmp(scanners{sesIdx(i)},scanners{sesIdx(j)})
                continue
            end

            X = [squeeze(data(:,i,:)); ...
                 squeeze(data(:,j,:))]';

            batch = [ones(nSub,1); 2*ones(nSub,1)];

            [~,est] = combat(X,batch,[],1);

            % Convert gamma to original metric scale
            mu = [mean(X(:,1:nSub),2), ...
                  mean(X(:,nSub+1:end),2)];

            residuals = [X(:,1:nSub)-mu(:,1), ...
                         X(:,nSub+1:end)-mu(:,2)];

            scale = sqrt(mean(residuals.^2,2))';

            gamma = {est.gamma_hat .* scale, ...
                     est.gamma_star .* scale};

            delta = {est.delta_hat,est.delta_star};

            for h = 1:2
                parameters{1,h,p} = [parameters{1,h,p}; gamma{h}(:)];
                parameters{2,h,p} = [parameters{2,h,p}; delta{h}(:)];
            end

        end
    end
end


%% Figure

fig = figure('Color','w','Units','centimeters', ...
             'Position',[2 2 36 16]);

t = tiledlayout(2,7,'TileSpacing','compact','Padding','compact');

colors = colororder;
hatColor = colors(1,:);
starColor = colors(2,:);

iqrColor = [0.4 0.4 0.4];
medianColor = [0.5 0.5 0.5];

for p = 1:numel(props)

    for m = 1:2

        ax = nexttile(t,p+(m-1)*7);
        hold(ax,'on')

        hat = parameters{m,1,p};
        star = parameters{m,2,p};

        violinplot({star,hat},[], ...
            'ViolinColor',{starColor,hatColor}, ...
            'ViolinAlpha',{0.3,0.3}, ...
            'ShowData',false, ...
            'ShowNotches',false, ...
            'ShowMean',false, ...
            'ShowBox',false, ...
            'ShowMedian',false, ...
            'ShowWhiskers',false, ...
            'QuartileStyle','none');

        % IQR and median
        q = [prctile(hat,[25 50 75]); ...
             prctile(star,[25 50 75])];

        for h = 1:2
            x = 0.92 + (h-1)*0.16;

            plot(ax,[x x],q(h,[1 3]), ...
                'Color',iqrColor,'LineWidth',2);

            scatter(ax,x,q(h,2),25,'w','filled', ...
                'MarkerEdgeColor',medianColor, ...
                'LineWidth',1);
        end

        if m == 1
            title(ax,propNames{p},'FontSize',18);
        else
            ylim(ax,[0 3])
        end

        xticks(ax,[])

        if p == 1
            labels = {'\gamma','\delta'};
            yl = ylabel(ax,labels{m},'FontSize',15);
            yl.Units = 'normalized';
            yl.Position(1) = -0.30;

            if m == 1
                text(ax,-0.26,1.06,'(\mum^2/ms)', ...
                    'Units','normalized', ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom', ...
                    'FontSize',14);
            end
        end

        set(ax,'FontName','Arial','FontSize',18, ...
            'LineWidth',1,'TickDir','out','Box','off');

    end
end


%% Legend

ax = nexttile(t,7,[2 1]);
hold(ax,'on')
axis(ax,'off')

h(1) = patch(ax,nan,nan,hatColor, ...
    'FaceAlpha',0.3,'EdgeColor',[0.5 0.5 0.5]);

h(2) = patch(ax,nan,nan,starColor, ...
    'FaceAlpha',0.3,'EdgeColor',[0.5 0.5 0.5]);

h(3) = plot(ax,nan,nan,'Color',iqrColor,'LineWidth',2);

h(4) = scatter(ax,nan,nan,25,'w','filled', ...
    'MarkerEdgeColor',medianColor,'LineWidth',1);

legend(ax,h,{'hat','star','IQR','Median'}, ...
    'Location','northwest','FontName','Arial', ...
    'FontSize',18,'Box','off');


%% Export

exportgraphics(fig,'ComBat_batch_parameters.eps');

end