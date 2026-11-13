function QC_violin

% This script aims to replicate a panel in Supplementary Figure 1 of the following article:
% Taguma, D., Yokoi, I., Kinjo, T., Tsuchida, S., Miyata, T., Matsuda, T., Lerma-Usabiaga, G. & Takemura, H. (2026) {title}. {journal info}
%
% Daiki Taguma, NIPS SCBM/SOKENDAI

baseDir = '/fileserver/project/proj_gari/code4paper/TractometryData';

protocols = {'Protocol1','Protocol2'};
subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05'};

% Subject-specific color and marker
colors = colororder;
colors = colors(1:numel(subjects),:);

markers = {'o','s','^','d','p'};

%% Add path to external tools
addpath(genpath('../../ExternalTools'))

violinColor = [0.82 0.82 0.82];
violinEdgeColor = [0.30 0.30 0.30];
violinWidth = 0.30;

jitterWidth = 0.13;
markerSize = 32;

rng(1);

Protocol = {};
Subject = {};

AverageAbsMotion = [];
AverageRelMotion = [];
OutlierSlicesPercent = [];

%% ============================================================
% Load QC data
% =============================================================

for p = 1:numel(protocols)

    for s = 1:numel(subjects)

        subjectDir = fullfile( ...
            baseDir, ...
            protocols{p}, ...
            subjects{s});

        sessions = dir(fullfile(subjectDir, 'ses-*'));
        sessions = sessions([sessions.isdir]);

        for j = 1:numel(sessions)

            session = sessions(j).name;
            sessionDir = fullfile(subjectDir, session);

            movementFile = fullfile( ...
                sessionDir, ...
                'eddy_movement_rms');

            outlierFile = fullfile( ...
                sessionDir, ...
                'eddy_outlier_map');

            if ~isfile(movementFile) || ~isfile(outlierFile)
                continue
            end

            M = readmatrix( ...
                movementFile, ...
                'FileType','text');

            O = readmatrix( ...
                outlierFile, ...
                'FileType','text', ...
                'NumHeaderLines',1);

            M = M(~all(isnan(M),2),:);
            M = M(:,~all(isnan(M),1));

            O = O(~all(isnan(O),2),:);
            O = O(:,~all(isnan(O),1));

            Protocol{end+1,1} = protocols{p};
            Subject{end+1,1} = subjects{s};

            AverageAbsMotion(end+1,1) = ...
                mean(M(:,1),'omitnan');

            AverageRelMotion(end+1,1) = ...
                mean(M(:,2),'omitnan');

            OutlierSlicesPercent(end+1,1) = ...
                100 * sum(O(:)>0) / numel(O);

        end
    end
end


%% ============================================================
% Create table
% =============================================================

T = table( ...
    Protocol, ...
    Subject, ...
    AverageAbsMotion, ...
    AverageRelMotion, ...
    OutlierSlicesPercent);


%% ============================================================
% Metrics
% =============================================================

metrics = {
    T.AverageAbsMotion
    T.AverageRelMotion
    T.OutlierSlicesPercent
};

yLabels = {
    'Average abs. motion (mm)'
    'Average rel. motion (mm)'
    'Outlier slices (%)'
};


%% ============================================================
% Figure
% =============================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[3 3 29 9]);

tl = tiledlayout( ...
    fig, ...
    1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

axesHandles = gobjects(3,1);


%% ============================================================
% Plot
% =============================================================

for k = 1:3

    ax = nexttile(tl);
    axesHandles(k) = ax;

    hold(ax,'on');

    values = metrics{k};


    %% Violin plot

    v = violinplot( ...
        values, ...
        [], ...
        'Width',violinWidth, ...
        'ViolinColor',violinColor, ...
        'ViolinAlpha',1, ...
        'EdgeColor',violinEdgeColor, ...
        'ShowData',false, ...
        'ShowNotches',false, ...
        'ShowMean',false, ...
        'ShowBox',false, ...
        'ShowMedian',false, ...
        'ShowWhiskers',false, ...
        'QuartileStyle','shadow');

    if ~isempty(v(1).ViolinPlotQ)
        v(1).ViolinPlotQ.Visible = 'off';
    end


    %% Median + IQR

    med = median(values,'omitnan');
    q = prctile(values,[25 75]);

    plot( ...
        ax, ...
        [1 1], ...
        q, ...
        'k-', ...
        'LineWidth',2);

    plot( ...
        ax, ...
        [0.88 1.12], ...
        [med med], ...
        'k-', ...
        'LineWidth',2);


    %% Individual datasets
    % Color + shape = subject
    % All markers are open

    for i = 1:height(T)

        s = find( ...
            strcmp(subjects,T.Subject{i}), ...
            1);

        x = 1 + ...
            (rand - 0.5) * 2 * jitterWidth;

        scatter( ...
            ax, ...
            x, ...
            values(i), ...
            markerSize, ...
            'Marker',markers{s}, ...
            'MarkerEdgeColor',colors(s,:), ...
            'MarkerFaceColor','none', ...
            'LineWidth',1.2);

    end


    %% Axis formatting

    xlim(ax,[0.52 1.48]);

    xticks(ax,1);
    xticklabels(ax,{''});

    ylabel(ax,yLabels{k});

    set( ...
        ax, ...
        'FontName','Arial', ...
        'FontSize',15, ...
        'LineWidth',1, ...
        'TickDir','out', ...
        'Box','off');

    % Show y-axis tick labels with one decimal place
    % e.g. 1.0, 2.0, 3.0
    ytickformat(ax,'%.1f');

end


%% ============================================================
% Legend: subject only
% =============================================================

ax = axesHandles(3);

h = gobjects(numel(subjects),1);

for s = 1:numel(subjects)

    h(s) = scatter( ...
        ax, ...
        nan,nan, ...
        markerSize, ...
        'Marker',markers{s}, ...
        'MarkerEdgeColor',colors(s,:), ...
        'MarkerFaceColor','none', ...
        'LineWidth',1.2);

end

legend( ...
    ax, ...
    h, ...
    subjects, ...
    'Location','eastoutside', ...
    'Box','off', ...
    'FontName','Arial', ...
    'FontSize',14);


%% ============================================================
% Export
% =============================================================

exportgraphics(fig, 'QC_violin.eps');

end