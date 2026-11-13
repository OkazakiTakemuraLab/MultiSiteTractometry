function ICCheatmap_T1wB0

% This script aims to replicate a panel in Supplementary Figure 16 of the following article:
% Taguma, D., Yokoi, I., Kinjo, T., Tsuchida, S., Miyata, T., Matsuda, T., Lerma-Usabiaga, G. & Takemura, H. (2026) {title}. {journal info}
%
% Daiki Taguma, NIPS SCBM/SOKENDAI

%% Define base directories and subject information
baseDirPairs = {
    '../../TractometryData/Protocol1/'
    '../../TractometryData/Protocol2/'
};

subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05'};

%% Define sessions
sessions_p = {'ses-VerioA1','ses-VerioA2','ses-VerioB1','ses-VerioB2','ses-PrismaFitA1'};

% Heatmap labels
sessions = [ ...
    strcat('Protocol1-', sessions_p), ...
    strcat('Protocol2-', sessions_p) ...
];

prop = 't1w_over_b0.csv';
propCols = 1:20;

nFeat = numel(propCols);
nSub  = numel(subjects);
nSes  = numel(sessions);

%% Add external tools to the path
addpath(genpath('../../ExternalTools'))

%% 1) Read all files once: data(sub, ses, feat)
data = nan(nSub, nSes, nFeat);

for sesIdx = 1:nSes
    if startsWith(sessions{sesIdx}, 'Protocol1-')
        baseDir = baseDirPairs{1};
        ses     = extractAfter(sessions{sesIdx}, 'Protocol1-');
    else
        baseDir = baseDirPairs{2};
        ses     = extractAfter(sessions{sesIdx}, 'Protocol2-');
    end

    for subIdx = 1:nSub
        subjectID = subjects{subIdx};
        filePath  = fullfile(baseDir, subjectID, ses, prop);

        A = readmatrix(filePath);
        data(subIdx, sesIdx, :) = mean(A(1:100, propCols), 1);
    end
end

%% 2) Compute session-to-session ICC (Original)
R_orig = nan(nSes, nSes);

for i = 1:nSes
    R_orig(i,i) = 1;
    for j = i+1:nSes
        Xi = squeeze(data(:, i, :));  % nSub × nFeat
        Xj = squeeze(data(:, j, :));  % nSub × nFeat

        ns = size(Xi,1);

        all_r = nan(ns,1);
        type='A-1'; alpha=0.05; r0=0;
        for k = 1:ns
            T = [Xi(k,:); Xj(k,:)]';  % nFeat × 2
            [r, ~, ~, ~, ~, ~, ~] = ICC(T, type, alpha, r0);
            all_r(k) = r;
        end
        R_orig(i,j) = mean(all_r,'omitnan');
    end
end

R_orig = triu(R_orig) + triu(R_orig,1)';

%% 3) Compute session-to-session ICC (ComBat-harmonized)
% If the two sessions are under the same condition
% (= same protocol and same scanner), ComBat is not applied.
R_cb = nan(nSes, nSes);

for i = 1:nSes
    R_cb(i,i) = 1;
    for j = i+1:nSes
        Xi = squeeze(data(:, i, :));  % nSub × nFeat
        Xj = squeeze(data(:, j, :));  % nSub × nFeat

        % ----- Check whether the two sessions are the same condition -----
        cond_i = get_condition_label(sessions{i});
        cond_j = get_condition_label(sessions{j});

        all_r = nan(nSub,1);
        type='A-1'; alpha=0.05; r0=0;

        if strcmp(cond_i, cond_j)
            % Same condition -> do NOT apply ComBat
            for k = 1:nSub
                T = [Xi(k,:); Xj(k,:)]';  % nFeat × 2
                [r, ~, ~, ~, ~, ~, ~] = ICC(T, type, alpha, r0);
                all_r(k) = r;
            end
        else
            % Different condition -> apply ComBat
            select_data  = [Xi; Xj]';                        % nFeat × (2*nSub)
            batch_labels = [ones(nSub,1); 2*ones(nSub,1)];  % (2*nSub) × 1
            mod = [];

            data_harmonized = combat(select_data, batch_labels, mod, 1);
            Xh = data_harmonized';  % (2*nSub) × nFeat

            for k = 1:nSub
                T = [Xh(k,:); Xh(k+nSub,:)]';  % nFeat × 2
                [r, ~, ~, ~, ~, ~, ~] = ICC(T, type, alpha, r0);
                all_r(k) = r;
            end
        end

        R_cb(i,j) = mean(all_r,'omitnan');
    end
end

R_cb = triu(R_cb) + triu(R_cb,1)';

%% 4) Plot: left=Original, right=ComBat-harmonized
fig = figure('Position',[100,100,1530,600], 'Renderer','opengl');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Left: Original
nexttile;
h1 = heatmap(sessions, sessions, R_orig);
h1.Title = 'ICC Heatmap for T1w/b0 (Original)';
h1.XLabel = 'Session';
h1.YLabel = 'Session';
h1.Colormap = jet;
h1.ColorLimits = [0 1];
h1.CellLabelFormat = '%.2f';

% Right: ComBat
nexttile;
h2 = heatmap(sessions, sessions, R_cb);
h2.Title = 'ICC Heatmap for T1w/b0 (ComBat-harmonized)';
h2.XLabel = 'Session';
h2.YLabel = 'Session';
h2.Colormap = jet;
h2.ColorLimits = [0 1];
h2.CellLabelFormat = '%.2f';

%% output figure as eps
exportgraphics(fig, 'ICCheatmap_T1wB0.eps');

end

%% ===== local function =====
function cond_label = get_condition_label(session_name)
% Convert session label into a "condition" label defined by
% protocol + scanner identity, while ignoring repeated scan number.
%
% Examples:
%   Protocol1-ses-VerioA1 -> Protocol1-VerioA
%   Protocol1-ses-VerioA2 -> Protocol1-VerioA
%   Protocol2-ses-VerioB1 -> Protocol2-VerioB

parts = split(session_name, '-');
protocol = parts{1};      % Protocol1 or Protocol2
scanner  = parts{3};      % VerioA1, VerioA2, VerioB1, ...

% Remove the trailing scan number
scanner_base = regexprep(scanner, '\d+$', '');

cond_label = [protocol '-' scanner_base];
end