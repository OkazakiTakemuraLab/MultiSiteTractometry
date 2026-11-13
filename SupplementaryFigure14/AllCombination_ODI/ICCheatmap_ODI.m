function ICCheatmap_ODI

% This script aims to replicate a panel in Supplementary Figure 14 of the following article:
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
sessions_p = {'ses-VerioA1','ses-VerioA2','ses-VerioB1','ses-VerioB2','ses-PrismaFitA1','ses-PrismaFitA2'};

% Heatmap labels
sessions = [ ...
    strcat('Protocol1-', sessions_p), ...
    strcat('Protocol2-', sessions_p) ...
];

prop = 'odi.csv';
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

        ns = size(Xi,1); % = nSub

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
% Same condition (same protocol AND same scanner) -> DO NOT apply ComBat.

R_cb = nan(nSes, nSes);

% helper: get protocol ("Protocol1"/"Protocol2") and scanner name ("VerioA", etc.)
getProt = @(s) regexp(s, '^(Protocol\d+)-', 'tokens', 'once');
getSes  = @(s) regexp(s, '^(?:Protocol\d+)-(.+)$', 'tokens', 'once'); % raw session string after ProtocolX-
getScan = @(ses) regexprep(ses, '\d+$', ''); % remove trailing run number(s), e.g., VerioA1 -> VerioA

for i = 1:nSes
    R_cb(i,i) = 1;

    prot_i = getProt(sessions{i});  prot_i = prot_i{1};
    ses_i  = getSes(sessions{i});   ses_i  = ses_i{1};
    scan_i = getScan(ses_i);

    for j = i+1:nSes
        prot_j = getProt(sessions{j});  prot_j = prot_j{1};
        ses_j  = getSes(sessions{j});   ses_j  = ses_j{1};
        scan_j = getScan(ses_j);

        Xi = squeeze(data(:, i, :));  % nSub × nFeat
        Xj = squeeze(data(:, j, :));  % nSub × nFeat

        sameCondition = strcmp(prot_i, prot_j) && strcmp(scan_i, scan_j);

        type='A-1'; alpha=0.05; r0=0;

        if sameCondition
            % ---- NO ComBat: compute ICC directly from raw data ----
            all_r = nan(nSub,1);
            for k = 1:nSub
                T = [Xi(k,:); Xj(k,:)]';  % nFeat × 2
                [r, ~, ~, ~, ~, ~, ~] = ICC(T, type, alpha, r0);
                all_r(k) = r;
            end
            R_cb(i,j) = mean(all_r,'omitnan');

        else
            % ---- Apply ComBat only when protocol or scanner differs ----
            select_data   = [Xi; Xj]';                      % nFeat × (2*nSub)
            batch_labels  = [ones(nSub,1); 2*ones(nSub,1)];  % (2*nSub)×1
            mod = [];
            data_harmonized = combat(select_data, batch_labels, mod, 1);
            Xh = data_harmonized';  % (2*nSub) × nFeat

            all_r = nan(nSub,1);
            for k = 1:nSub
                T = [Xh(k,:) ; Xh(k+nSub,:)];               % 2 × nFeat
                [r, ~, ~, ~, ~, ~, ~] = ICC(T', type, alpha, r0); % nFeat × 2
                all_r(k) = r;
            end
            R_cb(i,j) = mean(all_r,'omitnan');
        end
    end
end

R_cb = triu(R_cb) + triu(R_cb,1)';

%% 4) Plot: left=Original, right=ComBat-harmonized
fig = figure('Position',[100,100,1530,600], 'Renderer','opengl');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Left: Original
nexttile;
h1 = heatmap(sessions, sessions, R_orig);
h1.Title = 'ICC Heatmap for ODI (Original)';
h1.XLabel = 'Session';
h1.YLabel = 'Session';
h1.Colormap = jet;
h1.ColorLimits = [0 1];
h1.CellLabelFormat = '%.2f';

% Right: ComBat
nexttile;
h2 = heatmap(sessions, sessions, R_cb);
h2.Title = 'ICC Heatmap for ODI (ComBat-harmonized)';
h2.XLabel = 'Session';
h2.YLabel = 'Session';
h2.Colormap = jet;
h2.ColorLimits = [0 1];
h2.CellLabelFormat = '%.2f';

%% output figure as eps
exportgraphics(fig, 'ICCheatmap_ODI.eps');

end
