function TractProfile_MD_sub02

% This script aims to replicate a panel in Figure 3 of the following article:
% Taguma, D., Yokoi, I., Kinjo, T., Tsuchida, S., Miyata, T., Matsuda, T., Lerma-Usabiaga, G. & Takemura, H. (2026) {title}. {journal info}
%
% Daiki Taguma, NIPS SCBM/SOKENDAI

%% Define the base directory and subject information
baseDir1 = '../../TractometryData/Protocol1/';
baseDir2 = '../../TractometryData/Protocol2/';
subjects = {'sub-02'};

%% Define session information and property file
sessions = {'ses-VerioA1', 'ses-VerioA2', 'ses-VerioB1', 'ses-VerioB2', 'ses-PrismaFitA1', 'ses-PrismaFitA2', 'ses-PrismaFitB1', 'ses-PrismaFitB2'};
prop = 'md.csv';
md(100, 2*numel(sessions)) = 0;
for subIdx = 1:numel(subjects)
    subjectID = subjects{subIdx};
    for sesIdx = 1:numel(sessions)
        sessionID = sessions{sesIdx};
        % Load prop data
        filePath = fullfile(baseDir1, subjectID, sessionID, prop);
        prop_data = readmatrix(filePath);
        md(1:100, sesIdx) = prop_data(1:100, 3);
    end
end
for subIdx = 1:numel(subjects)
    subjectID = subjects{subIdx};
    for sesIdx = 1:numel(sessions)
        sessionID = sessions{sesIdx};
        % Load prop data
        filePath = fullfile(baseDir2, subjectID, sessionID, prop);
        prop_data = readmatrix(filePath);
        md(1:100, sesIdx + numel(sessions)) = prop_data(1:100, 3);
    end
end

%% Plot the tract profile
fig = figure;
hold on
plot([1:1:100],md(1:100,1),'color','#CE7E8B','LineWidth',2);
plot([1:1:100],md(1:100,2),'color','#CE7E8B','LineWidth',2, 'LineStyle', '--');
plot([1:1:100],md(1:100,3),'color','#962676','LineWidth',2);
plot([1:1:100],md(1:100,4),'color','#962676','LineWidth',2, 'LineStyle', '--');
plot([1:1:100],md(1:100,5),'color','#FFD600','LineWidth',2);
plot([1:1:100],md(1:100,6),'color','#FFD600','LineWidth',2,'LineStyle', '--');
plot([1:1:100],md(1:100,7),'color','#F8A339','LineWidth',2);
plot([1:1:100],md(1:100,8),'color','#F8A339','LineWidth',2, 'LineStyle', '--');
plot([1:1:100],md(1:100,9),'color','#CE7E8B','LineWidth',2,'LineStyle', '-.');
plot([1:1:100],md(1:100,10),'color','#CE7E8B','LineWidth',2, 'LineStyle', ':');
plot([1:1:100],md(1:100,11),'color','#962676','LineWidth',2,'LineStyle', '-.');
plot([1:1:100],md(1:100,12),'color','#962676','LineWidth',2, 'LineStyle', ':');
plot([1:1:100],md(1:100,13),'color','#FFD600','LineWidth',2,'LineStyle', '-.');
plot([1:1:100],md(1:100,14),'color','#FFD600','LineWidth',2, 'LineStyle', ':');
plot([1:1:100],md(1:100,15),'color','#F8A339','LineWidth',2,'LineStyle', '-.');
plot([1:1:100],md(1:100,16),'color','#F8A339','LineWidth',2, 'LineStyle', ':');
mean_Protocol1_VerioA1 = mean(md(1:100,1));
mean_Protocol1_VerioA2 = mean(md(1:100,2));
mean_Protocol1_VerioB1 = mean(md(1:100,3));
mean_Protocol1_VerioB2 = mean(md(1:100,4));
mean_Protocol1_PrismaFitA1 = mean(md(1:100,5));
mean_Protocol1_PrismaFitA2 = mean(md(1:100,6));
mean_Protocol1_PrismaFitB1 = mean(md(1:100,7));
mean_Protocol1_PrismaFitB2 = mean(md(1:100,8));
mean_Protocol2_VerioA1 = mean(md(1:100,9));
mean_Protocol2_VerioA2 = mean(md(1:100,10));
mean_Protocol2_VerioB1 = mean(md(1:100,11));
mean_Protocol2_VerioB2 = mean(md(1:100,12));
mean_Protocol2_PrismaFitA1 = mean(md(1:100,13));
mean_Protocol2_PrismaFitA2 = mean(md(1:100,14));
mean_Protocol2_PrismaFitB1 = mean(md(1:100,15));
mean_Protocol2_PrismaFitB2 = mean(md(1:100,16));
mean_Protocol1_VerioA1
mean_Protocol1_VerioA2
mean_Protocol1_VerioB1
mean_Protocol1_VerioB2
mean_Protocol1_PrismaFitA1
mean_Protocol1_PrismaFitA2
mean_Protocol1_PrismaFitB1
mean_Protocol1_PrismaFitB2
mean_Protocol2_VerioA1
mean_Protocol2_VerioA2
mean_Protocol2_VerioB1
mean_Protocol2_VerioB2
mean_Protocol2_PrismaFitA1
mean_Protocol2_PrismaFitA2
mean_Protocol2_PrismaFitB1
mean_Protocol2_PrismaFitB2

%% Set up plot appearance
set(gca,'tickdir','out', ...
    'box','off', ...
    'ylim',[0.6 1.2],'ytick',[0.6 0.9 1.2])
ylabel('MD','fontsize',16);
xlabel('Position','fontsize',16);

%% output figure as eps
exportgraphics(fig, 'tractprofile_MD.eps');
end