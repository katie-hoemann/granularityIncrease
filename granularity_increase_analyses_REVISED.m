clear;
clc;

%% set parameters
removeNegativeICCs = 2; % set to 1 to set negative ICC values to 0; set to 2 to drop (NaN) those values
imputeData = 0; % set to 1 to impute values for missing days
dropFirstICC = 0; % set to 1 to drop first (non-NaN) ICC value
dropNegPPs = 0; % set to 1 drop participants with more than a few (>2) negative ICCs

%% load data file, along with word file that includes raw norms
dataFile = '18ARIEOD_daily_filtered.xlsx'; % does not include excluded subjects, days with < 6 completed prompts and late surveys
wordFile = 'words18.csv'; 
rawData = importdata(dataFile);
allData = rawData.data.NoLateSurveys;
subjectIDlist = unique(allData(:,1)); % grab subject IDs from first column of data file
words = readtable(wordFile); 
wordList = rawData.colheaders.NoLateSurveys(5:end)';  % grab sampled words from top row of data file

%% set valence categories for sampled words
for i_word = 1:height(words) % define valence categories
if words.Valence(i_word) > 5 % derived based on the database mean for Warriner et al (2013)
    valCat(i_word) = {'Positive'};
    positive(i_word) = 1;
    valence(i_word) = 1;
else
    valCat(i_word) = {'Negative'};
    positive(i_word) = 0;
    valence(i_word) = 2;
end
end 
words = [words valCat']; % append table with category assignments
words.Properties.VariableNames(5:end) = {'ValCat'}; % label new variables
labels = [positive']; % create matrix for logical indexing in ICC commands

%% grab data for each subject and run through calculations
for i_subject = 1:length(subjectIDlist)
    subjectData = [];
    subjectID = subjectIDlist(i_subject);
    index1 = find(allData(:,1)==subjectID);
    subjectData = allData(index1,2:end);
    dayIDlist = unique(subjectData(:,1));
    % remove invalid values and missing data
    missingData = isnan(subjectData); 
    missingData2 = any(missingData,2); 
    subjectData = subjectData(~missingData2,:);
    % calculate experience sampling statistics
    numPrompts(i_subject) = size(subjectData,1);
    numDays(i_subject) = length(dayIDlist);
    maxDays(i_subject) = max(dayIDlist);
    missingDays(i_subject) = maxDays(i_subject)-numDays(i_subject);
    numPrompts_day(i_subject) = numPrompts(i_subject)/numDays(i_subject);
    % calculate measures of experienced affect
    mPositive(i_subject) = mean(mean(subjectData(:,(valence==1))));
    mNegative(i_subject) = mean(mean(subjectData(:,(valence==2))));
    sdPositive(i_subject) = mean(std(subjectData(:,(valence==1))));
    sdNegative(i_subject) = mean(std(subjectData(:,(valence==2)))); 
    for i_day = 1:length(dayIDlist)
        dayData = [];
        dayID = dayIDlist(i_day);
        index2 = find(subjectData(:,1)==dayID);
        dayData = subjectData(index2,4:end);
        %% remove invalid values and missing data
        missingData = isnan(dayData); 
        missingData2 = any(missingData,2); 
        dayData = dayData(~missingData2,:); 
        %% compute ICCs
        rawICC(i_day,1) = ICC(dayData(:,(labels(:,1)==0)),'A-k'); % negative valence ICC
        rawICC(i_day,2) = ICC(dayData(:,(labels(:,1)==1)),'A-k'); % positive valence ICC
        %% count negative ICC value(s)
        if rawICC(i_day,:) < 0
            negICCs(i_day,:) = 1;
        else
            negICCs(i_day,:) = 0;
        end
        %% set lower bound of ICCs to 0 (no negative values)
        if removeNegativeICCs == 1
            rawICC(rawICC<0) = 0;
        elseif removeNegativeICCs == 2
            rawICC(rawICC<0) = NaN;
        end
        %% Fisher transform ICCs (noted as zICCs)
        rtoZ(i_day,1) = 0.5*log((1+rawICC(i_day,1))/(1-rawICC(i_day,1)));
        rtoZ(i_day,2) = 0.5*log((1+rawICC(i_day,2))/(1-rawICC(i_day,2)));
        %% invert zICCs for intuitive directionality
        zInv(i_day,1) = rtoZ(i_day,1)*-1;
        zInv(i_day,2) = rtoZ(i_day,2)*-1;
    end
    %% impute missing data, if desired
    if imputeData == 1
        zInv = fillmissing(zInv,'movmean',3,1);
    end
    %% add subject results to summary table
    ppID = ['PP' num2str(subjectID)];
    if i_subject == 1
        day = (1:1:max(allData(:,2)))';
        day = array2table(day,'VariableNames',{'Day'});
        zInv_n = tablecompile_start(day,dayIDlist,zInv(:,1),ppID,'Day','ppDay','ppDay');
        zInv_p = tablecompile_start(day,dayIDlist,zInv(:,2),ppID,'Day','ppDay','ppDay');
    else
        zInv_n = tablecompile_iter(zInv_n,dayIDlist,zInv(:,1),ppID,'Day','ppDay','ppDay');
        zInv_p = tablecompile_iter(zInv_p,dayIDlist,zInv(:,2),ppID,'Day','ppDay','ppDay');
    end
    %% get negative ICC value stats
    numNegICCs(i_subject) = sum(negICCs);
    propNegICCs(i_subject) = sum(negICCs)/numDays(i_subject);
    clear rawICC rtoZ zInv negICCs
end

%% compile all variables of interest, loading additional data as needed
% experience sampling stats and measures of experienced affect
samplingData = [numPrompts' numDays' maxDays' missingDays' numPrompts_day' mPositive' mNegative' sdPositive' sdNegative'];
samplingColumns = {'numPrompts','numDays','maxDays','missingDays','numPrompts_day','mPos','mNeg','sdPos','sdNeg'};

% text-based descriptions of events from EOD surveys
rawData_text = importdata('EOD_event_description_variables.xlsx');
textData = rawData_text.data(:,3:end);
textColumns = rawData_text.colheaders(3:end);

% other related variables
rawData_other = importdata('Other_granularity_related_variables.xlsx');
otherData = rawData_other.data(:,2:end);
otherColumns = rawData_other.colheaders(2:end);

% compile everything for given # of starting/ending days
allData = [subjectIDlist samplingData textData otherData];
allColumns = ['PPID' samplingColumns textColumns otherColumns];
allData_Table = array2table(allData,'VariableNames',allColumns);

%% convert tables to matrices
zInv_n_array = real(table2array(zInv_n(:,2:end)))';
zInv_p_array = real(table2array(zInv_p(:,2:end)))';
% impute missing data, if desired
if imputeData == 1
    zInv_n_array = fillmissing(zInv_n_array,'movmean',3,2);
    zInv_p_array = fillmissing(zInv_p_array,'movmean',3,2);
end

%% get participant-level slopes
for i_subject = 1:length(subjectIDlist)
    % negative granularity
    dailyGran_n = zInv_n_array(i_subject,:);
    dailyGran_n = dailyGran_n(~isnan(dailyGran_n));
    if dropFirstICC == 1
        dailyGran_n = dailyGran_n(2:end);
    end
    dailyGran_n_S = zscore(dailyGran_n)';
    numDays_n(i_subject) = numel(dailyGran_n);
    day_n_S = zscore(1:1:numDays_n(i_subject))';
    granCoefs_n(i_subject,:) = regress(dailyGran_n_S,[ones(numDays_n(i_subject),1) day_n_S]);
    % positive granularity
    dailyGran_p = zInv_p_array(i_subject,:);
    dailyGran_p = dailyGran_p(~isnan(dailyGran_p));
    if dropFirstICC == 1
        dailyGran_p = dailyGran_p(2:end);
    end
    dailyGran_p_S = zscore(dailyGran_p)';
    numDays_p(i_subject) = numel(dailyGran_p);
    day_p_S = zscore(1:1:numDays_p(i_subject))';
    granCoefs_p(i_subject,:) = regress(dailyGran_p_S,[ones(numDays_p(i_subject),1) day_p_S]);
    clear dailyGran_n dailyGran_n_S dailyGran_p dailyGran_p_S day_n_S day_p_S
end

%% drop participants with (above threshold) negative ICCs
if dropNegPPs == 1
    index3 = numNegICCs>2;
    granCoefs_p(index3,:) = [];
    granCoefs_n(index3,:) = [];
    allData(index3,:) = [];
end

%% run one-sample t-tests to see if slopes are different than zero
[~,p_n,~,stats_n] = ttest(granCoefs_n(:,2));
[~,p_p,~,stats_p] = ttest(granCoefs_p(:,2));
m_n = mean(granCoefs_n(:,2));
m_p = mean(granCoefs_p(:,2));
sd_n = std(granCoefs_n(:,2));
sd_p = std(granCoefs_p(:,2));
d_n = m_n/sd_n; % effect size estimate
d_p = m_p/sd_p; % effect size estimate

%% run regressions predicting change in granularity
% predictors: numDays, numPrompts_day, meanAll_prompt, affect_prompt, mPos, mNeg, restRSA_p2m
predictors = [allData(:,3) allData(:,6) allData(:,18) allData(:,27) allData(:,7) allData(:,8) allData(:,65)];
predictors_S = zscore(predictors,[],1);
granCoefs_n_S = zscore(granCoefs_n(:,2));
granCoefs_p_S = zscore(granCoefs_p(:,2));

% fit Bayesian model & get stats
priorMdl = bayeslm(size(predictors,2),'ModelType','semiconjugate');
[posteriorMdl_n,estBeta_n,~,~,~,summary_n] = estimate(priorMdl,predictors_S,granCoefs_n_S,'Display',false);
[posteriorMdl_p,estBeta_p,~,~,~,summary_p] = estimate(priorMdl,predictors_S,granCoefs_p_S,'Display',false);
CI95_n = table2array(summary_n(1:end-1,3));
CI95_p = table2array(summary_p(1:end-1,3));
positive_n = table2array(summary_n(1:end-1,4));
positive_p = table2array(summary_p(1:end-1,4)); 
n_Table = array2table(transpose(posteriorMdl_n.BetaDraws),'VariableNames',{'intercept','nDays','nPrompts','mWords','mAffectWords','mPosAffect','mNegAffect','restingRSA'});
p_Table = array2table(transpose(posteriorMdl_p.BetaDraws),'VariableNames',{'intercept','nDays','nPrompts','mWords','mAffectWords','mPosAffect','mNegAffect','restingRSA'});

% summary table
summary = [estBeta_n positive_n estBeta_p positive_p];
summary_Table = array2table(summary,'VariableNames',{'B_n','p_n','B_p','p_p'},'RowNames',{'Inter','nDays','nPrompts','mAll','mAffect','mPos','mNeg','rRSA'});

%% generate figures
% negative granularity
% plot participant-level change
figure;
toPlot = granCoefs_n(:,2);
toPlot(find(isnan(toPlot))) = 0;
bar(toPlot,'FaceColor',rgb('LightGray'));
ylim([-1 1]);
yticks([-1:0.2:1]);
ylabel({'change in negative emotional granularity'});
xlim([0 length(subjectIDlist)+1]);
xticks([1:1:length(subjectIDlist)]);
xticklabels(num2cell(subjectIDlist));
xlabel({'participant number'});

% plot posterior distributions of estimated regression coefficients
figure;
hold on;
v = violinplot(n_Table,[],'ShowData',false,'ViolinColor',rgb('OrangeRed'));
starHeight = .90;
sigFactors = find(positive_n>.95 | positive_n<.05);
trendFactors = setdiff(find(positive_n>.90 | positive_n<.10),sigFactors);
plot(sigFactors,starHeight*ones(numel(sigFactors),1),'k*');
plot(trendFactors,starHeight*ones(numel(trendFactors),1),'k+');
xtickangle(30);
ylim([-1 1]);
ylabel({'posterior mean (estimated ?)'});
hline(0,'k:');
hold off;

% positive granularity
% plot participant-level change
figure;
toPlot = granCoefs_p(:,2);
toPlot(find(isnan(toPlot))) = 0;
bar(toPlot,'FaceColor',rgb('LightGray'));
ylim([-1 1]);
yticks([-1:0.2:1]);
ylabel({'change in positive emotional granularity'});
xlim([0 length(subjectIDlist)+1]);
xticks([1:1:length(subjectIDlist)]);
xticklabels(num2cell(subjectIDlist));
xlabel({'participant number'});

% plot posterior distributions of estimated regression coefficients
figure;
hold on;
v = violinplot(p_Table,[],'ShowData',false,'ViolinColor',rgb('MediumBlue'));
starHeight = .90;
sigFactors = find(positive_p>.95 | positive_p<.05);
trendFactors = setdiff(find(positive_p>.90 | positive_p<.10),sigFactors);
plot(sigFactors,starHeight*ones(numel(sigFactors),1),'k*');
plot(trendFactors,starHeight*ones(numel(trendFactors),1),'k+');
xtickangle(30);
ylim([-1 1]);
ylabel({'posterior mean (estimated ?)'});
hline(0,'k:');
hold off;