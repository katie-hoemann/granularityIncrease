clear all;
clc;

%% set parameters
useFiltered = 0; % set to 1 to filter out excluded subjects, days with < 6 completed prompts and late surveys
upperBoundDays = 4; % set to the maximum number of days to include as starting/ending days

%% load data file of emotion intensity ratings, along with word file that includes raw norms
rawData_ratings = importdata('18ARIEOD_daily_filtered.xlsx');
if useFiltered == 1
    ratingsData = rawData_ratings.data.NoLateSurveys; 
else
    ratingsData = rawData_ratings.data.SelectColumns;
end
words = readtable('words18.csv'); 
wordList = rawData_ratings.colheaders.NoLateSurveys(5:end)';  % grab sampled words from top row of data file

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

%% perform calculations across different # of starting/ending days
for i_days = 1:upperBoundDays
    %% grab data for each subject and run through calculations
    subjectIDlist = unique(ratingsData(:,1)); % grab subject IDs from first column of data file
    for i_subject = 1:length(subjectIDlist)
        subjectData = [];
        subjectID = subjectIDlist(i_subject);
        index1 = find(ratingsData(:,1)==subjectID);
        subjectData = ratingsData(index1,2:end);
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
        %% run through calculations for starting # experience sampling days
        start = subjectData(:,1) <= i_days;
        subjectData_start = subjectData(start==1,4:end);
        % compute ICCs - positive, negative, valence average; combinations of valence x arousal
        rawICC_start(i_subject,1) = ICC(subjectData_start(:,(labels(:,1)==0)),'A-k'); % negative valence ICC
        rawICC_start(i_subject,2) = ICC(subjectData_start(:,(labels(:,1)==1)),'A-k'); % positive valence ICC
        rawICC_start(i_subject,3) = (rawICC_start(i_subject,1)+rawICC_start(i_subject,2))/2; % valence mean ICC
        % calculate granularity from ICCs
        gran_start(i_subject,1) = 1-rawICC_start(i_subject,1);
        gran_start(i_subject,2) = 1-rawICC_start(i_subject,2);
        gran_start(i_subject,3) = 1-rawICC_start(i_subject,3);
        % Fisher transform ICCs (noted as zICCs)
        rtoZ_start(i_subject,1) = 0.5*log((1+rawICC_start(i_subject,1))/(1-rawICC_start(i_subject,1)));
        rtoZ_start(i_subject,2) = 0.5*log((1+rawICC_start(i_subject,2))/(1-rawICC_start(i_subject,2)));
        rtoZ_start(i_subject,3) = 0.5*log((1+rawICC_start(i_subject,3))/(1-rawICC_start(i_subject,3)));
        % invert zICCs for intuitive directionality
        zInv_start(i_subject,1) = rtoZ_start(i_subject,1)*-1;
        zInv_start(i_subject,2) = rtoZ_start(i_subject,2)*-1;
        zInv_start(i_subject,3) = rtoZ_start(i_subject,3)*-1;
        %% run through calculations for ending # experience sampling days
        ending = subjectData(:,1) >= (max(subjectData(:,1))-(i_days-1));
        subjectData_end = subjectData(ending==1,4:end);
        % compute ICCs - positive, negative, valence average; combinations of valence x arousal
        rawICC_end(i_subject,1) = ICC(subjectData_end(:,(labels(:,1)==0)),'A-k'); % negative valence ICC
        rawICC_end(i_subject,2) = ICC(subjectData_end(:,(labels(:,1)==1)),'A-k'); % positive valence ICC
        rawICC_end(i_subject,3) = (rawICC_end(i_subject,1)+rawICC_end(i_subject,2))/2; % valence mean ICC
        % calculate granularity from ICCs
        gran_end(i_subject,1) = 1-rawICC_end(i_subject,1);
        gran_end(i_subject,2) = 1-rawICC_end(i_subject,2);
        gran_end(i_subject,3) = 1-rawICC_end(i_subject,3);
        % Fisher transform ICCs (noted as zICCs)
        rtoZ_end(i_subject,1) = 0.5*log((1+rawICC_end(i_subject,1))/(1-rawICC_end(i_subject,1)));
        rtoZ_end(i_subject,2) = 0.5*log((1+rawICC_end(i_subject,2))/(1-rawICC_end(i_subject,2)));
        rtoZ_end(i_subject,3) = 0.5*log((1+rawICC_end(i_subject,3))/(1-rawICC_end(i_subject,3)));
        % invert zICCs for intuitive directionality
        zInv_end(i_subject,1) = rtoZ_end(i_subject,1)*-1;
        zInv_end(i_subject,2) = rtoZ_end(i_subject,2)*-1;
        zInv_end(i_subject,3) = rtoZ_end(i_subject,3)*-1;
    end

    %% compile all variables of interest, loading additional data as needed
    % experience sampling stats and measures of experienced affect
    samplingData = [numPrompts' numDays' maxDays' missingDays' numPrompts_day' mPositive' mNegative' sdPositive' sdNegative'];
    samplingColumns = {'numPrompts','numDays','maxDays','missingDays','numPrompts_day','mPos','mNeg','sdPos','sdNeg'};
    if useFiltered == 0
        excludedIDs = [4; 8; 16; 31; 35; 42; 44; 48; 53; 60; 65]; % subjects excluded from N=50 analyses
        excludedSubjects = zeros(length(subjectIDlist),1);
        for i_row = 1:length(subjectIDlist)
            if any(excludedIDs == subjectIDlist(i_row))
                excludedSubjects(i_row) = 1;
            end
        end
        samplingData(excludedSubjects==1,:) = []; % drop excluded subjects
        subjectIDlist(excludedSubjects==1) = []; 
        zInv_start(excludedSubjects==1,:) = [];
        zInv_end(excludedSubjects==1,:) = [];
    end

    % text-based descriptions of events from EOD surveys
    rawData_text = importdata('EOD_event_description_variables.xlsx');
    textData = rawData_text.data(:,3:end);
    textColumns = rawData_text.colheaders(3:end);

    % other related variables
    rawData_other = importdata('Other_granularity_related_variables.xlsx');
    otherData = rawData_other.data(:,2:end);
    otherColumns = rawData_other.colheaders(2:end);

    % column headers for starting and ending granularity values
    startColumns = {'NegV_start','PosV_start','MeanV_start'};
    endColumns = {'NegV_end','PosV_end','MeanV_end'};

    % compile everything
    allData = real([subjectIDlist samplingData textData otherData zInv_start(:,:) zInv_end(:,:)]);
    allColumns = ['PPID' samplingColumns textColumns otherColumns startColumns endColumns];
    allData_Table = array2table(allData,'VariableNames',allColumns);

    %% test whether granularity improves between start and end
    [~,p_m(i_days),ci_m(:,i_days),stats_m] = ttest(allData(:,73),allData(:,70)); % overall granularity
    [~,p_n(i_days),ci_n(:,i_days),stats_n] = ttest(allData(:,71),allData(:,68)); % negative granularity
    [~,p_p(i_days),ci_p(:,i_days),stats_p] = ttest(allData(:,72),allData(:,69)); % positive granularity
    t_m(i_days) = stats_m.tstat;
    t_n(i_days) = stats_n.tstat;
    t_p(i_days) = stats_p.tstat;
    df_m(i_days) = stats_m.df;
    df_n(i_days) = stats_n.df;
    df_p(i_days) = stats_p.df;
    p_corrected_m(i_days) = p_m(i_days)*upperBoundDays;
    p_corrected_n(i_days) = p_n(i_days)*upperBoundDays;
    p_corrected_p(i_days) = p_p(i_days)*upperBoundDays;
    % run correlations between start and end values
    r_m(i_days) = corr(allData(:,73),allData(:,70),'rows','complete'); % overall granularity
    r_n(i_days) = corr(allData(:,71),allData(:,68),'rows','complete'); % negative granularity
    r_p(i_days) = corr(allData(:,72),allData(:,69),'rows','complete'); % positive granularity
    % get descriptive stats
    startGranM_m(i_days) = nanmean(allData(:,70));
    endGranM_m(i_days) = nanmean(allData(:,73));
    startGranM_n(i_days) = nanmean(allData(:,68));
    endGranM_n(i_days) = nanmean(allData(:,71));
    startGranM_p(i_days) = nanmean(allData(:,69));
    endGranM_p(i_days) = nanmean(allData(:,72));
    startGranSD_m(i_days) = nanstd(allData(:,70));
    endGranSD_m(i_days) = nanstd(allData(:,73));
    startGranSD_n(i_days) = nanstd(allData(:,68));
    endGranSD_n(i_days) = nanstd(allData(:,71));
    startGranSD_p(i_days) = nanstd(allData(:,69));
    endGranSD_p(i_days) = nanstd(allData(:,72));
    % get effect size estimates
    d_num_m = endGranM_m(i_days)-startGranM_m(i_days);
    d_num_n = endGranM_n(i_days)-startGranM_n(i_days);
    d_num_p = endGranM_p(i_days)-startGranM_p(i_days);
    d_denom_m = stats_m.sd;
    d_denom_n = stats_n.sd;
    d_denom_p = stats_p.sd;
    d_m(i_days) = d_num_m/d_denom_m;
    d_n(i_days) = d_num_n/d_denom_n;
    d_p(i_days) = d_num_p/d_denom_p;
    
    % derive change scores
    granChange_m = allData(:,73)-allData(:,70); % overall granularity
    granChange_n = allData(:,71)-allData(:,68); % negative granularity
    granChange_p = allData(:,72)-allData(:,69); % positive granularity
    
    % get change score statistics
    min_granChange_m(i_days) = min(granChange_m);
    min_granChange_n(i_days) = min(granChange_n);
    min_granChange_p(i_days) = min(granChange_p);
    max_granChange_m(i_days) = max(granChange_m);
    max_granChange_n(i_days) = max(granChange_n);
    max_granChange_p(i_days) = max(granChange_p);
    
    % plot participant-level change in overall emotional granularity
    figure;
    toPlot = granChange_m;
    toPlot(find(isnan(toPlot))) = 0;
    bar(toPlot,'FaceColor',rgb('LightGray'));
    ylim([-2 2]);
    yticks([-2:0.2:2]);
    ylabel({'change in emotional granularity'});
    xlim([0 length(subjectIDlist)+1]);
    xticks([1:1:length(subjectIDlist)]);
    xticklabels(num2cell(subjectIDlist));
    xlabel({'participant number'});
    title({['day ' num2str(i_days)]});

    %% run regressions predicting change in granularity
    % standardize change scores
    granChange_m = (granChange_m - nanmean(granChange_m))/nanstd(granChange_m); 
    granChange_n = (granChange_n - nanmean(granChange_n))/nanstd(granChange_n); 
    granChange_p = (granChange_p - nanmean(granChange_p))/nanstd(granChange_p);

    % predictors: numDays, numPrompts_day, meanAll_prompt, affect_prompt, mPos, mNeg, restRSA_p2m
    predictors = [allData(:,3) allData(:,6) allData(:,18) allData(:,27) allData(:,7) allData(:,8) allData(:,65)];
    for i_variable = 1:size(predictors,2) % ensure predictors are normally distributed and correct if not
        [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(predictors(:,i_variable));
        if h(i_variable) == 1
            fracRank = FractionalRankings(predictors(:,i_variable));
            fracRank = fracRank/(max(fracRank));
            fracRank(fracRank == 1) = .999;
            predictors(:,i_variable) = norminv(fracRank,mean(predictors(:,i_variable)),std(predictors(:,i_variable)));
        end
        predictors(:,i_variable) = zscore(predictors(:,i_variable)); % standardize variable
    end

    % fit Bayesian model & get stats
    priorMdl = bayeslm(size(predictors,2));
    [~,estBeta_m(:,i_days),~,~,~,posteriorMdl_m] = estimate(priorMdl,predictors,granChange_m,'Display',false);
    [~,estBeta_n(:,i_days),~,~,~,posteriorMdl_n] = estimate(priorMdl,predictors,granChange_n,'Display',false);
    [~,estBeta_p(:,i_days),~,~,~,posteriorMdl_p] = estimate(priorMdl,predictors,granChange_p,'Display',false);
    numObs_m(i_days) = length(granChange_m)-sum(isnan(granChange_m));
    numObs_n(i_days) = length(granChange_n)-sum(isnan(granChange_n));
    numObs_p(i_days) = length(granChange_p)-sum(isnan(granChange_p));
    CI95_m(:,:,i_days) = table2array(posteriorMdl_m(1:end-1,3));
    CI95_n(:,:,i_days) = table2array(posteriorMdl_n(1:end-1,3));
    CI95_p(:,:,i_days) = table2array(posteriorMdl_p(1:end-1,3));
    positive_m(:,i_days) = table2array(posteriorMdl_m(1:end-1,4));
    positive_n(:,i_days) = table2array(posteriorMdl_n(1:end-1,4));
    positive_p(:,i_days) = table2array(posteriorMdl_p(1:end-1,4));    
end

%% get summary statistics across all iterations
medianBeta_m = median(estBeta_m,2);
medianBeta_n = median(estBeta_n,2);
medianBeta_p = median(estBeta_p,2);
medianPositive_m = median(positive_m,2);
medianPositive_n = median(positive_n,2);
medianPositive_p = median(positive_p,2);
for i_variable = 1:size(positive_m,1)
    numberSig_m(i_variable) = sum(positive_m(i_variable,:)<.05)+sum(positive_m(i_variable,:)>.95);
    numberSig_n(i_variable) = sum(positive_n(i_variable,:)<.05)+sum(positive_n(i_variable,:)>.95);
    numberSig_p(i_variable) = sum(positive_p(i_variable,:)<.05)+sum(positive_p(i_variable,:)>.95);
    numberTrend_m(i_variable) = sum(positive_m(i_variable,:)<.10)+sum(positive_m(i_variable,:)>.90);
    numberTrend_n(i_variable) = sum(positive_n(i_variable,:)<.10)+sum(positive_n(i_variable,:)>.90);
    numberTrend_p(i_variable) = sum(positive_p(i_variable,:)<.10)+sum(positive_p(i_variable,:)>.90);
end

%% create layered figure of results across all # starting/ending days
% overall granularity
figure;
hold on;
for i_days = 1:upperBoundDays % plot 95% credible intervals with transparency
    for i_variable = 1:size(estBeta_m,1)
        patch('XData',[(i_variable-.25) (i_variable+.25) (i_variable+.25) (i_variable-.25)],...
            'YData',[CI95_m(i_variable,1,i_days) CI95_m(i_variable,1,i_days) CI95_m(i_variable,2,i_days) CI95_m(i_variable,2,i_days)],...
            'FaceColor',rgb('LightGray'),'FaceAlpha',0.3,'EdgeColor',rgb('DarkGray'));
    end
end
boxplot(transpose(estBeta_m));
xticklabels({'intercept','nDays','nPrompts/day','mWords/prompt','mAffect/prompt','mPosAffect','mNegAffect','restingRSA'});
xtickangle(30);
ylim([-1 1]);
ylabel({'posterior mean (estimated beta)'});
hline(0,'k:');
%title({'change in overall emotional granularity'});
hold off;

% negative granularity
figure;
hold on;
for i_days = 1:upperBoundDays % plot 95% credible intervals with transparency
    for i_variable = 1:size(estBeta_n,1)
        patch('XData',[(i_variable-.25) (i_variable+.25) (i_variable+.25) (i_variable-.25)],...
            'YData',[CI95_n(i_variable,1,i_days) CI95_n(i_variable,1,i_days) CI95_n(i_variable,2,i_days) CI95_n(i_variable,2,i_days)],...
            'FaceColor',rgb('LightGray'),'FaceAlpha',0.3,'EdgeColor',rgb('DarkGray'));
    end
end
boxplot(transpose(estBeta_n));
xticklabels({'intercept','nDays','nPrompts/day','mWords/prompt','mAffect/prompt','mPosAffect','mNegAffect','restingRSA'});
xtickangle(30);
ylim([-1 1]);
ylabel({'posterior mean (estimated beta)'});
hline(0,'k:');
%title({'change in negative emotional granularity'});
hold off;

% positive granularity
figure;
hold on;
for i_days = 1:upperBoundDays % plot 95% credible intervals with transparency
    for i_variable = 1:size(estBeta_p,1)
        patch('XData',[(i_variable-.25) (i_variable+.25) (i_variable+.25) (i_variable-.25)],...
            'YData',[CI95_p(i_variable,1,i_days) CI95_p(i_variable,1,i_days) CI95_p(i_variable,2,i_days) CI95_p(i_variable,2,i_days)],...
            'FaceColor',rgb('LightGray'),'FaceAlpha',0.3,'EdgeColor',rgb('DarkGray'));
    end
end
boxplot(transpose(estBeta_p));
xticklabels({'intercept','nDays','nPrompts/day','mWords/prompt','mAffect/prompt','mPosAffect','mNegAffect','restingRSA'});
xtickangle(30);
ylim([-1 1]);
ylabel({'posterior mean (estimated beta)'});
hline(0,'k:');
%title({'change in positive emotional granularity'});
hold off;





