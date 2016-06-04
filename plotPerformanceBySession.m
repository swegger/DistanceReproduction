function out = plotPerformanceBySession(slist,varargin)
%% plotPerformanceBySession
%
%   Plots behavioral performance and parameters fit to the BLS model by
%   session number for a list of subjects.
%
%%


%% Defaults
PlotOpts_default.titles = {'RS1G','RS2G'};
PlotOpts_default.RelativeFigSize = [1/5 1/2 3/5 3/5];
PlotOpts_default.colors = [0 0 1; 1 0 0; 0.6 0.6 0.6; 0 0 0];
PlotOpts_default.nbins = 20;
PlotOpts_default.Hold = 'On';

%% Parse input
Parser = inputParser;

addRequired(Parser,'slist')     % List of subjects to analyze
addParameter(Parser,'runs',NaN)     % Runs to analyze
addParameter(Parser,'interval_N',1:2)    % Trial types to analyze
addParameter(Parser,'sampleIntervals',[600 800 1000])   % Sample intervals to measure sensitivity to
addParameter(Parser,'outlier',Inf)      % Number of standard deviations away from the mean to exclude data as outlier
addParameter(Parser,'trialWin',[1 Inf]) % Trials in each run to analyze ([start# end#])
addParameter(Parser,'PlotTypes',{'Errors','modelParams'})      % Data to plot
addParameter(Parser,'PlotOpts',PlotOpts_default)    % Plotting options
addParameter(Parser,'Save','No')           % Saving options
addParameter(Parser,'MinMaxTp',[-Inf Inf])      % Minimum and maximum values of t_p to keep
addParameter(Parser,'ConflictType','equal')     % For experiments with cue conflict
addParameter(Parser,'DiffTolerance',2/60)       % Tolerance for difference in sample times before calling it conflict
addParameter(Parser,'LapseRejection','none')
addParameter(Parser,'ModelParams',NaN)

parse(Parser,slist,varargin{:})

slist = Parser.Results.slist;
runs = Parser.Results.runs;
interval_N = Parser.Results.interval_N;
sampleIntervals = Parser.Results.sampleIntervals;
outlier = Parser.Results.outlier;
trialWin = Parser.Results.trialWin;
PlotTypes = Parser.Results.PlotTypes;
PlotOpts = Parser.Results.PlotOpts;
Save = Parser.Results.Save;
MinMaxTp = Parser.Results.MinMaxTp;
ConflictType = Parser.Results.ConflictType;
DiffTolerance = Parser.Results.DiffTolerance;
LapseRejection = Parser.Results.LapseRejection;
ModelParams = Parser.Results.ModelParams;


% Global variables
for i = 1:length(interval_N)
    intervalNstr{i} = num2str(interval_N(i));
end

% Open global figures
if strcmp(PlotOpts.Hold,'On')
    % Set up errors figure for all subjects
    if any(strcmp('Errors',PlotTypes))
        fHandle.Errors = figure('Name',[' Errors']);
        aHandle.Errors(1) = subplot(1,3,1);
        aHandle.Errors(2) = subplot(1,3,2);
        aHandle.Errors(3) = subplot(1,3,3);
    end
    
    % Set up BiasVariance figure for all subjects
    if any(strcmp('BiasVariance',PlotTypes))
        fHandle.BiasVariance = figure('Name',[' Bias and Variance']);
        aHandle.BiasVariance(1) = subplot(1,2,1);
        aHandle.BiasVariance(2) = subplot(1,2,2);
    end
    
    % Set up RMSE figure for all subjects
    if any(strcmp('RMSE',PlotTypes))
        fHandle.RMSE = figure('Name',[slist{i} ' RMSE']);
    end
    
    if any(strcmp('modelParams',PlotTypes))
        fHandle.modelParams = figure('Name',[' Model Parameters']);
        aHandle.modelParams(1) = subplot(1,4,1);
        aHandle.modelParams(2) = subplot(1,4,2);
        aHandle.modelParams(3) = subplot(1,4,3);
        aHandle.modelParams(4) = subplot(1,4,4);
    end
    
    if any(strcmp('d''',PlotTypes))
        fHandle.dprime = figure('Name',[' d''']);
        subplotind = 0;
        temp = meshgrid(sampleIntervals);
        subplots = sum(sum(triu(temp) == 0));
        for j = 1:length(sampleIntervals)
            for k = 1:length(sampleIntervals)
                if j < k
                    subplotind = subplotind + 1;
                    aHandle.dprime(subplotind) = subplot(1,subplots,subplotind);
                end
            end
        end
    end
end

% Cycle through each subject and plot 
for i = 1:length(slist)
    
    % Load subject data
    d = LoadMWorksProj('RS2G_psychophysics',slist{i},'SaveLocation','LocalHeader');
    
    % Find the runs to analyze
    if isnan(runs)
        runstemp = 1:d.runs;
    else
        runstemp = runs;
        if any(runstemp > d.runs)
            warning('Asked for a run number outside of the runs analyzed. Analyzing up to the last run.')
            runstemp = runstemp(runstemp <= d.runs);
        end
    end
    
    % Errors plotting
    if any(strcmp('Errors',PlotTypes))
        ApplyOffsetCorrection.On = 'Yes';
        ApplyOffsetCorrection.Source = 'SessionAverage';
        for run = runstemp
            [biassquared{i}(run,:), variance{i}(run,:), rmse{i}(run,:)] = RS2G_psychophysicsErrorMetrics(d,'runs',run,'interval_N',interval_N,'outlier',outlier,'trialWin',trialWin,...
                                                                          'MinMaxTp',MinMaxTp,'ConflictType',ConflictType,'DiffTolerance',DiffTolerance,'ApplyOffsetCorrection',ApplyOffsetCorrection,...
                                                                          'LapseRejection',LapseRejection,'ModelParams',ModelParams);
        end
        
        % Plot the data
        if strcmp(PlotOpts.Hold,'On')
            figure(fHandle.Errors)
        elseif strcmp(PlotOpts.Hold,'Off')
            figure('Name',[slist{i} ' Errors'])
        else
            error(['Plot Hold option ' PlotOpts.Hold ' not recognized!'])
        end
        
        if strcmp(PlotOpts.Hold,'On')
            axes(aHandle.Errors(1))
        elseif strcmp(PlotOpts.Hold,'Off')
            subplot(1,3,1)
        end
        h1 = plot(runstemp,rmse{i},'.','MarkerSize',10);
        hold on
        h2 = plot(runstemp,rmse{i},'LineWidth',2);
        for j = 1:length(interval_N)
            set(h1(j),'Color',PlotOpts.colors(j,:))
            set(h2(j),'Color',PlotOpts.colors(j,:))
        end
        ax(1,:) = axis;
        xlabel('Run')
        ylabel('RMSE (ms)')
        
        if strcmp(PlotOpts.Hold,'On')
            axes(aHandle.Errors(2))
        elseif strcmp(PlotOpts.Hold,'Off')
            subplot(1,3,2)
        end
        h1 = plot(runstemp,sqrt(biassquared{i}),'.','MarkerSize',10);
        hold on
        h2 = plot(runstemp,sqrt(biassquared{i}),'LineWidth',2);
        for j = 1:length(interval_N)
            set(h1(j),'Color',PlotOpts.colors(j,:))
            set(h2(j),'Color',PlotOpts.colors(j,:))
        end
        ax(2,:) = axis;
        xlabel('Run')
        ylabel('Bias (ms)')
        
        if strcmp(PlotOpts.Hold,'On')
            axes(aHandle.Errors(3))
        elseif strcmp(PlotOpts.Hold,'Off')
            subplot(1,3,3)
        end
        h1 = plot(runstemp,sqrt(variance{i}),'.','MarkerSize',10);
        hold on
        h2 = plot(runstemp,sqrt(variance{i}),'LineWidth',2);
        for j = 1:length(interval_N)
            set(h1(j),'Color',PlotOpts.colors(j,:))
            set(h2(j),'Color',PlotOpts.colors(j,:))
        end
        ax(3,:) = axis;
        xlabel('Run')
        ylabel('sqrt(VAR)')
        legend(intervalNstr)
        
        for j = 1:3
            if strcmp(PlotOpts.Hold,'On')
                axes(aHandle.Errors(j))
            elseif strcmp(PlotOpts.Hold,'Off')
                subplot(1,3,j)
            end
            axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
        end
    end
    
    % BiasVariance plotting
    if any(strcmp('BiasVariance',PlotTypes))
        ApplyOffsetCorrection.On = 'Yes';
        ApplyOffsetCorrection.Source = 'SessionAverage';
        for run = runstemp
            [biassquared{i}(run,:), variance{i}(run,:)] = RS2G_psychophysicsErrorMetrics(d,'runs',run,'interval_N',interval_N,'outlier',outlier,'trialWin',trialWin,'MinMaxTp',MinMaxTp,'ConflictType',ConflictType,'DiffTolerance',DiffTolerance,'ApplyOffsetCorrection',ApplyOffsetCorrection);
        end
        
        
        % Plot the data
        if strcmp(PlotOpts.Hold,'On')
            figure(fHandle.BiasVariance)
        elseif strcmp(PlotOpts.Hold,'Off')
            figure('Name',[slist{i} ' Bias and Variance'])
        else
            error(['Plot Hold option ' PlotOpts.Hold ' not recognized!'])
        end
        if strcmp(PlotOpts.Hold,'On')
            axes(aHandle.BiasVariance(1))
        elseif strcmp(PlotOpts.Hold,'Off')
            subplot(1,2,1)
        end
        h1 = plot(runstemp,sqrt(biassquared{i}),'.','MarkerSize',10);
        hold on
        h2 = plot(runstemp,sqrt(biassquared{i}),'LineWidth',2);
        for j = 1:length(interval_N)
            set(h1(j),'Color',PlotOpts.colors(j,:))
            set(h2(j),'Color',PlotOpts.colors(j,:))
        end
        xlabel('Run')
        ylabel('Bias (ms)')
        
        if strcmp(PlotOpts.Hold,'On')
            axes(aHandle.BiasVariance(2))
        elseif strcmp(PlotOpts.Hold,'Off')
            subplot(1,2,2)
        end
        h1 = plot(runstemp,sqrt(variance{i}),'.','MarkerSize',10);
        hold on
        h2 = plot(runstemp,sqrt(variance{i}),'LineWidth',2);
        for j = 1:length(interval_N)
            set(h1(j),'Color',PlotOpts.colors(j,:))
            set(h2(j),'Color',PlotOpts.colors(j,:))
        end
        xlabel('Run')
        ylabel('sqrt(VAR)')
        legend(intervalNstr)
    end
    
    % RMSE plotting
    if any(strcmp('RMSE',PlotTypes))
        ApplyOffsetCorrection.On = 'Yes';
        ApplyOffsetCorrection.Source = 'SessionAverage';
        for run = runstemp
            [~, ~, rmse{i}(run,:)] = RS2G_psychophysicsErrorMetrics(d,'runs',run,'interval_N',interval_N,'outlier',outlier,'trialWin',trialWin,'MinMaxTp',MinMaxTp,'ConflictType',ConflictType,'DiffTolerance',DiffTolerance,'ApplyOffsetCorrection',ApplyOffsetCorrection);
        end
        
        % Plot the data
        if strcmp(PlotOpts.Hold,'On')
            figure(fHandle.RMSE)
        elseif strcmp(PlotOpts.Hold,'Off')
            figure('Name',[slist{i} ' RMSE'])
        else
            error(['Plot Hold option ' PlotOpts.Hold ' not recognized!'])
        end
        h1 = plot(runstemp,rmse{i},'.','MarkerSize',10);
        hold on
        h2 = plot(runstemp,rmse{i},'LineWidth',2);
        for j = 1:length(interval_N)
            set(h1(j),'Color',PlotOpts.colors(j,:))
            set(h2(j),'Color',PlotOpts.colors(j,:))
        end
        xlabel('Run')
        ylabel('RMSE (ms)')
        legend(intervalNstr)
    end
    
    % Model parameters
    if any(strcmp('modelParams',PlotTypes))
        for run = runstemp
            wm{i}(run) = d.modelParams{run}.wm;
            wp{i}(run) = d.modelParams{run}.wp;
            b{i}(run) = d.modelParams{run}.b;
            if isfield(d.modelParams{run},'lapse')
                lapse{i}(run) = d.modelParams{run}.lapse;
            else
                lapse{i}(run) = NaN;
            end
        end
        
        % Plot the data
        if strcmp(PlotOpts.Hold,'On')
            figure(fHandle.modelParams)
        elseif strcmp(PlotOpts.Hold,'Off')
            figure('Name',[slist{i} ' Model Parameters'])
        else
            error(['Plot Hold option ' PlotOpts.Hold ' not recognized!'])
        end
        if strcmp(PlotOpts.Hold,'On')
            axes(aHandle.modelParams(1))
        elseif strcmp(PlotOpts.Hold,'Off')
            subplot(1,4,1)
        end
        h1 = plot(runstemp,wm{i},'.','MarkerSize',10,'Color',[0 0 0]);
        hold on
        h2 = plot(runstemp,wm{i},'LineWidth',2,'Color',[0 0 0]);
        xlabel('Run')
        ylabel('wm')
        
        if strcmp(PlotOpts.Hold,'On')
            axes(aHandle.modelParams(2))
        elseif strcmp(PlotOpts.Hold,'Off')
            subplot(1,4,2)
        end
        h1 = plot(runstemp,wp{i},'.','MarkerSize',10,'Color',[0 0 0]);
        hold on
        h2 = plot(runstemp,wp{i},'LineWidth',2,'Color',[0 0 0]);
        xlabel('Run')
        ylabel('wp')
        
        if strcmp(PlotOpts.Hold,'On')
            axes(aHandle.modelParams(3))
        elseif strcmp(PlotOpts.Hold,'Off')
            subplot(1,4,3)
        end
        h1 = plot(runstemp,b{i},'.','MarkerSize',10,'Color',[0 0 0]);
        hold on
        h2 = plot(runstemp,b{i},'LineWidth',2,'Color',[0 0 0]);
        xlabel('Run')
        ylabel('b')
        
        if strcmp(PlotOpts.Hold,'On')
            axes(aHandle.modelParams(4))
        elseif strcmp(PlotOpts.Hold,'Off')
            subplot(1,4,4)
        end
        h1 = plot(runstemp,lapse{i},'.','MarkerSize',10,'Color',[0 0 0]);
        hold on
        h2 = plot(runstemp,lapse{i},'LineWidth',2,'Color',[0 0 0]);
        xlabel('Run')
        ylabel('prob(lapse)')
    end
    
    % d'
    if any(strcmp('d''',PlotTypes))
        for run = runstemp
            [dprime{i}(:,:,:,run)] = RS2G_psychophysics_dprime(d,'runs',run,'interval_N',interval_N,'outlier',outlier,'trialWin',trialWin,'MinMaxTp',MinMaxTp,'ConflictType',ConflictType,'DiffTolerance',DiffTolerance,'sampleIntervals',sampleIntervals);
        end
        
        % Plot the data
        if strcmp(PlotOpts.Hold,'On')
            figure(fHandle.dprime)
        elseif strcmp(PlotOpts.Hold,'Off')
            figure('Name',[slist{i} ' d'''])
        else
            error(['Plot Hold option ' PlotOpts.Hold ' not recognized!'])
        end
        subplotind = 0;
        temp = meshgrid(sampleIntervals);
        subplots = sum(sum(triu(temp) == 0));
        for j = 1:length(sampleIntervals)
            for k = 1:length(sampleIntervals)
                if j < k
                    subplotind = subplotind + 1;
                    if strcmp(PlotOpts.Hold,'On')
                        axes(aHandle.dprime(subplotind))
                    elseif strcmp(PlotOpts.Hold,'Off')
                        subplot(1,subplots,subplotind)
                    end
                    for jj = 1:length(interval_N)
                        h1 = plot(runstemp,squeeze(dprime{i}(j,k,jj,:)),'.','MarkerSize',10);
                        hold on
                        h2 = plot(runstemp,squeeze(dprime{i}(j,k,jj,:)),'LineWidth',2);
                        set(h1,'Color',PlotOpts.colors(jj,:))
                        set(h2,'Color',PlotOpts.colors(jj,:))
                        title(['ts = ' num2str(sampleIntervals(j)) ' and ' num2str(sampleIntervals(k))])
                    end
                end
            end
        end
        legend(intervalNstr)
    end
    
end

    
out = [];