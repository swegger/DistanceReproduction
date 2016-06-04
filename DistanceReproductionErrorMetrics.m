function [biassquared, variance, rmse, BIASs, VARs] = DistanceReproductionErrorMetrics(d,varargin)
%% DistanceReproductionErrorMetrics
%
%   [biassquared, variance, rmse, BIASs, VARs] = DistanceReproductionErrorMetrics(varargin)   
%
%%

%% Defaults
bootstrap_default.nbootstraps = NaN;
bootstrap_default.nsamps = NaN;
PlotOpts_default.titles = {'S1','S2'};
PlotOpts_default.RelativeFigSize = [1/5 1/2 2/5 2/5];
PlotOpts_default.colors = [0 0 1; 1 0 0; 0.6 0.6 0.6; 0 0 0];
ApplyOffsetCorrection_default.On = 'Yes';
ApplyOffsetCorrection_default.Source = 'SessionAverage';

%% Parse input
Parser = inputParser;

addRequired(Parser,'d')     % Data structure
addParameter(Parser,'runs',NaN)     % Runs to analyze
addParameter(Parser,'Distance_N',1:2)    % Trial types to analyze
addParameter(Parser,'outlier',Inf)      % Number of standard deviations away from the mean to exclude data as outlier
addParameter(Parser,'trialWin',[1 Inf]) % Trials in each run to analyze ([start# end#])
addParameter(Parser,'Plot','none')      % Data to plot
addParameter(Parser,'PlotOpts',PlotOpts_default)    % Plotting options
addParameter(Parser,'Save','No')           % Saving options
addParameter(Parser,'Bootstrap',bootstrap_default)  % Bootstrap of bias/variance
addParameter(Parser,'MinMaxDp',[-Inf Inf])      % Minimum and maximum values of t_p to keep
addParameter(Parser,'ConflictType','equal')     % For experiments with cue conflict
addParameter(Parser,'DiffTolerance',2/60)       % Tolerance for difference in sample times before calling it conflict
addParameter(Parser,'ApplyOffsetCorrection',ApplyOffsetCorrection_default)  % Apply offset correction if measure of 
addParameter(Parser,'LapseRejection','none')
addParameter(Parser,'ModelParams',NaN)

parse(Parser,d,varargin{:})

d = Parser.Results.d;
runs = Parser.Results.runs;
Distance_N = Parser.Results.Distance_N;
outlier = Parser.Results.outlier;
trialWin = Parser.Results.trialWin;
Plot = Parser.Results.Plot;
PlotOpts = Parser.Results.PlotOpts;
Save = Parser.Results.Save;
Bootstrap = Parser.Results.Bootstrap;
MinMaxDp = Parser.Results.MinMaxDp;
ConflictType = Parser.Results.ConflictType;
DiffTolerance = Parser.Results.DiffTolerance;
ApplyOffsetCorrection = Parser.Results.ApplyOffsetCorrection;
LapseRejection = Parser.Results.LapseRejection;
ModelParams = Parser.Results.ModelParams;

% Check to see if run information was provided
if isnan(runs)
    runs = 1:d.runs;            % Defaults to pool data from all runs
end

% Set m to Distance_Ns
m = 1:length(Distance_N);

% Pull out bootstrap parameters
if any(isnan([Bootstrap.nbootstraps Bootstrap.nsamps]))
    bootflg = 0;
else
    bootflg = 1;
    nbootstraps = Bootstrap.nbootstraps;
    nsamps = Bootstrap.nsamps;
end

if isnan(ModelParams)
    for j = runs
        wm(j) = d.modelParams{j}.wm;
        wp(j) = d.modelParams{j}.wp;
        b(j) = d.modelParams{j}.b;
        if isfield(d.modelParams{j},'lapse')
            lapse(j) = d.modelParams{j}.lapse;
        else
            lapse(j) = NaN;
        end
    end
    WM = mean(wm);
    WP = mean(wp);
    B = mean(b);
    LAPSE = nanmean(lapse);
    if isnan(LAPSE)
        LAPSE = 0;
    end
end

%% Analyze data

for i = m
    % Grab the appropriate data
    [ds1{i}, ds2{i}, ds{i}, dp{i}, ~, ~, Trials{i}, correct{i}, ~, ~, ~, N{i}] = DistanceReproduction_pooldata(d,'runs',runs,'Distance_N',Distance_N(i),'trialWin',trialWin);
    
    % Find data associated with the desired conflict type
    switch ConflictType
        case 'all'
            ds{i} = ds{i};
            Trials{i} = Trials{i};
        case 'equal'
            ds{i} = ds{i}(abs(ds1{i} - ds2{i}) <= DiffTolerance);
            Trials{i} = Trials{i}(abs(ds1{i} - ds2{i}) <= DiffTolerance);
            dp{i} = dp{i}(abs(ds1{i} - ds2{i}) <= DiffTolerance);
        case 'ts1 > ts2'
            ds{i} = ds{i}(ds1{i} > ds2{i} & abs(ds1{i}-ds2{i}) > DiffTolerance);
            Trials{i} = Trials{i}(ds1{i} > ds2{i} & abs(ds1{i}-ds2{i}) > DiffTolerance);
            dp{i} = dp{i}(ds1{i} > ds2{i} & abs(ds1{i}-ds2{i}) > DiffTolerance);
        case 'ts1 < ts2'
            ds{i} = ds{i}(ds1{i} < ds2{i} & abs(ds1{i}-ds2{i}) > DiffTolerance);
            Trials{i} = Trials{i}(ds1{i} < ds2{i} & abs(ds1{i}-ds2{i}) > DiffTolerance);
            dp{i} = dp{i}(ds1{i} < ds2{i} & abs(ds1{i}-ds2{i}) > DiffTolerance);
    end
    tss = unique(ds{i});
    
    % Estimate Weber fraction
    disp('Estimating Weber fraction')
    if isempty(dp{i})
        weber(i) = 0.1;
    elseif sum(dp{i} > 0)
        fun = @(w)(sum(w.^2*tss.^2) - var(dp{i}(dp{i} > 0 & dp{i} >= MinMaxDp(1) & dp{i} <= tss(end)+MinMaxDp(2))));
        weber(i) = lsqnonlin(fun,0.1);
    else
        weber(i) = NaN;
    end
    
    % Find the offset parameter
    switch ApplyOffsetCorrection.On
        case {'Yes','yes','YES','y','Y',1}
            switch ApplyOffsetCorrection.Source
                case 'SessionAverage'
                    for j = runs
                        bs(j) = d.modelParams{j}.b;
                    end
                    b = mean(bs);
                    
                case 'UserSupplied'
                    b = ApplyOffsetCorrection.b;
            end
        case {'No','no','NO','n','N',0}
            b = 0;
        otherwise
            error(['ApplyOffsetCorrection on flag parameter ' ApplyOffsetCorrection.On ' not recognized!'])
    end
    
    % Sort data by sample time
    errors{i} = [];
    tsIn{i} = [];
    tpIn{i} = [];
    trialsIn{i} = [];
    for ii = 1:length(tss)
        ts_sorted{i}{ii} = ds{i}(ds{i} == tss(ii));
        tp_sorted{i}{ii} = dp{i}(ds{i} == tss(ii)) - b;
        Trials_sorted{i}{ii} = Trials{i}(ds{i} == tss(ii));
        
        mtp(ii,i) = nanmean(tp_sorted{i}{ii}(tp_sorted{i}{ii} >= MinMaxDp(1) & tp_sorted{i}{ii} <= tss(ii)+MinMaxDp(2)));
        
        stdtp(ii,i) = weber(i)*tss(ii);
        
        % Find outliers and recalculate mean and variance
        ins{i}{ii} = find(abs(tp_sorted{i}{ii} - mtp(ii,i)) < outlier*stdtp(ii,i) & tp_sorted{i}{ii} >= MinMaxDp(1) & tp_sorted{i}{ii} <= tss(ii)+MinMaxDp(2));
        
        ts_in{i}{ii} = ts_sorted{i}{ii}(ins{i}{ii});
        tp_in{i}{ii} = tp_sorted{i}{ii}(ins{i}{ii});
        trials_in{i}{ii} = Trials_sorted{i}{ii}(ins{i}{ii});
        tsIn{i} = [tsIn{i}; ts_sorted{i}{ii}(ins{i}{ii})];
        tpIn{i} = [tpIn{i}; tp_sorted{i}{ii}(ins{i}{ii})];
        trialsIn{i} = [trialsIn{i}; Trials_sorted{i}{ii}(ins{i}{ii})];
        
        mtp_in(ii,i) = nanmean(tp_sorted{i}{ii}(ins{i}{ii}));
        
        stdtp_in(ii,i) = nanstd(tp_sorted{i}{ii}(ins{i}{ii}),1);
        
        errors{i} = [errors{i}; tss(ii) - tp_in{i}{ii}];
        
    end
    
    switch LapseRejection
        case 'On'
            [~, ~, ~, ~, like] = prob_tp_take_ts_wm_wp(tpIn{i}-B,tsIn{i},WM,WP,i);
            loglikeLapse = log(LAPSE/(LapseSupport(2)-LapseSupport(1)));
            lapseTrials{i} = log((1-LAPSE)*like) < loglikeLapse;
            
            for ii = 1:length(tss)
                mtp_in(ii,i) = mean(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
                stdtp_in(ii,i) = std(tpIn{i}(tsIn{i} == tss(ii) & ~lapseTrials{i}));
            end
            errors{i} = tsIn{i}(~lapseTrials{i}) - tpIn{i}(~lapseTrials{i});
            
        case 'none'
            lapseTrials{i} = false(size(tsIn{i}));
            
        otherwise
            error(['Lapse rejection option ' LapseRejection ' not recognized!'])
    end
    % Test if production times for longest interval are significantly longer than shortest interval
    tpshort = tp_in{i}{1};
    tplong = tp_in{i}{end};
    
    [h, pval(i)] = ttest2(tplong,tpshort,'tail','right');
    
end         % Distance_N loop

% Calculate RMSE, Bias and Variance
for i = m
    rmse(i) = sqrt(mean((errors{i}+B).^2));
    biassquared(i) = mean((mtp_in(:,i) - B - tss).^2);
    variance(i) = mean(stdtp_in(:,i).^2);
    
    if bootflg
        % Bootstrap bias and variance
        for j = 1:nbootstraps
            es = [];
            for k = 1:length(tss)
                tempP = tpIn{i}(~lapseTrials{i} & tsIn{i} == tss(k));
                es = [es; tempP - tss(k)];
                inds = ceil(length(tempP)*rand(nsamps,1));
                tempm(k,:) = mean(tempP(inds));
                tempstd(k,:) = std(tempP(inds));
%                 inds = ceil(length(tp_in{i}{k})*rand(nsamps,1));
%                 tempm(k,:) = mean(tp_in{i}{k}(inds));
%                 tempstd(k,:) = std(tp_in{i}{k}(inds));
            end
            RMSEs(j,i) = sqrt(mean(es.^2));
            BIASs(j,i) = mean((tempm/G - B - tss).^2);
            VARs(j,i) = mean(tempstd.^2);
        end
    else
        BIASs = NaN(1,length(m));
        VARs = NaN(1,length(m));
    end
end



%% PLOTTING
switch Plot
    case {'Yes','yes','YES','y','Y',1}
        titles = PlotOpts.titles;
        RelativeFigSize = PlotOpts.RelativeFigSize;
        colors = PlotOpts.colors;
        scrsz = get(groot,'ScreenSize');
        maxrmse = max(max(rmse));
        
        figure('Name',[d.sname ' bias vs. sqrt(variance)'],'Position',[scrsz(3) scrsz(4) scrsz(3) scrsz(4)].*RelativeFigSize)
        for i = m
            h(i) = plot(0:0.01:rmse(i),sqrt(rmse(i)^2-(0:0.01:rmse(i)).^2),'Color',colors(i,:));
            hold on
            if bootflg
                plot(sqrt(BIASs(:,i)),sqrt(VARs(:,i)),'.','Color',colors(i,:) + 0.7*(colors(i,:)==0))
            end
            h(i) = plot(sqrt(biassquared(i)),sqrt(variance(i)),'.','Color',colors(i,:),'MarkerSize',20);
        end
        axis([0 maxrmse+50 0 maxrmse+50])
        axis square
        xlabel('bias')
        ylabel('sqrt(variance)')
        legend(h,titles)
        
    case {'No','no','NO','N','n','none','None',0}
    otherwise
        error('Plot option not recognized!');
        
end
        
%% Saving
switch Save
    case {'default','YES','Yes','yes','y','Y',1}
        % Save the results to the default location
        filename = [d.spath '/' d.sname '_' d.projname '_PoolAnalysisResults' datestr(now,'yyyymmdd')];
        disp(['Saving results to ' filename])
        save(filename);
        
    case {'No','no','NO','n','N',0}
        % Don't save the results
    otherwise
        save(Save);     % Save to file specified by the Save string
        
end
