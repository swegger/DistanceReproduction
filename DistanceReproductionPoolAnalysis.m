function [mdp_in, stddp_in, bias, variance, rmse, WM, WP, B, pval, weber, dsIn,...
    dpIn, trialsIn, ds_in, dp_in, Trials_sorted, estb, ta, G, Llikelihood,...
    BIASs, VARs, lapse, lapseTrials, simbias, simv, RMSEs, LmodelEvidence,...
    notFitLlikelihood, notFitLmodelEvidence, SIGP] = DistanceReproductionPoolAnalysis(d,varargin)
%% DistanceReproductionPoolAnalysis
%
%     
%
%%

%% Defaults
Fit_default.fittype = 'none';       % No fit by default
Fit_default.trialtypes = 1;
Fit_default.CrossValidation.Type = 'None';
Fit_default.ModelEvidence.method = 'None';
Fit_default.modelUsed = 1;
bootstrap_default.nbootstraps = NaN;
bootstrap_default.nsamps = NaN;
PlotOpts_default.titles = {'$N = 1$','$N = 2$'};
PlotOpts_default.RelativeFigSize = [1/5 1/2 3/5 1/3];
PlotOpts_default.colors = [0 0 1; 1 0 0; 0.6 0.6 0.6; 0 0 0];
DAexpectation_default.method = 'none';

%% Parse input
Parser = inputParser;

addRequired(Parser,'d')     % Data structure
addParameter(Parser,'runs',NaN)     % Runs to analyze
addParameter(Parser,'Distance_N',1:2)    % Trial types to analyze
addParameter(Parser,'outlier',Inf)      % Number of standard deviations away from the mean to exclude data as outlier
addParameter(Parser,'trialWin',[1 Inf]) % Trials in each run to analyze ([start# end#])
addParameter(Parser,'Fit',Fit_default)  % Fitting options
addParameter(Parser,'Plot','none')      % Data to plot
addParameter(Parser,'PlotOpts',PlotOpts_default)    % Plotting options
addParameter(Parser,'Save','No')           % Saving options
addParameter(Parser,'WeberFractionCheck',0) % Check fit of weber fraction against data
addParameter(Parser,'Bootstrap',bootstrap_default)  % Bootstrap of bias/variance
addParameter(Parser,'MinMaxDp',[-Inf Inf])      % Minimum and maximum values of t_p to keep
addParameter(Parser,'ConflictType','equal')     % For experiments with cue conflict
addParameter(Parser,'DiffTolerance',2/60)       % Tolerance for difference in sample distances before calling it conflict
addParameter(Parser,'DAexpectation',DAexpectation_default)  % For controlling the calculation of the expected value of aim distances under a model
addParameter(Parser,'OutlierRejectionRounds',3)
addParameter(Parser,'viewDistance',310) % Programmed viewing distance; for converting deg to mm

parse(Parser,d,varargin{:})

d = Parser.Results.d;
runs = Parser.Results.runs;
Distance_N = Parser.Results.Distance_N;
outlier = Parser.Results.outlier;
trialWin = Parser.Results.trialWin;
Fit = Parser.Results.Fit;
Plot = Parser.Results.Plot;
PlotOpts = Parser.Results.PlotOpts;
Save = Parser.Results.Save;
WeberFractionCheck = Parser.Results.WeberFractionCheck;
Bootstrap = Parser.Results.Bootstrap;
MinMaxDp = Parser.Results.MinMaxDp;
ConflictType = Parser.Results.ConflictType;
DiffTolerance = Parser.Results.DiffTolerance;
DAexpectation = Parser.Results.DAexpectation;
OutlierRejectionRounds = Parser.Results.OutlierRejectionRounds;
viewDistance = Parser.Results.viewDistance;

% Check to see if run information was provided
if isnan(runs)
    runs = 1:d.runs;            % Defaults to pool data from all runs
end

% Set m to Distance_Ns
m = Distance_N;

% Pull out bootstrap parameters
if any(isnan([Bootstrap.nbootstraps Bootstrap.nsamps]))
    bootflg = 0;
else
    bootflg = 1;
    nbootstraps = Bootstrap.nbootstraps;
    nsamps = Bootstrap.nsamps;
end

% Determine if fitting options have fields for CrossValidation and
% ModelEvdience
if ~isfield(Fit,'CrossValidation')
    Fit.CrossValidation.Type = 'None';
end
if ~isfield(Fit,'ModelEvidence')
    Fit.ModelEvidence.method = 'None';
end
if ~isfield(Fit,'modelUsed')
    Fit.modelUsed = 1;      % Use the first model fit by default
end

%% Analyze data

for i = m
    % Grab the appropriate data
    [ds1{i}, ds2{i}, ds{i}, dp{i}, ~, ~, Trials{i}, correct{i}, ~, ~, ~, N{i}] = DistanceReproduction_pooldata(d,'runs',runs,'Distance_N',i,'trialWin',trialWin);
    
    % Convert from deg to mm
    ds1{i} = viewDistance*tand(ds1{i});
    ds2{i} = viewDistance*tand(ds2{i});
    ds{i} = viewDistance*tand(ds{i});
    dp{i} = viewDistance*tand(dp{i});
    
    % Find data associated with the desired conflict type
    switch ConflictType
        case 'all'
            ds{i} = ds{i};
            Trials{i} = Trials{i};
        case 'equal'
            ds{i} = ds{i}(abs(ds1{i} - ds2{i}) <= DiffTolerance);
            Trials{i} = Trials{i}(abs(ds1{i} - ds2{i}) <= DiffTolerance);
            dp{i} = dp{i}(abs(ds1{i} - ds2{i}) <= DiffTolerance);
        case 'ds1 > ds2'
            ds{i} = ds{i}(ds1{i} > ds2{i} & abs(ds1{i}-ds2{i}) > DiffTolerance);
            Trials{i} = Trials{i}(ds1{i} > ds2{i} & abs(ds1{i}-ds2{i}) > DiffTolerance);
            dp{i} = dp{i}(ds1{i} > ds2{i} & abs(ds1{i}-ds2{i}) > DiffTolerance);
        case 'ds1 < ds2'
            ds{i} = ds{i}(ds1{i} < ds2{i} & abs(ds1{i}-ds2{i}) > DiffTolerance);
            Trials{i} = Trials{i}(ds1{i} < ds2{i} & abs(ds1{i}-ds2{i}) > DiffTolerance);
            dp{i} = dp{i}(ds1{i} < ds2{i} & abs(ds1{i}-ds2{i}) > DiffTolerance);
    end
    dss = unique(ds{i});
    
    % Estimate Weber fraction
    disp('Estimating Weber fraction')
    if isempty(dp{i})
        weber(i) = 0.1;
    elseif sum(dp{i} > 0)
        fun = @(w)(sum(w.^2*dss.^2) - var(dp{i}(dp{i} > 0 & dp{i} >= MinMaxDp(1) & dp{i} <= dss(end)+MinMaxDp(2))));
        weber(i) = lsqnonlin(fun,0.1);
    else
        weber(i) = NaN;
    end
    
    % Sort data by sample distance
    errors{i} = [];
    dsIn{i} = [];
    dpIn{i} = [];
    trialsIn{i} = [];
    for ii = 1:length(dss)
        ds_sorted{i}{ii} = ds{i}(ds{i} == dss(ii));
        dp_sorted{i}{ii} = dp{i}(ds{i} == dss(ii));
        Trials_sorted{i}{ii} = Trials{i}(ds{i} == dss(ii));
        
        
        mdp_in(ii,i) = nanmean(dp_sorted{i}{ii}(dp_sorted{i}{ii} >= MinMaxDp(1) & dp_sorted{i}{ii} <= dss(ii)+MinMaxDp(2)));
        
        stddp_in(ii,i) = nanstd(dp_sorted{i}{ii}(dp_sorted{i}{ii} >= MinMaxDp(1) & dp_sorted{i}{ii} <= dss(ii)+MinMaxDp(2)),1); %weber(i)*dss(ii);
        for jj = 1:OutlierRejectionRounds
            mdp(ii,i) = mdp_in(ii,i);
            
            stddp(ii,i) = stddp_in(ii,i);
            
            % Find outliers and recalculate mean and variance
            ins{i}{ii} = find(abs(dp_sorted{i}{ii} - mdp(ii,i)) < outlier*stddp(ii,i) & dp_sorted{i}{ii} >= MinMaxDp(1) & dp_sorted{i}{ii} <= dss(ii)+MinMaxDp(2));
            
            
            mdp_in(ii,i) = nanmean(dp_sorted{i}{ii}(ins{i}{ii}));
            
            stddp_in(ii,i) = nanstd(dp_sorted{i}{ii}(ins{i}{ii}),1);
            
        end
        ds_in{i}{ii} = ds_sorted{i}{ii}(ins{i}{ii});
        dp_in{i}{ii} = dp_sorted{i}{ii}(ins{i}{ii});
        trials_in{i}{ii} = Trials_sorted{i}{ii}(ins{i}{ii});
        dsIn{i} = [dsIn{i}; ds_sorted{i}{ii}(ins{i}{ii})];
        dpIn{i} = [dpIn{i}; dp_sorted{i}{ii}(ins{i}{ii})];
        trialsIn{i} = [trialsIn{i}; Trials_sorted{i}{ii}(ins{i}{ii})];
        errors{i} = [errors{i}; dss(ii) - dp_in{i}{ii}];
        
        
        % Weber fraction fit check
        if WeberFractionCheck
            if ~isempty(dp{i})
                figure
                [n bins] = hist(dp_sorted{i}{ii},15);
                bar(bins,n/sum(n))
                hold on
                plot(bins,mvnpdf(bins',mean(dp_sorted{i}{ii}),std(dp_sorted{i}{ii}).^2)*(bins(2)-bins(1)),'r--')
                plot(bins,mvnpdf(bins',mean(dp_sorted{i}{ii}),dss(1).^2*weber(i).^2)*(bins(2)-bins(1)),'r')
                title('Go1')
            end
        end
    end
    
    
    % Test if production distances for longest interval are significantly longer than shortest interval
    dpshort = dp_in{i}{1};
    dplong = dp_in{i}{end};
    
    [h, pval(i)] = ttest2(dplong,dpshort,'tail','right');
    
end         % Sets loop

% Fit the data
if ~iscell(Fit.fittype)
    fittype{1} = Fit.fittype;
else
    fittype = Fit.fittype;
end
for fits = 1:length(fittype)
    if ~strcmp(fittype{fits},'BLS_wm_wp_sigp') && ~strcmp(fittype{fits},'aveMeas_wm_wp_sigp')
        SIGP(:,fits) = NaN;
    end
    switch fittype{fits}
        case 'BLSbiased'
            switch Fit.method
                case 'quad'
                    disp('Fitting BLS with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(dpIn{i}) - nanmean(dsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0];
                        end
                    end
                    dt = Fit.dx;
                    [WM(:,fits), WP(:,fits), B(:,fits), Llikelihood(:,fits)] = BLSbiased_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt],'N',num2cell(Fit.trialtypes));
                    LmodelEvidence(:,fits) = NaN;
                case 'quad_batch'
                    disp('Fitting BLS with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(dpIn{i}) - nanmean(dsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    [WM(:,fits), WP(:,fits), B(:,fits), Llikelihood(:,fits)] = BLSbiased_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt batchsize],'N',num2cell(Fit.trialtypes));
                    LmodelEvidence(:,fits) = NaN;
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
            lapse(:,fits) = 0;
            if fits == Fit.modelUsed;
                for i = 1:length(dsIn)
                    lapseTrials{i} = ~ones(size(dsIn{i}));
                end
            end
            
        case 'gBLS'
            switch Fit.method
                case 'quad'
                    disp('Fitting BLS with gain model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estg'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estg(i) = regress(dpIn{i},dsIn{i});
                                init(3) = estg(i);
                            case 'default'
                                init = [0.1 0.06 1];
                        end
                    end
                    dt = Fit.dx;
                    [WM(:,fits), WP(:,fits), G(:,fits), Llikelihood(:,fits)] = gBLS_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt],'N',num2cell(Fit.trialtypes));
                    LmodelEvidence(:,fits) = NaN;
                case 'quad_batch'
                    disp('Fitting BLS with gain model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estg'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estg(i) = regress(dpIn{i},dsIn{i});
                                init(3) = estg(i);
                            case 'default'
                                init = [0.1 0.06 1];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    [WM(:,fits), WP(:,fits), G(:,fits), Llikelihood(:,fits)] = gBLS_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt batchsize],'N',num2cell(Fit.trialtypes));
                    LmodelEvidence = NaN;
                otherwise
                    error('Fitting method not recognized for gBLS fitter!')
            end
            B(:,fits) = 0;
            lapse(:,fits) = 0;
            if fits == Fit.modelUsed;
                for i = 1:length(dsIn)
                    lapseTrials{i} = ~ones(size(dsIn{i}));
                end
            end
            
        case 'BLSbiasedLapse'
            switch Fit.method
                case 'quad'
                    % Fit the BLS model with bias and lapses to the data
                    disp('Fitting BLS with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN 0.05];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(dpIn{i}) - nanmean(dsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0 0.05];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxDp(1) max(dss)+MinMaxDp(2)];
                    [WM(:,fits), WP(:,fits), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = BLSbiasedLapse_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = BLSbiasedLapse_Validator(dsIn(~fitted),dpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(dss) max(dss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed;
                        for i = m
                            [~, ~, loglike, ~, like] = prob_dp_take_ds_wm_wp(dpIn{i}-mean(B(:,fits)),dsIn{i},mean(WM(:,fits)),mean(WP(:,fits)),i);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mdp and stddp
                    for i = m
                        for ii = 1:length(dss)
                            mdp_in(ii,i) = mean(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                            stddp_in(ii,i) = std(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = dsIn{i}(~lapseTrials{i}) - dpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    disp('Fitting BLS with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(dpIn{i}) - nanmean(dsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    LapseSupport = [MinMaxDp(1) max(dss)+MinMaxDp(2)];
                    [WM(:,fits), WP(:,fits), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = BLSbiasedLapse_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt batchsize],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = BLSbiasedLapse_Validator(dsIn(~fitted),dpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(dss) max(dss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed;
                        for i = m
                            [~, ~, loglike, ~, like] = prob_dp_take_ds_wm_wp(dpIn{i}-mean(B(:,fits)),dsIn{i},mean(WM(:,fits)),mean(WP(:,fits)),i);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mdp and stddp
                    for i = m
                        for ii = 1:length(dss)
                            mdp_in(ii,i) = mean(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                            stddp_in(ii,i) = std(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = dsIn{i}(~lapseTrials{i}) - dpIn{i}(~lapseTrials{i});
                    end
                    
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
        
        case 'BLS_wm_wp_sigp'
            switch Fit.method
                case 'quad'
                    % Fit the BLS model with bias and lapses to the data
                    disp('Fitting BLS with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN 0.05 0];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(dpIn{i}) - nanmean(dsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0 0.05 0];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxDp(1) max(dss)+MinMaxDp(2)];
                    [WM(:,fits), WP(:,fits), B(:,fits), lapse(:,fits), SIGP(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = BLS_wm_wp_sigp_b_lapse_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = BLS_wm_wp_sigp_b_lapse_Validator(dsIn(~fitted),dpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits)),mean(B(:,fits)),mean(lapse(:,fits)),mean(SIG(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(dss) max(dss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed;
                        for i = m
                            [~, ~, loglike, ~, like] = prob_dp_take_ds_wm_wp_sigp(dpIn{i}-mean(B(:,fits)),dsIn{i},mean(WM(:,fits)),mean(WP(:,fits)),mean(SIGP(:,fits)),i);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mdp and stddp
                    for i = m
                        for ii = 1:length(dss)
                            mdp_in(ii,i) = mean(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                            stddp_in(ii,i) = std(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = dsIn{i}(~lapseTrials{i}) - dpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    disp('Fitting BLS with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN 0];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(dpIn{i}) - nanmean(dsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0 0];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    LapseSupport = [MinMaxDp(1) max(dss)+MinMaxDp(2)];
                    [WM(:,fits), WP(:,fits), B(:,fits), lapse(:,fits). SIGP(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = BLS_wm_wp_sigp_b_lapse_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt batchsize],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = BLS_wm_wp_sigp_b_lapse_Validator(dsIn(~fitted),dpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits)),mean(B(:,fits)),mean(lapse(:,fits)),mean(SIG(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(dss) max(dss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed;
                        for i = m
                            [~, ~, loglike, ~, like] = prob_dp_take_ds_wm_wp_sigp(dpIn{i}-mean(B(:,fits)),dsIn{i},mean(WM(:,fits)),mean(WP(:,fits)),mean(SIGP(:,fits)),i);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mdp and stddp
                    for i = m
                        for ii = 1:length(dss)
                            mdp_in(ii,i) = mean(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                            stddp_in(ii,i) = std(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = dsIn{i}(~lapseTrials{i}) - dpIn{i}(~lapseTrials{i});
                    end
                    
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
        
        case 'MAPbiasedLapse'
            switch Fit.method
                case 'quad'
                    % Fit the BLS model with bias and lapses to the data
                    disp('Fitting MAP with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN 0.05];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(dpIn{i}) - nanmean(dsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0 0.05];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxDp(1) max(dss)+MinMaxDp(2)];
                    [WM(:,fits), WP(:,fits), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = MAPbiasedLapse_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = MAPbiasedLapse_Validator(dsIn(~fitted),dpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(dss) max(dss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed;
                        for i = m
                            [~, ~, loglike, ~, like] = prob_dp_take_ds_wm_wp(dpIn{i}-mean(B(:,fits)),dsIn{i},mean(WM(:,fits)),mean(WP(:,fits)),i);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mdp and stddp
                    for i = m
                        for ii = 1:length(dss)
                            mdp_in(ii,i) = mean(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                            stddp_in(ii,i) = std(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = dsIn{i}(~lapseTrials{i}) - dpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    disp('Fitting MAP with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(dpIn{i}) - nanmean(dsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    LapseSupport = [MinMaxDp(1) max(dss)+MinMaxDp(2)];
                    [WM(:,fits), WP(:,fits), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = MAPbiasedLapse_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt batchsize],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = MAPbiasedLapse_Validator(dsIn(~fitted),dpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(dss) max(dss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed;
                        for i = m
                            [~, ~, loglike, ~, like] = prob_dp_take_ds_wm_wp(dpIn{i}-mean(B(:,fits)),dsIn{i},mean(WM(:,fits)),mean(WP(:,fits)),i);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mdp and stddp
                    for i = m
                        for ii = 1:length(dss)
                            mdp_in(ii,i) = mean(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                            stddp_in(ii,i) = std(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = dsIn{i}(~lapseTrials{i}) - dpIn{i}(~lapseTrials{i});
                    end
                    
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
            
        case 'ObsActBiasedLapse'
            switch Fit.method
                case 'quad'
                    % Fit the BLS model with bias and lapses to the data
                    disp('Fitting BLS with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN 0.05];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(dpIn{i}) - nanmean(dsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0 0.05];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxDp(1) max(dss)+MinMaxDp(2)];
                    [WM(:,fits), WP(:,fits), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = ObsActBiasedLapse_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = ObsActBiasedLapse_Validator(dsIn(~fitted),dpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(dss) max(dss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed;
                        for i = m
                            [~, ~, loglike, ~, like] = prob_dp_take_ds_wm_wp(dpIn{i}-mean(B(:,fits)),dsIn{i},mean(WM(:,fits)),mean(WP(:,fits)),i);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mdp and stddp
                    for i = m
                        for ii = 1:length(dss)
                            mdp_in(ii,i) = mean(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                            stddp_in(ii,i) = std(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = dsIn{i}(~lapseTrials{i}) - dpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    disp('Fitting BLS with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(dpIn{i}) - nanmean(dsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    LapseSupport = [MinMaxDp(1) max(dss)+MinMaxDp(2)];
                    [WM(:,fits), WP(:,fits), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = ObsActBiasedLapse_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt batchsize],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = ObsActBiasedLapse_Validator(dsIn(~fitted),dpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(dss) max(dss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed;
                        for i = m
                            [~, ~, loglike, ~, like] = prob_dp_take_ds_wm_wp(dpIn{i}-mean(B(:,fits)),dsIn{i},mean(WM(:,fits)),mean(WP(:,fits)),i);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mdp and stddp
                    for i = m
                        for ii = 1:length(dss)
                            mdp_in(ii,i) = mean(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                            stddp_in(ii,i) = std(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = dsIn{i}(~lapseTrials{i}) - dpIn{i}(~lapseTrials{i});
                    end
                    
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
            
        case 'aveMeasurements'
            switch Fit.method
                case 'quad'
                    % Fit the aveMeas model with bias and lapses to the data
                    disp('Fitting measurement averaging model with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN 0.05];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(dpIn{i}) - nanmean(dsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0 0.05];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxDp(1) max(dss)+MinMaxDp(2)];
                    [WM(:,fits), WP(:,fits), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = aveMeasbiasedLapse_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = aveMeasbiasedLapse_Validator(dsIn(~fitted),dpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(dss) max(dss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed;
                        for i = m
                            estimator.type = 'weightedMean';
                            estiamtor.weights = ones(1,Fit.trialtypes(i))/Fit.trialtypes(i);
                            [~, ~, loglike, ~, like] = prob_dp_take_ds_wm_wp(dpIn{i}-mean(B(:,fits)),dsIn{i},mean(WM(:,fids)),mean(WP(:,fids)),i,'estimator',estimator);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mdp and stddp
                    for i = m
                        for ii = 1:length(dss)
                            mdp_in(ii,i) = mean(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                            stddp_in(ii,i) = std(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = dsIn{i}(~lapseTrials{i}) - dpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    disp('Fitting measurement averaging model with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(dpIn{i}) - nanmean(dsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    LapseSupport = [MinMaxDp(1) max(dss)+MinMaxDp(2)];
                    [WM(:,fits), WP(:,fits), B(:,fits), lapse(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = aveMeasbiasedLapse_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt batchsize],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = aveMeasbiasedLapse_Validator(dsIn(~fitted),dpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(dss) max(dss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed;
                        for i = m
                            estimator.type = 'weightedMean';
                            estiamtor.weights = ones(1,Fit.trialtypes(i))/Fit.trialtypes(i);
                            [~, ~, loglike, ~, like] = prob_dp_take_ds_wm_wp(dpIn{i}-mean(B(:,fits)),dsIn{i},mean(WM(:,fits)),mean(WP(:,fits)),i,'estimator',estimator);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mdp and stddp
                    for i = m
                        for ii = 1:length(dss)
                            mdp_in(ii,i) = mean(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                            stddp_in(ii,i) = std(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = dsIn{i}(~lapseTrials{i}) - dpIn{i}(~lapseTrials{i});
                    end
                    
                otherwise
                    error('Fitting method not recognized for BLSbiased fitter!')
            end
            G(:,fits) = 1;
        
        case 'aveMeas_wm_wp_sigp'
            switch Fit.method
                case 'quad'
                    % Fit the aveMeas model with bias and lapses to the data
                    disp('Fitting measurement averaging model with bias and lapses model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN 0.05 0];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(dpIn{i}) - nanmean(dsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0 0.05 0];
                        end
                    end
                    dt = Fit.dx;
                    LapseSupport = [MinMaxDp(1) max(dss)+MinMaxDp(2)];
                    [WM(:,fits), WP(:,fits), B(:,fits), lapse(:,fits), SIGP(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = aveMeas_wm_wp_sigp_b_lapse_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = aveMeas_wm_wp_sigp_b_lapse_Validator(dsIn(~fitted),dpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(dss) max(dss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed;
                        for i = m
                            estimator.type = 'weightedMean';
                            estiamtor.weights = ones(1,Fit.trialtypes(i))/Fit.trialtypes(i);
                            [~, ~, loglike, ~, like] = prob_dp_take_ds_wm_wp_sigp(dpIn{i}-mean(B(:,fits)),dsIn{i},mean(WM(:,fids)),mean(WP(:,fids)),mean(SIGP(:,fids)),i,'estimator',estimator);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mdp and stddp
                    for i = m
                        for ii = 1:length(dss)
                            mdp_in(ii,i) = mean(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                            stddp_in(ii,i) = std(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = dsIn{i}(~lapseTrials{i}) - dpIn{i}(~lapseTrials{i});
                    end
                    
                case 'quad_batch'
                    disp('Fitting measurement averaging model with bias model to the data using Simpsons quadrature...')
                    init = Fit.init;
                    if isstr(init)
                        switch init
                            case 'estb'
                                init = [0.1 0.06 NaN 0];      % Default initial search but with estimate of baseline
                                estb(i) = nanmean(dpIn{i}) - nanmean(dsIn{i});
                                init(3) = estb(i);
                            case 'default'
                                init = [0.1 0.06 0 0];
                        end
                    end
                    dt = Fit.dx;
                    if isfield(varargin{fitnum+1},'batchsize')
                        batchsize = varargin{fitnum+1}.batchsize;
                    else
                        batchsize = 1000000;
                    end
                    LapseSupport = [MinMaxDp(1) max(dss)+MinMaxDp(2)];
                    [WM(:,fits), WP(:,fits), B(:,fits), lapse(:,fits), SIGP(:,fits), Llikelihood(:,fits), LmodelEvidence(:,fits)] = aveMeas_wm_wp_sigp_b_lapse_fitter(dsIn(Fit.trialtypes),dpIn(Fit.trialtypes),'InitCond',init,'FitType',Fit.method,[min(dss) max(dss) dt batchsize],'N',num2cell(Fit.trialtypes),'LapseSupport',LapseSupport,...
                        'CrossValidation',Fit.CrossValidation,'ModelEvidence',Fit.ModelEvidence);
                    
                    % Find likelihood of model on left out condition, if it exists
                    for i = m
                        fitted(i) = any(i == Fit.trialtypes);
                    end
                    if any(~fitted)
                        [notFitLlikelihood(:,fits), notFitLmodelEvidence(:,fits)] = aveMeas_wm_wp_sigp_b_lapse_Validator(dsIn(~fitted),dpIn(~fitted),mean(WM(:,fits)),mean(WP(:,fits)),mean(B(:,fits)),mean(lapse(:,fits)),...
                            'N',num2cell(m(~fitted)),'LapseSupport',LapseSupport,'ModelEvidence',Fit.ModelEvidence,'FitType',Fit.method,[min(dss) max(dss) dt]); 
                    else
                        notFitLlikelihood(:,fits) = NaN;
                        notFitLmodelEvidence(:,fits) = NaN;
                    end
                    
                    % Identify lapse trials
                    if fits == Fit.modelUsed;
                        for i = m
                            estimator.type = 'weightedMean';
                            estiamtor.weights = ones(1,Fit.trialtypes(i))/Fit.trialtypes(i);
                            [~, ~, loglike, ~, like] = prob_dp_take_ds_wm_wp_sigp(dpIn{i}-mean(B(:,fits)),dsIn{i},mean(WM(:,fits)),mean(WP(:,fits)),mean(SIGP(:,fits)),i,'estimator',estimator);
                            loglikeLapse = log(mean(lapse(:,fits))/(LapseSupport(2)-LapseSupport(1)));
                            lapseTrials{i} = log((1-mean(lapse(:,fits)))*like) < loglikeLapse;
                        end
                    end
                    
                    % Recalculate mdp and stddp
                    for i = m
                        for ii = 1:length(dss)
                            mdp_in(ii,i) = mean(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                            stddp_in(ii,i) = std(dpIn{i}(dsIn{i} == dss(ii) & ~lapseTrials{i}));
                        end
                        errors{i} = dsIn{i}(~lapseTrials{i}) - dpIn{i}(~lapseTrials{i});
                    end
            end
            
        case 'none'
            WM(:,fits) = NaN(m(end),1);
            WP(:,fits) = NaN(m(end),1);
            B(:,fits) = 0;
            G(:,fits) = 1;
            Llikelihood(:,fits) = NaN(m(end),1);
            estb(:,fits) = NaN(m(end),1);
            lapse(:,fits) = 0;
            for i = 1:length(dsIn)
                lapseTrials{i} = ~ones(size(dsIn{i}));
            end
            LmodelEvidence(:,fits) = NaN;
            notFitLlikelihood(:,fits) = NaN;
            notFitLmodelEvidence(:,fits) = NaN;
            
        otherwise
            warning('Fit type argument not recognized, no fit performed.')
            WM(:,fits) = NaN(m(end),1);
            WP(:,fits) = NaN(m(end),1);
            B(:,fits) = 0;
            G(:,fits) = 1;
            Llikelihood(:,fits) = NaN(m(end),1);
            estb(:,fits) = NaN(m(end),1);
            lapse(:,fits) = 0;
            for i = 1:length(dsIn)
                lapseTrials{i} = ~ones(size(dsIn{i}));
            end
            LmodelEvidence(:,fits) = NaN;
            notFitLlikelihood(:,fits) = NaN;
            notFitLmodelEvidence(:,fits) = NaN;
            
    end
end

% Find the parameters to use for the remaining analyses
WMhat = mean(WM(:,Fit.modelUsed));
WPhat = mean(WP(:,Fit.modelUsed));
Bhat = mean(B(:,Fit.modelUsed));
Ghat = mean(G(:,Fit.modelUsed));
SIGPhat = mean(SIGP(:,Fit.modelUsed));

% Calculate RMSE, Bias and Variance
for i = m
    rmse(i) = sqrt(mean((errors{i}+Bhat).^2));
    bias(i) = mean((mdp_in(:,i) - Bhat - dss).^2);
    variance(i) = mean(stddp_in(:,i).^2);
    
    if bootflg
        % Bootstrap bias and variance
        for j = 1:nbootstraps
            es = [];
            for k = 1:length(dss)
                tempP = dpIn{i}(~lapseTrials{i} & dsIn{i} == dss(k));
                inds = ceil(length(tempP)*rand(nsamps,1));
                es = [es; tempP(inds)/Ghat - Bhat - dss(k)];
                tempm(k,:) = mean(tempP(inds));
                tempstd(k,:) = std(tempP(inds));
            end
            RMSEs(j,i) = sqrt(mean(es.^2));
            BIASs(j,i) = mean((tempm/Ghat - Bhat - dss).^2);
            VARs(j,i) = mean(tempstd.^2);
        end
    else
        RMSEs = NaN(1,length(m));
        BIASs = NaN(1,length(m));
        VARs = NaN(1,length(m));
    end
end

% Generate expected value of aim distances based on model fits
switch DAexpectation.method
    case 'none'
        ds_vec = dss;
        ta = nan(length(dss),length(m));
        simbias = nan(1,length(Distance_N));
        simv = nan(1,length(Distance_N));
    case 'numerical'
        if ~isfield(DAexpectation,'ds_vec') && ~isfield(DAexpectation,'dt')
            ds_vec = dss;
            ds_vec = ds_vec(:);
        elseif ~isfield(DAexpectation,'ds_vec')
            ds_vec = dss(1):DAexpectation.dt:dss(end);
            ds_vec = ds_vec(:);
        else
            ds_vec = DAexpectation.ds_vec;
            ds_vec = ds_vec(:);
        end
        if isfield(DAexpectation,'simtrials')
            simtrials = DAexpectation.simtrials;
        else
            simtrials = 1000;
        end
        for i = m
            switch fittype{Fit.modelUsed}
                case {'BLSbiased','BLS','BLSbiasedLapse'}
                    method_opts.dx = 0.01;
                    if i >= 4
                        ta(:,i) = NaN(size(ds_vec));
                    else
                        ta(:,i) = ta_expectation3(ds_vec,WMhat,i,dt,'Type','BLS','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(dss) max(dss)]);
                    end
                    [~, ~, simbias(i), simv(i)] = ta_expectation3(dss',WMhat,i,dt,'method','numerical','trials',simtrials,'wp',WPhat,'Support',[min(dss) max(dss)]);
                    
                case 'gBLS'
                    method_opts.dx = 0.01;
                    options.g = Ghat;
                    if i >= 4
                        ta(:,i) = NaN(size(ds_vec));
                    else
                        ta(:,i) = ta_expectation3(ds_vec,WMhat,i,dt,'Type','gBLS','method_options',method_opts,'options',options,'method',DAexpectation.method,'Support',[min(dss) max(dss)]);
                    end
                    
                    if i == 1
                        [~, ~, simbias(1), simv(1)] = ta_expectation3(dss',WMhat,1,dt,'method','numerical','trials',1000,'wp',WPhat,'Support',[min(dss) max(dss)]);
                    elseif i == 2
                        [~, ~, simbias(2), simv(2)] = ta_expectation3(dss',WMhat,2,dt,'method','numerical','trials',1000,'wp',WPhat,'Support',[min(dss) max(dss)]);
                    elseif i > 2
                        [~, ~, simbias(i), simv(i)] = ta_expectation3(dss',WMhat,i,dt,'method','numerical','trials',1000,'wp',WPhat,'Support',[min(dss) max(dss)]);
                    end
                case {'ObsAct','ObsActBiased','ObsActBiasedLapse'}
                    method_opts.dx = 0.01;
                    if i >= 4
                        [ta(:,i), ta_std(:,i)] = NaN(size(ds_vec));
                    else
                        [ta(:,i), ta_std(:,i)] = ta_expectation3(ds_vec,WMhat,i,dt,'Type','ObsAct','method_options',method_opts,'method','numerical','trials',simtrials,'wp',WPhat,'Support',[min(dss) max(dss)]);
                    end
                    [~, ~, simbias(i), simv(i)] = ta_expectation3(dss',WMhat,i,dt,'Type','ObsAct','method','numerical','trials',simtrials,'wp',WPhat,'Support',[min(dss) max(dss)]);
                
                case 'aveMeasurements'
                    method_opts.dx = 0.01;
                    ta(:,i) = ta_expectation3(ds_vec,WMhat,i,dt,'Type','aveMeasurements','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(dss) max(dss)]);
                    [~, ~, simbias(i), simv(i), simrmse(i)] = ta_expectation3(dss',WMhat,i,dt,'Type','aveMeasurements','method_options',method_opts,'method','numerical','trials',simtrials,'wp',WPhat,'Support',[min(dss) max(dss)]);
                    
                case 'MAP'
                    method_opts.dx = 0.01;
                    ta(:,i) = ta_expectation3(ds_vec,WMhat,i,dt,'Type','MAP','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'Support',[min(dss) max(dss)]);
                    [~, ~, simbias(i), simv(i), simrmse(i)] = ta_expectation3(dss',WMhat,i,dt,'Type','MAP','method_options',method_opts,'method','numerical','trials',simtrials,'wp',WPhat,'Support',[min(dss) max(dss)]);
                   
                case 'BLS_wm_wp_sigp'
                    method_opts.dx = 0.01;
                    if i >= 4
                        ta(:,i) = NaN(size(ds_vec));
                    else
                        ta(:,i) = ta_expectation3(ds_vec,WMhat,i,dt,'Type','BLS_wm_wp_sigp','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'sigp',SIGPhat,'Support',[min(dss) max(dss)]);
                    end
                    [~, ~, simbias(i), simv(i)] = ta_expectation3(dss',WMhat,i,dt,'Type','BLS_wm_wp_sigp','method','numerical','trials',simtrials,'wp',WPhat,'sigp',SIGPhat,'Support',[min(dss) max(dss)]);
                    
                    
                case 'aveMeas_wm_wp_sigp'
                    method_opts.dx = 0.01;
                    ta(:,i) = ta_expectation3(ds_vec,WMhat,i,dt,'Type','aveMeasurements','method_options',method_opts,'method','numerical','trials',simtrials,'wp',0,'sigp',SIGPhat,'Support',[min(dss) max(dss)]);
                    [~, ~, simbias(i), simv(i), simrmse(i)] = ta_expectation3(dss',WMhat,i,dt,'Type','aveMeasurements','method_options',method_opts,'method','numerical','trials',simtrials,'wp',WPhat,'sigp',SIGPhat,'Support',[min(dss) max(dss)]);
                    
                otherwise
                    ta = NaN;
                    simbias(i) = NaN;
                    simv(i) = NaN;
                    simrmse(i) = NaN;
            end
        end
        
    case 'analytical'
        if ~isfield(DAexpectation,'ds_vec') && ~isfield(DAexpectation,'dt')
            ds_vec = dss;
            ds_vec = ds_vec(:);
        elseif ~isfield(DAexpectation,'ds_vec')
            ds_vec = dss(1):DAexpectation.dt:dss(end);
            ds_vec = ds_vec(:);
        else
            ds_vec = DAexpectation.ds_vec;
            ds_vec = ds_vec(:);
        end
        if isfield(DAexpectation,'simtrials')
            simtrials = DAexpectation.simtrials;
        else
            simtrials = 1000;
        end
        for i = m
            switch fittype{Fit.modelUsed}
                case {'BLSbiased','BLS'}
                    method_opts.dx = 0.01;
                    if i >= 4
                        ta(:,i) = NaN(size(ds_vec));
                    else
                        ta(:,i) = ta_expectation3(ds_vec,WMhat,i,dt,'Type','BLS','method_options',method_opts,'method','analytical','Support',[min(dss) max(dss)]);
                    end
                    [~, ~, simbias(i), simv(i)] = ta_expectation3(dss,WMhat,i,dt,'method','numerical','trials',simtrials,'wp',WPhat,'Support',[min(dss) max(dss)]);
                    
                case 'gBLS'
                    method_opts.dx = 0.01;
                    options.g = Ghat;
                    if i >= 4
                        ta(:,i) = NaN(size(ds_vec));
                    else
                        ta(:,i) = ta_expectation3(ds_vec,WMhat,i,dt,'Type','gBLS','method_options',method_opts,'options',options,'method',DAexpectation.method,'Support',[min(dss) max(dss)]);
                    end
                    
                    if i == 1
                        [~, ~, simbias(1), simv(1)] = ta_expectation3(dss',WMhat,1,dt,'method','numerical','trials',1000,'wp',WPhat,'Support',[min(dss) max(dss)]);
                    elseif i == 2
                        [~, ~, simbias(2), simv(2)] = ta_expectation3(dss',WMhat,2,dt,'method','numerical','trials',1000,'wp',WPhat,'Support',[min(dss) max(dss)]);
                    elseif i > 2
                        [~, ~, simbias(i), simv(i)] = ta_expectation3(dss',WMhat,i,dt,'method','numerical','trials',1000,'wp',WPhat,'Support',[min(dss) max(dss)]);
                    end
                case 'ObsAct'
                    error('Not yet supported')
                    ta(:,i) = ta_expectation3(ds_vec,WMhat,1,'Type','ObsAct','wp',wp(i),'Support',[min(dss) max(dss)]);
                    
                    if i == 2
                        ta2 = ta_expectation(ds_vec,wm(1),2,'Type','ObsAct','wp',wp(1),'Support',[min(dss) max(dss)]);
                    elseif i > 2
                        error('Not yet supported')
                    end
                    
                case 'aveMeasurements'
                    method_opts.dx = 0.01;
                    if i >= 4
                        ta(:,i) = NaN(size(ds_vec));
                    else
                        ta(:,i) = ta_expectation3(ds_vec,WMhat,i,dt,'Type','aveMeasurements','method_options',method_opts,'method','analytical','Support',[min(dss) max(dss)]);
                    end
                    [~, ~, simbias(i), simv(i), simrmse(i)] = ta_expectation3(dss,WMhat,i,dt,'Type','aveMeasurements','method','numerical','trials',simtrials,'wp',WPhat,'Support',[min(dss) max(dss)]);
                    
                otherwise
                    ta = NaN;
                    simbias(i) = NaN;
                    simv(i) = NaN;
                    
            end
        end
        
    otherwise
        error(['Model expectation method ' DAexpectation.method ' not recognized!']);
end

%% PLOTTING
switch Plot
    case {'All','all'}
        titles = PlotOpts.titles;
        RelativeFigSize = PlotOpts.RelativeFigSize;
        colors = PlotOpts.colors;
        
        % Dependence on sample distance
        scrsz = get(groot,'ScreenSize');
        figure('Name',[d.sname ' dependence on sample distance'],'Position',[scrsz(3) scrsz(4) scrsz(3) scrsz(4)].*RelativeFigSize)
        maxrmse = max(max(rmse));
        plotind = 1;
        allds = [];
        alldp = [];
        for i = m
            for ii = 1:length(ds_in{i})
                allds = [allds; ds_in{i}{ii}];
                alldp = [alldp; dp_in{i}{ii}];
            end
        end
        ax = [min(alldp)-0.2*min(alldp) max(alldp)+0.2*max(alldp) min(alldp)-0.2*min(alldp) max(alldp)+0.2*max(alldp)];
        for i = m
            subplot(1,length(m),plotind)
            plotind = plotind+1;
            plot(dsIn{i},dpIn{i},'o','Color',colors(i,:)+(1 - colors(i,:))/1.5)
            hold all
            plot(dsIn{i}(lapseTrials{i}),dpIn{i}(lapseTrials{i}),'.','Color',[0 0 0])
            text(dss(end),dss(1)-100,['p = ' num2str(pval(i))]);
            axis(ax)
            xlabel('$d_s$')
            ylabel('$d_p$')
        end
        
        plotind = 1;
        for i = m
            subplot(1,length(m),plotind)
            plotind = plotind+1;
            errorbar(dss,mdp_in(:,i),stddp_in(:,i),'.','Color',colors(i,:),'MarkerSize',20,'LineWidth',2)
            hold all
            if ~strcmp('none',fittype{Fit.modelUsed}) && any(i == Distance_N)
                plot(ds_vec,ta(:,i)+Bhat,'Color',colors(i,:))
            end
            axis(ax)
            plotUnity;
            mymakeaxis(gca,'xytitle',titles{i},'interpreter','latex')
        end
        
        figure('Name',[d.sname ' bias vs. sqrt(variance)'])
        for i = m
            h(i) = plot(0:0.01:rmse(i),sqrt(rmse(i)^2-(0:0.01:rmse(i)).^2),'Color',colors(i,:));
            hold on
            if bootflg
                plot(sqrt(BIASs(:,i)),sqrt(VARs(:,i)),'.','Color',colors(i,:) + 0.7*(colors(i,:)==0))
            end
            if ~strcmp('none',fittype{Fit.modelUsed})
                h(i) = plot(sqrt(simbias(i)),sqrt(simv(i)),'o','Color',colors(i,:));
            end
            h(i) = plot(sqrt(bias(i)),sqrt(variance(i)),'.','Color',colors(i,:),'MarkerSize',20);
        end
        axis([0 maxrmse+0.5 0 maxrmse+0.5])
        axis square
        xlabel('BIAS')
        ylabel('VAR$^{\frac{1}{2}}$')
        legend(h,titles,'interpreter','latex')
        mymakeaxis(gca,'interpreter','latex')
        
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
