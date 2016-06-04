function DistanceReproduction_fittingOpenMind(SubjectN)
%% DistanceReproduction_fittingOpenMind
%
%   Fits the BLSbiased model to a set of subjects and (optional) saves
%   results.
%
%%

%% Subject map
Subjects = {'JW','MD','SM','SS','VD'};
runmap = {[2,4,5],[2:8],[3:5],[2:6],[2:5]};

%% Variables
% Gobal Variables
trialWin = [100 Inf];
interval_N = 1:2;
MinMax = [0 1000];
outlier = Inf;
FitAll = 1;
FitRuns = 1;
SaveFlg = 1;
dt = 10;

% Fit parameters
fparams.fittype = {'BLSbiasedLapse','aveMeasurements'};     % Specifies which models to fit to the data
fparams.modelUsed = 1;                                      % Specifies which model to use in subsequent model-based analysis
fparams.method = 'quad';                                    % Integration method in model fitting
fparams.init = 'estb';                                      % Initialization values/method for model fitting
fparams.dx = dt;                                            % Step size of integration (in ms)
fparams.trialtypes = [1 2];                                 % Trial types to fit

% Cross validation control
fparams.CrossValidation.Type = 'LNOCV';                              % Cross validation method
fparams.CrossValidation.N = 100;                                       % Left out trials for validation

% Model evidence control
fparams.ModelEvidence.method = 'none';
fparams.ModelEvidence.paramLimits = [0.0001 1;...
                             0.0001 1;...
                             -200 200;...
                             0 1];
fparams.ModelEvidence.integrationMethod = 'quad';
fparams.ModelEvidence.integrationOptions.dx = 0.1;
fparams.ModelEvidence.OpenMind = 1;

% Bias/Variance bootstrap parameters
bootparams.nbootstraps = 100;
bootparams.nsamps = 500;

% Parameters for calculating expected aim times
TAexpectation.method = 'numerical';
TAexpectation.trialtypes = [1 2];
TAexpectation.ts_vec = (550:10:1050)';
TAexpectation.simtrials = 10000;

runs = runmap{SubjectN};

% Run through each subject and fit
disp(['Subject ' Subjects{SubjectN}])

% Load the data
d = load([Subjects{SubjectN} '_DistanceReproduction']);


if SaveFlg
    SaveFileBase = ['/home/swegger/Projects/DistanceReproduction/' Subjects{SubjectN} '/' Subjects{SubjectN} '_BLSbiasedFitResults'];
    SaveParam = [SaveFileBase datestr(now,'yyyymmdd')];
else
    SaveParam = 'No';
end


% Fit all the data
if FitAll
    [mtp, stdtp, bias, variance, rmse, wm, wp, b, pval, weber, tsIn, tpIn, trialsIn, ts_in, tp_in, Trials_sorted, estb, ta, G, Llikelihood, BIASs, VARs, lapse, lapseTrials] = RS2G_psychophysicsPoolAnalysis(d,'runs',runs,'interval_N',interval_N,'Fit',fparams,'Plot','No',...
        'outlier',outlier,'ConflictType','equal','MinMax',MinMax,'Bootstrap',bootparams,...
        'TAexpectation',TAexpectation,'trialWin',trialWin,'Save',SaveParam);
end


