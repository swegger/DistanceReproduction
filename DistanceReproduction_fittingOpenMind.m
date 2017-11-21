function DistanceReproduction_fittingOpenMind(SubjectN,varargin)
%% DistanceReproduction_fittingOpenMind
%
%   Fits the BLSbiased model to a set of subjects and (optional) saves
%   results.
%
%%

%% Defaults

%% Gobal Variables

% Subject map
Subjects = {'JW','MD','SM','SS','SWE','TT','VD'};
runmap = {[2,4,5:8],[2:10],[3:6],[2:9,11],[2:10],[2:10],[2:5]};

trialWin = [100 Inf];
Distance_N = 1:2;
MinMax = [2 30];
outlier = Inf;
FitAll = 1;
SaveFlg = 1;
dt = 5;
viewDistance = 310;
fixPos = 13;

% Fit parameters
fparams_default.fittype = {'BLS_wm_wp_sigp','aveMeas_wm_wp_sigp'};     % Specifies which models to fit to the data
fparams_default.modelUsed = 1;                                      % Specifies which model to use in subsequent model-based analysis
fparams_default.method = 'quad';                                    % Integration method in model fitting
fparams_default.init = 'estb';                                      % Initialization values/method for model fitting
fparams_default.dx = dt;                                            % Step size of integration (in ms)
fparams.trialtypes = [1 2];                                 % Trial types to fit

% Cross validation control
fparams_default.CrossValidation.Type = 'LNOCV';                              % Cross validation method
fparams_default.CrossValidation.N = 100;                                       % Left out trials for validation

% Model evidence control
fparams_default.ModelEvidence.method = 'none';
fparams_default.ModelEvidence.paramLimits = [0.0001 1;...
                             0.0001 1;...
                             -200 200;...
                             0 1;...
                             0 100];
fparams_default.ModelEvidence.integrationMethod = 'quad';
fparams_default.ModelEvidence.integrationOptions.dx = 0.5;
fparams_default.ModelEvidence.OpenMind = 1;

% Bias/Variance bootstrap parameters
bootparams.nbootstraps = 100;
bootparams.nsamps = 500;

% Parameters for calculating expected aim times
DAexpectation.method = 'numerical';
DAexpectation.trialtypes = [1 2];
DAexpectation.ds_vec = viewDistance*(tand(fixPos) - tand(fixPos - (13:0.1:19)'));
DAexpectation.simtrials = 10000;

runs = runmap{SubjectN};

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'SubjectN')
addParameter(Parser,'SaveOpts',SaveOpts_default)
addParameter(Parser,'fparams',fparams_default)
addParameter(Parser,'TestSave',false)

parse(Parser,SubjectN,varargin{:})

SubjectN = Parser.Results.SubjectN;
SaveOpts = Parser.Results.SaveOpts;
fparams = Parser.Results.fparams;
TestSave = Parser.Results.TestSave;

fparams = asign_fparams(fparams);
if SaveOpts.SaveFlg && ~isfield(SaveOpts,'SaveFile')
    SaveOpts.SaveFile = ['/om/user/swegger/Projects/DistanceReproduction/'...
        Subjects{SubjectN} '/' Subjects{SubjectN} '_' ...
        fparams.fittype{fparams.modelUsed} '_ObsAct' num2str(fparams.ObsAct)...
        '_' datestr(now,'yyyymmdd')];
end


%%
% Load the data
d = load([Subjects{SubjectN} '_DistanceReproduction']);


% if SaveFlg
%     SaveFileBase = ['/home/swegger/Projects/DistanceReproduction/' Subjects{SubjectN} '/' Subjects{SubjectN} '_BLSbiasedFitResults'];
%     SaveParam = [SaveFileBase datestr(now,'yyyymmdd')];
% else
%     SaveParam = 'No';
% end
if SaveOpts.SaveFlg
    SaveParam = SaveOpts.SaveFile;
else
    SaveParam = 'No';
end

% Fit all the data
if FitAll
    [mtp, stdtp, bias, variance, rmse, wm, wp, b, pval, weber, tsIn, tpIn,...
        trialsIn, ts_in, tp_in, Trials_sorted, estb, ta, G, Llikelihood,...
        BIASs, VARs, lapse, lapseTrials] = DistanceReproductionPoolAnalysis(d,...
        'runs',runs,'Distance_N',Distance_N,'Fit',fparams,'Plot','No',...
        'outlier',outlier,'ConflictType','equal','MinMax',MinMax,'Bootstrap',bootparams,...
        'DAexpectation',DAexpectation,'trialWin',trialWin,'Save',SaveParam);
end


%% Functions

    function fparams = asign_fparams(fparams)
    if ~isfield(fparams,'fittype')
        fparams.fittype = {'SubOptMemBias','BLSbiasedLapse'};
    end
    if ~isfield(fparams,'modelUsed')
        fparams.modelUsed = 1;                                      % Specifies which model to use in subsequent model-based analysis
    end
    if ~isfield(fparams,'method')
        fparams.method = 'quad';                                    % Integration method in model fitting
    end
    if ~isfield(fparams,'init')
        fparams.init = 'estb';                                      % Initialization values/method for model fitting
    end
    if ~isfield(fparams,'dx')
        fparams.dx = 5;                                            % Step size of integration (in ms)
    end
    if ~isfield(fparams,'trialtypes')
        fparams.trialtypes = [1 2];                                 % Trial types to fit
    end
    if ~isfield(fparams,'ObsAct')
        fparams.ObsAct = 0;
    end

    % Cross validation control
    if ~isfield(fparams,'CrossValidation')
        fparams.CrossValidation.Type = 'LNOCV';                              % Cross validation method
        fparams.CrossValidation.N = 100;                                       % Left out trials for validation
    end

    % Model evidence control
    if ~isfield(fparams,'ModelEvidence')
        fparams.ModelEvidence.method = 'none';
        fparams.ModelEvidence.paramLimits = [0.0001 1;...
            0.0001 1;...
            -200 200;...
            0 1];
        fparams.ModelEvidence.integrationMethod = 'quad';
        fparams.ModelEvidence.integrationOptions.dx = 0.1;
        fparams.ModelEvidence.OpenMind = 1;
    end