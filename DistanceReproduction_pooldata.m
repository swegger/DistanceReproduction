function [ds1, ds2, ds, dp, dsMin, dsMax, Trials, correct, target_x, target_y, AssignedTrials, N] = DistanceReproduction_pooldata(d,varargin)
%% RSCCGmwk_pooldata
%
%   [ts tp tsmin tsmax] = RSCCGmwk_pooldata(d,...)
%   Pools the sample and production data across a set of conditions and
%   runs.
%
%   [ts ...] = RSCCGmwk_pooldata(d,'runs',run_vec)
%   Pools the data across all conditions for the runs specified by run_vec.
%
%   [ts ...] = RSCCGmwk_pooldata(d,'runs',run_vec,'trialtype',type_vec)
%   Pools data for the trial type in type_vec across all runs specified
%   by run_vec
%
%   [ts ...] = RSCCGmwk_pooldata(d,'trialWin',trialWin)
%   Pools data across all runs and trial types for trials trailWin(1)
%   through trialWin(2)
%
%%

% Defaults
target_fun = @(target_x,target_y)(target_x == target_x);

% Parse inputs
p = inputParser;
addRequired(p,'d');
addParameter(p,'flash_loc',-1:1)
addParameter(p,'runs',NaN)
addParameter(p,'trialWin',[1 Inf])
addParameter(p,'TargetLocationFunction',target_fun)
addParameter(p,'CorrectValues',[0 1])
addParameter(p,'Distance_N',[1 2])

parse(p,d,varargin{:})

d = p.Results.d;
flash_loc = p.Results.flash_loc;
runs = p.Results.runs;
trialWin = p.Results.trialWin;
CorrectValues = p.Results.CorrectValues;
TargetLocationFunction = p.Results.TargetLocationFunction;
Distance_N = p.Results.Distance_N;

if isnan(runs)
    runs = 1:d.runs;
end
    

% Initialize
ds1 = [];
ds2 = [];
ds = [];
dp = [];
Trials = [];
correct = [];
target_x = [];
target_y = [];
AssignedTrials = [];
N = [];

% Get data
for j = runs
    
    if length(d.dp{j}) == length(d.ds{j}) % Checks to make sure number trials matches for tr and tp
        
        % Align data by trial and flash location, starting from trial number "trialWin(1)" and ending at trial number "trialWin(2)" (or to end if trialWin(2)== Inf)
        flash_locSort = zeros(size(d.flash_loc{j},1),1);
        for k = flash_loc
            flash_locSort(d.flash_loc{j}(:,2) == k) = 1;
        end
        x = d.target_x{j}(:,2);
        y = d.target_y{j}(:,2);
        targetSort = TargetLocationFunction(x,y);
        correctSort = zeros(size(d.flash_loc{j},1),1);
        for k = CorrectValues
            correctSort(d.correct{j}(:,2) == k) = 1;
        end
        Distance_N_sort = any(repmat(d.Distance_N{j}(:,2),1,length(Distance_N)) == repmat(Distance_N(:)',size(d.Distance_N{j}(:,2),1),1),2);
        trials = d.flash_loc{j}(d.flash_loc{j}(:,1) >= trialWin(1) & d.flash_loc{j}(:,1) <= trialWin(2) & targetSort & flash_locSort & correctSort & Distance_N_sort,1);
        
        [~, ~, inds2] = intersect(trials,d.ds{j}(:,1));
        ds1 = [ds1; double(d.ds1{j}(inds2,2))];
        ds2 = [ds2; double(d.ds2{j}(inds2,2))];
        ds = [ds; double(d.ds{j}(inds2,2))];
        dp = [dp; double(d.dp{j}(inds2,2))];
        Trials = [Trials; trials(:)];
        correct = [correct; double(d.correct{j}(inds2,2))];
        target_x = [target_x; double(d.target_x{j}(inds2,2))];
        target_y = [target_y; double(d.target_y{j}(inds2,2))];
        AssignedTrials = [AssignedTrials; d.trialNumber{j}(inds2,2)];
        N = [N; d.Distance_N{j}(inds2,2)];
        
    else
        error(['Number of tr and tp trials does not agree for subject ' snames ' run ' num2str(j)])
    end
end

% Find minimum and maximum sample times
dsMin = min(ds);
dsMax = max(ds);