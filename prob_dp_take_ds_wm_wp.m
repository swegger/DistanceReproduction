function [pdist, modelprob, loglike, dss, like] = prob_dp_take_ds_wm_wp(dp,ds,wm,wp,N,varargin)
% prob_dp_take_ds_wm_wp
%
%   [prob, loglike] = prob_dp_take_ds_wm_wp(dp,ds,wm,wp)
%   
%   Determines the distribution of production times as a function of sample
%   times under the BLS model with a given wm and wp. Allows for conflict
%   in sample times.
%
%%

% Defaults
estimator_default.type = 'BLS';

% Parse inputs
p = inputParser;

addRequired(p,'dp');
addRequired(p,'ds');
addRequired(p,'wm');
addRequired(p,'wp');
addRequired(p,'N');
addParameter(p,'dpvec',[]);
addParameter(p,'dss',[]);
addParameter(p,'dsMinMax',[]);
addParameter(p,'estimator',estimator_default)

parse(p,dp,ds,wm,wp,N,varargin{:})

dp = p.Results.dp;
ds = p.Results.ds;
wm = p.Results.wm;
wp = p.Results.wp;
N = p.Results.N;
dpvec = p.Results.dpvec;
dss = p.Results.dss;
dsMinMax = p.Results.dsMinMax;
estimator = p.Results.estimator;

if isempty(dpvec)
    % Generate a vector of points to find the probability of a production time
    dpvec = linspace(min(dp),max(dp),100);
end

if isempty(dss)
    dss = unique(ds,'rows');
end

if isempty(dsMinMax)
    dsMinMax = [min(ds(:)) max(ds(:))];
end
    

if size(ds,2) == 1 && N > 1
    % If number of columns of ds is 1 and N > 1, assume no conflict
    ds = repmat(ds,1,N);
    dss = repmat(dss,1,N);
elseif size(ds,2) ~= N
    error('The number of columns in ds must be equal to N or 1.')
end

% Determine the distribution of production times for each combination of sample intervals
for i = 1:length(dss)
    inds = find(sum(ds == repmat(dss(i,:),size(ds,1),1),2) == N);
    n = hist(dp(inds),dpvec);
    pdist(:,i) = n/sum(n);
end

% Determine expected distributions under each model for each combination of sample intervals
for i = 1:length(dss)
    dstemp = dss(i,:);
    modelprob(:,i) = p_dp(dstemp,dpvec(:),wm,wp,N,dsMinMax(1),dsMinMax(2),estimator);
end

% Determine the likelihood of the model given the data
for i = 1:size(ds,1)
    like(i,:) = p_dp(ds(i,:),dp(i),wm,wp,N,dsMinMax(1),dsMinMax(2),estimator);
end
loglike = sum(log(like(:)));

%% Functions
function p = p_dp(ds,dp,wm,wp,N,dsmin,dsmax,estimator)
%% For generating the probability of each dp
% [y1, y2] = ndgrid(ds,dp);
% Y = [y1(:) y2(:)];

Y{1} = ds;
Y{2} = dp;

functionHandle = @(tm,Y)(integrand(tm,Y,wm,wp,dsmin,dsmax,estimator));
options.dx = 0.02;
p = ndintegrate(functionHandle,repmat([mean(ds)-5*mean(ds)*wm mean(ds)+5*mean(ds)*wm],N,1),'method','quad','options',options,'ExtraVariables',Y);


function out = integrand(dm,Y,wm,wp,dsmin,dsmax,estimator)
%% For integration of across measurements in generating the p(dp|ds)

% Unpack extra variables
% ds = Y(:,1);
% dp = Y(:,2);
ds = Y{1};
dp = Y{2};

% Reshape the elements so that the computations can be performed correctly
% DS = repmat(permute(ds,[3 2 1]),size(tm,1),size(tm,2),1);
% TP = repmat(permute(dp,[3 2 1]),size(tm,1),1,1);
DS = repmat(ds,size(dm,1),1);
DP = repmat(permute(dp,[3 2 1]),size(dm,1),1,1);
method_opts.type = 'quad';
method_opts.dx = 1;
fBLS = ScalarBayesEstimators(dm,wm,dsmin,dsmax,'method',method_opts,'estimator',estimator);
fBLS = repmat(fBLS,1,1,size(DP,3));

% Find the probability for each combination of variables
if size(dm,2) == 1
    p_dp_take_fBLS = (1./sqrt(2*pi*wp^2*fBLS.^2)).*exp(-(DP-fBLS).^2./(2*wp^2*fBLS.^2));
    p_tm_take_ds = (1./sqrt(2*pi*wm^2*DS.^2)).*exp(-sum((dm-DS).^2,2)./(2*wm^2*DS.^2));
elseif size(dm,2) == 2
%     TPmod = TP(:,:,1:end/2);
%     fBLSmod = fBLS(:,:,1:end/2);
%     p_dp_take_fBLS = (1./sqrt(2*pi*wp^2*fBLSmod.^2)).*exp(-(TPmod-fBLSmod).^2./(2*wp^2*fBLSmod.^2));
    p_dp_take_fBLS = (1./sqrt(2*pi*wp^2*fBLS.^2)).*exp(-(DP-fBLS).^2./(2*wp^2*fBLS.^2));
%     DS1 = DS(:,1,1:2:end);
%     DS2 = DS(:,2,2:2:end);
%     p_tm_take_ds = (1./sqrt(2*pi*wm^2*DS1.^2)).*exp(-(repmat(tm(:,1),1,1,size(DS1,3))-DS1).^2./(2*wm^2*DS1.^2)) .* (1./sqrt(2*pi*wm^2*DS2.^2)).*exp(-(repmat(tm(:,2),1,1,size(DS2,3))-DS2).^2./(2*wm^2*DS2.^2)); %(1./(2*pi*wm^2*DS1.*DS2)).*exp(-(repmat(tm(:,1),1,1,size(DS1,3))-DS1).^2./(2*wm^2*DS1.^2) - (repmat(tm(:,2),1,1,size(DS2,3))-DS2).^2./(2*wm^2*DS2.^2));
    p_tm_take_ds = (1./sqrt(2*pi*wm^2*DS(:,1).^2)).*exp(-(dm(:,1)-DS(:,1)).^2./(2*wm^2*DS(:,1).^2)) .* (1./sqrt(2*pi*wm^2*DS(:,2).^2)).*exp(-(dm(:,2)-DS(:,2)).^2./(2*wm^2*DS(:,2).^2));
end

out = permute(p_dp_take_fBLS.*repmat(p_tm_take_ds,1,1,size(DP,3)),[1 3 2]);
