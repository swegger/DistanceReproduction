function [squaredBias, variance, RMSE] = BootStrapBiasVariance(s,p,varargin)
%% BootStrapBiasVariance
%
%   [squaredBias, variance, RMSE] = BootStrapBiasVariance(s,p)
%
%   Finds a bootstrap distribution of bias and variance of production as a
%   function of sample time for several conditions. Inputs s and p are
%   required to be cells to accomodate multiple conditions.
%
%   ...  = BootStrapBiasVariance(s,p,'nbootstraps',N)
%   
%   Performs N bootstraps of the data. Defaults to 100.
%
%   ... = BootStrapBiasVariance(...,'nsamps',N)
%
%   Performs each bootstrap on N samples from each condition and sample
%   number. Defaults to the number of samples in each condition and sample.
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'s');
addRequired(Parser,'p');
addParameter(Parser,'nbootstraps',100);
addParameter(Parser,'nsamps',NaN);
addParameter(Parser,'Mask',[]);
addParameter(Parser,'sUnique',NaN);
addParameter(Parser,'b',0);
addParameter(Parser,'g',1);

parse(Parser,s,p,varargin{:});

s = Parser.Results.s;
p = Parser.Results.p;
nbootstraps = Parser.Results.nbootstraps;
nsamps = Parser.Results.nsamps;
Mask = Parser.Results.Mask;
sUnique = Parser.Results.sUnique;
b = Parser.Results.b;
g = Parser.Results.g;

% Find number of conditions
m = length(s);

% Determine Mask if not supplied
if ~iscell(Mask)
    for i = 1:m
        Mask{i} = false(size(s{i}));
    end
end

% Find unique values of s if not supplied
if isnan(sUnique)
    S = cat(find(size(s{1}) == 1),s{:});        % Concatentate across conditions
    sUnique = unique(S);                        % unique values of S
end


%% Calculate RMSE, Bias and Variance
for i = 1:m
    % Bootstrap bias and variance
    for j = 1:nbootstraps
        es = [];
        for k = 1:length(sUnique)
            tempP = p{i}(~Mask{i} & s{i} == sUnique(k));
            es = [es; tempP - sUnique(k)];
            if isnan(nsamps)
                nsampstemp = length(tempP);
            else
                nsampestemp = nsamps;
            end
            inds = ceil(length(tempP)*rand(nsampstemp,1));
            tempm(k,:) = mean(tempP(inds));
            tempstd(k,:) = std(tempP(inds));
        end
        RMSE(j,i) = sqrt(mean(es.^2));
        squaredBias(j,i) = mean((tempm/g - b - sUnique).^2);
        variance(j,i) = mean(tempstd.^2);
    end
end
