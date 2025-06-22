function [H, band_H] = StateSpaceEntropyRatePerSource(X, Fs, downsampling, band, varmomax)

% STATESPACEENTROPYRATEPERSOURCE Computes entropy rate per source using a state-space model.
%
% Inputs:
%   X           	- 3D data matrix of size [D, T, M] (sources x time points x trials).
%   Fs          	- Sampling frequency in Hz.
%   downsampling 	- 'yes' or 'no', whether to downsample to a max of 200 Hz (default: 'yes').
%   band        	- Kx2 matrix specifying frequency bands for decomposition (optional).
%   varmomax    	- Maximum VAR model order (default: 20).
%
% Outputs:
%   H           	- Entropy rate per source (1xD vector).
%   band_H      	- Entropy rate decomposed into frequency bands per source (optional).
%
% Pedro Mediano, Luka Dubravica, Modified 2025


%% Parameter checks and initialisation

if isempty(X) || ~isnumeric(X)
  	error('Input must be a numeric 1D, 2D, or 3D data matrix.');
end
if nargin < 2 || isempty(Fs)
  	error('Not enough input arguments. You must provide a sampling frequency.');
end
if nargin < 3 || isempty(downsampling), downsampling = 'yes'; end
if nargin < 4 || isempty(band),         band         = [];    end
if nargin < 5 || isempty(varmomax),     varmomax     = 20;    end

[D, T, M] = size(X);
if T < D
  	warning('Data has more dimensions than time points. Consider transposing it.');
end

if strcmp(downsampling, 'yes') && Fs > 200
  	k = floor(Fs/200);
  	Fs = Fs/k;
  	X = X(:,1:k:end,:);
end

if nargout < 2 && ~isempty(band)
  	error(['Not enough output arguments. When using frequency decomposition, ', ...
        'please use\\[H, band_H] = StateSpaceEntropyRatePerSource(...)']);
end


%% Load Octave and MVGC

% Load Octave packages if relevant
if exist('OCTAVE_VERSION', 'builtin') && ~exist('fcdf', 'file')
  	pkg('load', 'statistics');
end

% Silently add MVGC to current path
p = mfilename('fullpath');
addpath(strrep(p, 'StateSpaceEntropyRatePerSource', 'private/mvgc_v2.0'));
evalc('mvgc_startup;');


%% Loop over sources and compute entropy rate

% Standardize the data
X = demean(X, true);

H = nan(1, D);  % Initialize vector to store entropy rate per source
band_H = nan(size(band,1), D); % Initialize band_H to NaN per source per band

for d = 1:D
  	y = X(d, :, :); % Extract data for the current source

  	% Entropy function for a multivariate normal
  	% MVGC logdet() is faster, safer and more accurate than log(det())
  	H_fun = @(C) 0.5 * logdet(2 * pi * exp(1) * C);

  	% Select VAR model order
  	[varmoaic, varmobic, varmohqc, varmolrt] = tsdata_to_varmo(y, varmomax, 'LWR', [], false);

  	% Select SS model order
  	if varmohqc < 1
    	% Force the algorithm to fit a SS model even if varmo gives up
    	varmohqc = 1;
    	ssmo = 2;
  	else
    	[ssmo, ~] = tsdata_to_sssvc(y, 2*varmohqc, [], []);
  	end

	% Fit SS model
	% 2*AIC is Bauer's recommendation for past/future window... personally I prefer
	% 2*HQC (HQC = Hannan-Quinn IC, generally sits between AIC and BIC)
	failed = false;
	try
		[A, C, K, V] = tsdata_to_ss(y, 2*varmohqc, ssmo);
		assert(all(size(V) == [1, 1]));  % Ensure model is 1D
		info = ss_info(A, C, K, V, 0);
		% Ignore negative multi-info error (bit 7)
		info.error = bitset(info.error, 7, 0);
		if info.error
			failed = true;
		end
	catch
		failed = true;
	end

	% Calculate the entropy rate for the current source
	if ~failed
		H(d) = H_fun(V);
	else
		H(d) = nan;
	end

  	%% Decompose entropy rate into bands, if requested
	if isempty(band) || failed
		continue;  % Skip to next source
	end

	% Compute CPSD
	fres = 1000; % Resolution in frequency space
	S = ss_to_cpsd(A,C,K,V,fres);
	H_freq = shiftdim(arrayfun(@(i) H_fun(S(:,:,i)),1:(fres+1)), -1);

	for j = 1:size(band,1)
		% Ensure band limits are finite
		band_j = band(j,:);
		if isinf(band_j(2))
			band_j(2) = floor(Fs/2.0);
		end
		% Store entropy rate for the current band
		band_H(j,d) = bandlimit(H_freq, 3, Fs, band_j);
	end

end
