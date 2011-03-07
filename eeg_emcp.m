% eeg_emcp() - Basic eye movement correction with raw average subtraction
%
% Usage:
%   >> EEG = eeg_emcp(EEG, 'parameter1', value1, ...
%                     'parameter2', value2, ...
%                     'parametern', valuen);
%   >> [EEG, corrArray, com] = eeg_emcp(EEG, 'parameter1', value1, ...
%                                       'parameter2', value2, ...
%                                       'parametern', valuen);
%
% Inputs:
%   EEG           - EEGLAB EEG structure
%   'eeg'         - vector EEG channels
%   'eog'         - vector EOG channels
%
% Optional inputs:
%   'avgsub'      - flag subtract raw average {default true}
%
% Outputs:
%   EEG           - EEGLAB EEG structure
%   corrArray     - propagation coefficient matrix used for correction
%   com           - history string
%
% References:
%   Gratton, G., Coles, M. G. H., & Donchin, E. (1983). A new method for
%   off-line removal of ocular artifact. Electroencephalography and
%   Clinical Neurophysiology, 55, 468?484.
%
% Note:
%  Bipolarize EOG channels before correction where applicable. No blink
%  detection is (yet?) implemented.
%
% Author: Andreas Widmann, University of Leipzig, 2011

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2011 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG, corrArray, com] = eeg_emcp(EEG, varargin)

com = '';
Arg = struct(varargin{:});

[nChans, nPnts, nTrials] = size(EEG.data);

% Raw average subtraction
if isfield(Arg, 'avgsub') && Arg.avgsub && nTrials > 1
    data = EEG.data - repmat(mean(EEG.data, 3), [1 1 nTrials]);
else
    data = EEG.data;
end

% Reshape data
data = reshape(data, [nChans nPnts * nTrials])';

% Covariance matrix
covArray = cov(data);

% Correction coefficients
corrArray = eye(nChans);
corrArray(Arg.eog, Arg.eeg) = -inv(covArray(Arg.eog, Arg.eog)) * covArray(Arg.eog, Arg.eeg);

% Correction
EEG.data = reshape(EEG.data, [nChans nPnts * nTrials])';
EEG.data = EEG.data * corrArray;
EEG.data = reshape(EEG.data', [nChans nPnts nTrials]);
