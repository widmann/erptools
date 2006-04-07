% eeg_rejdelta() - Mark trials with signal changes exceeding threshold
%                  for rejection
%
% Usage:
%   >> EEG = eeg_rejdelta(EEG, 'key1', value1, 'key2', value2, ...
%                              'keyn', valuen);
%
% Inputs:
%   EEG       - EEGLAB EEG structure
%   'thresh'  - threshold
%
% Optional inputs:
%   'chans'   - channels to take into consideration for rejection
%   'pnts'    - samples to take into consideration for rejection
%
% Output:
%   EEG       - EEGLAB EEG structure
%
% Author: Andreas Widmann, University of Leipzig, 2006
%
% See also: pop_rejepoch

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2006 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

% $Id$

function EEG = eeg_rejdelta(EEG, varargin)

args = struct(varargin{:});

if ~isfield(args, 'thresh') || isempty(args.thresh)
    error('Not enough input arguments.')
end
if ~isfield(args, 'chans') || isempty(args.chans)
    args.chans = 1:EEG.nbchan;
end
if ~isfield(args, 'pnts') || isempty(args.pnts)
    args.pnts = 1:EEG.pnts;
end

if ~isfield(EEG, 'reject') || ~isfield(EEG.reject, 'rejthreshE') || isempty(EEG.reject.rejthreshE),
    EEG.reject.rejthreshE = zeros(EEG.nbchan, EEG.trials);
end

EEG.reject.rejthreshE(args.chans, :) = EEG.reject.rejthreshE(args.chans, :) | squeeze(max(EEG.data(args.chans, args.pnts, :), [], 2) - min(EEG.data(args.chans, args.pnts, :), [], 2) > args.thresh);
EEG.reject.rejthresh = any(EEG.reject.rejthreshE);

disp([num2str(length(EEG.reject.rejthresh(EEG.reject.rejthresh))) '/' num2str(EEG.trials) ' trials marked for rejection.'])
