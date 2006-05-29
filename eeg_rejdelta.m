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
%   'win'     - vector time window to take into consideration for
%               rejection (ms) {default [EEG.xmin EEG.xmax] * 1000}
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

Arg = struct(varargin{:});

if ~isfield(Arg, 'thresh') || isempty(Arg.thresh)
    error('Not enough input arguments.')
end
if ~isfield(Arg, 'chans') || isempty(Arg.chans)
    Arg.chans = 1:EEG.nbchan;
end
if ~isfield(Arg, 'win') || isempty(Arg.win)
    Arg.win = [EEG.xmin EEG.xmax] * 1000;
end
Arg.win = round((Arg.win / 1000 - EEG.xmin) * EEG.srate + 1);
Arg.pnts = Arg.win(1):Arg.win(2);

if ~isfield(EEG, 'reject') || ~isfield(EEG.reject, 'rejthreshE') || isempty(EEG.reject.rejthreshE),
    EEG.reject.rejthreshE = zeros(EEG.nbchan, EEG.trials);
end

EEG.reject.rejthreshE(Arg.chans, :) = EEG.reject.rejthreshE(Arg.chans, :) | squeeze(max(EEG.data(Arg.chans, Arg.pnts, :), [], 2) - min(EEG.data(Arg.chans, Arg.pnts, :), [], 2) > Arg.thresh);
EEG.reject.rejthresh = any(EEG.reject.rejthreshE);

disp([num2str(length(EEG.reject.rejthresh(EEG.reject.rejthresh))) '/' num2str(EEG.trials) ' trials marked for rejection.'])
