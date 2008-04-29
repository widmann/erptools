% eeg_cdemod() - Complex demodulation
%
% Usage:
%   >> [EEG, com] = eeg_cdemod(EEG, 'key1', value1, 'key2', value2, 
%                                   'keyn', valuen); 
%
% Inputs:
%   EEG       - EEGLAB EEG structure
%   'freq'    - scalar carrier frequency 
%
% Optional inputs:
%   'lowpass' - scalar lowpass frequency (Hz) {default 5}
%   'forder'  - scalar low-pass filter order (mandatory even)
%               {default 4} 
%
% Outputs:
%   EEG       - EEGLAB EEG structure
%   com       - history string
%
% References:
%   [1] Müller, M. M., Keil, A., Kissler, J., Gruber, T. (1999).
%       Suppression of the auditory middle-latency response and evoked
%       gamma-band response in a paired-click paradigm. Experimental
%       Brain Research, 136, 474-479.
%   [2] Regan, D. (1989). Human brain electrophysiology: evoked
%       potentials and evoked magnetic fields in science and medicine.
%       NY:Elsevier.
%
% Author: Andreas Widmann, University of Leipzig, 2007

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2007 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

function [EEG, com] = eeg_cdemod(EEG, varargin)

com = '';
Arg = struct(varargin{:});

% Frequency
if ~isfield(Arg, 'freq') || isempty(Arg.freq)
    error('Not enough input arguments.')
end

% Filter order
if ~isfield(Arg, 'forder') || isempty(Arg.forder)
    Arg.forder = 4;
elseif mod(Arg.forder, 2)
    error('Filter order must be even.') % filtfilt
end

% Low-pass
if ~isfield(Arg, 'lowpass') || isempty(Arg.lowpass)
    Arg.lowpass = 5;
end

% Filter coefficients
[b, a] = butter(Arg.forder / 2, Arg.lowpass / (EEG.srate / 2));

% Demodulation
carArray = exp(2 * pi * i * Arg.freq * (0:EEG.pnts - 1) / EEG.srate);
for iChan = 1:EEG.nbchan
    for iTrial = 1:EEG.trials
        x = double(EEG.data(iChan, :, iTrial)) .* carArray;
        x = filtfilt(b, a, x);
        EEG.data(iChan, :, iTrial) = 2 * sqrt(real(x).^2 + imag(x).^2);
    end
end

% History string
inName = inputname(1);
if isempty(inName)
    inName = 'EEG';
end
com = [inName ' = eeg_cdemod(' inName ', ' arg2str(Arg) ');'];
