% pop_bdfread() - Read BioSemi BDF file
%
% Usage:
%   >> [EEG, com, Hdr] = pop_bdfread; % pop-up window mode
%   >> [EEG, com, Hdr] = pop_bdfread('parameter1', value1, ...
%                                    'parameter2', value2, ...
%                                    'parametern', valuen);
%
% Optional inputs:
%   'filename'    - string filename
%   'pathname'    - string pathname
%   'chans'       - vector of integers channels to read
%   'statusChan'  - scalar integer status channel {default channel of
%                   type 'Triggers and Status' in header}
%   'holdValue'   - scalar integer hold value {default 0}
%
% Outputs:
%   EEG           - EEGLAB EEG structure
%   com           - history string
%   Hdr           - header structure
%
% Author: Andreas Widmann, University of Leipzig, 2006

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

% $Id: pop_bdfread.m 9 2009-12-22 15:11:14Z widmann $

function [EEG, com, Hdr] = pop_bdfread(varargin)

EEG = [];
com = '';
Hdr = [];

if nargin < 1
    [Arg.filename Arg.pathname] = uigetfile2('*.bdf;*.BDF', 'Select BioSemi BDF file -- pop_bdfread()');
    if Arg.filename == 0, return, end
else
    Arg = struct(varargin{:});
end

% Open file
if ~isfield(Arg, 'filename') || isempty(Arg.filename)
    error('Not enough input arguments.')
end
if ~isfield(Arg, 'pathname')
    Arg.pathname = [];
end
[fid, message] = fopen(fullfile(Arg.pathname, Arg.filename));
if fid == -1
    error(message)
end
disp(['pop_bdfread(): loading file ' fullfile(Arg.pathname, Arg.filename) '...'])

% Arguments
if ~isfield(Arg, 'holdValue') || isempty(Arg.holdValue)
    Arg.holdValue = 0;
end

% Header
Hdr.fileFormat = fread(fid, 8, '*char')';
Hdr.subjectId = fread(fid, 80, '*char')';
Hdr.recordId = fread(fid, 80, '*char')';
Hdr.startDate = fread(fid, 8, '*char')';
Hdr.startTime = fread(fid, 8, '*char')';
Hdr.nBytes = str2double(fread(fid, 8, '*char')');
Hdr.dataFormat = fread(fid, 44, '*char')';
Hdr.nBlocks = str2double(fread(fid, 8, '*char')');
Hdr.blockDur = str2double(fread(fid, 8, '*char')');
Hdr.nChans = str2double(fread(fid, 4, '*char')');

% Channel headers
Hdr.labels = cellstr(fread(fid, [16 Hdr.nChans], '*char')');
Hdr.type = cellstr(fread(fid, [80 Hdr.nChans], '*char')');
Hdr.unit = cellstr(fread(fid, [8 Hdr.nChans], '*char')');
Hdr.physMin = str2num(fread(fid, [8 Hdr.nChans], '*char')');
Hdr.physMax = str2num(fread(fid, [8 Hdr.nChans], '*char')');
Hdr.digMin = str2num(fread(fid, [8 Hdr.nChans], '*char')');
Hdr.digMax = str2num(fread(fid, [8 Hdr.nChans], '*char')');
Hdr.filter = cellstr(fread(fid, [80 Hdr.nChans], '*char')');
Hdr.nSamples = str2num(fread(fid, [8 Hdr.nChans], '*char')');
Hdr.reserved = cellstr(fread(fid, [32 Hdr.nChans], '*char')');

% Consistent sampling rate
Hdr.nSamples = unique(Hdr.nSamples);
if length(Hdr.nSamples) > 1
    error('Inconsistent number of samples per channel.')
end

% Scale
Hdr.gain = single(diag((Hdr.physMax - Hdr.physMin) ./ (Hdr.digMax - Hdr.digMin)));
Hdr.offset = Hdr.physMin - Hdr.gain * Hdr.digMin;
Hdr.offset = single(Hdr.offset(:, ones(1, Hdr.nSamples)));
int24max = single(2 ^ 23 - 1);
int24sign = single(2 ^ 24);
castArray = single(pow2([0 8 16]));

% Status channel
if ~isfield(Arg, 'statusChan')
    Arg.statusChan = strmatch('Triggers and Status', Hdr.type);
end
Hdr.isStatus = ismember(1 : Hdr.nChans, Arg.statusChan);
Hdr.isData = ~Hdr.isStatus;

% Channels
if isfield(Arg, 'chans') && ~isempty(Arg.chans)
    Hdr.isData = Hdr.isData & ismember(1 : Hdr.nChans, Arg.chans);
end

% EEG structure
try
    EEG = eeg_emptyset;
catch
end
EEG.setname = 'BDF file';
EEG.filename = Arg.filename;
EEG.filepath = Arg.pathname;
EEG.pnts = Hdr.nSamples * Hdr.nBlocks;
EEG.nbchan = length(find(Hdr.isData));
EEG.trials = 1;
EEG.srate = Hdr.nSamples / Hdr.blockDur;
EEG.xmin = 0;
EEG.xmax = (EEG.pnts - 1) / EEG.srate;
[EEG.chanlocs(1 : EEG.nbchan).labels] = deal(Hdr.labels{Hdr.isData});
[EEG.chanlocs(1 : EEG.nbchan).type] = deal(Hdr.type{Hdr.isData});
temp = num2cell(1 : EEG.nbchan);
[EEG.chanlocs(1 : EEG.nbchan).datachan] = deal(temp{:});
EEG.comments = ['Original file: ' fullfile(Arg.pathname, Arg.filename)];
EEG.ref = 'common';

% Allocate memory
EEG.data = zeros(EEG.nbchan, EEG.pnts, 'single');
status = zeros(length(find(Hdr.isStatus)), EEG.pnts, 'single');

% Initialize progress indicator
nSteps = 20;
step = 0;
fprintf(1, 'pop_bdfread(): |');
strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '|   0%%']);
tic

for iBlock = 1 : Hdr.nBlocks

    % Read, cast, and reshape data
    buf = reshape(castArray * single(fread(fid, [3, Hdr.nSamples * Hdr.nChans], '*uint8')), [Hdr.nSamples Hdr.nChans])';

    % Unsigned to signed
    isNeg = Hdr.isData(ones(1, Hdr.nSamples), :)' & buf > int24max;
    buf(isNeg) = buf(isNeg) - int24sign;

    % Scale
    buf = Hdr.gain * buf + Hdr.offset;

    % Assign
    EEG.data(:, Hdr.nSamples * (iBlock - 1) + 1 : Hdr.nSamples * (iBlock - 1) + Hdr.nSamples) = buf(Hdr.isData, :);
    status(:, Hdr.nSamples * (iBlock - 1) + 1 : Hdr.nSamples * (iBlock - 1) + Hdr.nSamples) = buf(Hdr.isStatus, :);

    % Update progress indicator
    [step, strLength] = mywaitbar(iBlock, Hdr.nBlocks, step, nSteps, strLength);
    
end

% Close file
fclose(fid);

% Deinitialize progress indicator
fprintf(1, '\n')

% Event structure
if size(status, 1) == 1

    status = double(status);

    % CMS
    cmsArray = ~bitget(status, 21);
    if any(cmsArray)
        disp('Warning: CMS out of range; removing affected region(s) from data.')
        EEG.data(:, cmsArray) = [];
        EEG.pnts = size(EEG.data, 2);
        EEG.xmax = (EEG.pnts - 1) / EEG.srate;
        status(cmsArray) = [];
        
        % Insert boundary
        temp = diff([0 cmsArray]) == -1;
        cmsArray = temp(~cmsArray);
    end

    % Boundaries
    dcArray = bitget(status, 17);
    dcArray(1) = 1; % Mandatory boundary event at file start
    dcArray(~diff([0 dcArray])) = 0;
    latArray = num2cell(find(dcArray | cmsArray));
    typeArray(1 : length(latArray)) = {'boundary'};

    % Triggers
    trigArray = rem(status, 2 ^ 16) - Arg.holdValue;
    trigArray(~diff([0 trigArray])) = 0;
    if any(trigArray ~= 0)
        latArray = [latArray num2cell(find(trigArray))];
        typeArray = [typeArray strtrim(cellstr(num2str(trigArray(trigArray ~= 0)', '%d')))'];
    end

    % Sort by latency
    [foo, order] = sort([latArray{:}]);
    [EEG.event(1 : length(typeArray)).type] = deal(typeArray{order});
    [EEG.event.latency] = deal(latArray{order});
    
    % Urevent structure
    EEG.urevent = EEG.event;
    temp = num2cell(1 : length(EEG.urevent));
    [EEG.event.urevent] = deal(temp{:});

end

% History string
com = ['EEG = pop_bdfread(' arg2str(Arg) ');'];

end

function [step, strLength] = mywaitbar(compl, total, step, nSteps, strLength)

progStrArray = '/-\|';
tmp = floor(compl / total * nSteps);
if tmp > step
    fprintf(1, [repmat('\b', 1, strLength) '%s'], repmat('=', 1, tmp - step))
    step = tmp;
    ete = ceil(toc / step * (nSteps - step));
    strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '%s %3d%%, ETE %02d:%02d'], progStrArray(mod(step - 1, 4) + 1), floor(step * 100 / nSteps), floor(ete / 60), mod(ete, 60));
end

end

function str = arg2str(arg)

if isstruct(arg)
    str = '';
    fieldArray = fieldnames(arg);
    arg = struct2cell(arg);
    for iArg = 1:length(arg)
        if ischar(arg{iArg})
            str = [str '''' fieldArray{iArg} ''', ''' arg{iArg} ''''];
        elseif isnumeric(arg{iArg}) || islogical(arg{iArg})
            str = [str '''' fieldArray{iArg} ''', ' mat2str(arg{iArg})];
        elseif iscell(arg{iArg})
            str = [str '''' fieldArray{iArg} ''', ' arg2str(arg{iArg})];
        end
        if iArg < length(arg)
            str = [str ', '];
        end
    end
elseif iscell(arg)
    str = '{';
    for iArg = 1:length(arg)
        if ischar(arg{iArg})
            str = [str '''' arg{iArg} ''''];
        elseif isnumeric(arg{iArg}) || islogical(arg{iArg})
            str = [str mat2str(arg{iArg})];
        elseif iscell(arg{iArg})
            str = [str arg2str(arg{iArg})];
        end
        if iArg < length(arg)
            str = [str ' '];
        else
            str = [str '}'];
        end
    end
end

end
