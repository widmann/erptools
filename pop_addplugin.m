% pop_addplugin() - Load EEGLAB plugin
%
% Usage:
%   >> pop_addplugin % pop-up window mode
%   >> pop_addplugin(pathname)
%
% Optional inputs:
%   pathname  - Plugin path
%
% Note:
%   The getstrs subfunction is copy-pasted from the eeglab main function.
%   Copyright belongs to the authors of this function. pop_addplugin might
%   silently break when strings from eeglab main function change in future
%   eeglab version.
%   
% Author: Andreas Widmann, University of Leipzig, 2011
%
% See also:
%   eeglab

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

function pop_addplugin(pathname)

% Arguments
if nargin < 1 || isempty(pathname)
    pathname = uigetdir;
end
if exist(pathname, 'dir') ~= 7
    error('No such directory.')
end

% Find EEGLAB window
eegFigWin = findobj('Tag', 'EEGLAB');
if isempty(eegFigWin)
    error('EEGLAB figure window not found.')
end

% Strings
[trystrs, catchstrs] = getstrs;

% Find plugin function
PluginFunc = dir(fullfile(pathname, 'eegplugin*.m'));
if isempty(PluginFunc)
    error('No eegplugin function found.')
end

% Add to path
addpath(pathname)

% Execute plugin function
for iPluginFunc = 1:length(PluginFunc)
    pluginFuncHandle = str2func(PluginFunc(iPluginFunc).name(1:end - 2));
    vers = pluginFuncHandle(eegFigWin, trystrs, catchstrs);
    disp(['EEGLAB: adding "' vers '" plugin (see >> help ' PluginFunc(iPluginFunc).name(1:end - 2) ')' ]);
end

end

function [trystrs, catchstrs] = getstrs

% checking strings
% ----------------
e_try             = 'try,';
e_catch           = 'catch, eeglab_error; LASTCOM= ''''; clear EEGTMP ALLEEGTMP STUDYTMP; end;';
nocheck           = e_try;
ret               = 'if ~isempty(LASTCOM), if LASTCOM(1) == -1, LASTCOM = ''''; return; end; end;';
check             = ['[EEG LASTCOM] = eeg_checkset(EEG, ''data'');' ret ' eegh(LASTCOM);' e_try];
checkcont         = ['[EEG LASTCOM] = eeg_checkset(EEG, ''contdata'');' ret ' eegh(LASTCOM);' e_try];
checkica          = ['[EEG LASTCOM] = eeg_checkset(EEG, ''ica'');' ret ' eegh(LASTCOM);' e_try];
checkepoch        = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'');' ret ' eegh(LASTCOM);' e_try];
checkevent        = ['[EEG LASTCOM] = eeg_checkset(EEG, ''event'');' ret ' eegh(LASTCOM);' e_try];
checkbesa         = ['[EEG LASTCOM] = eeg_checkset(EEG, ''besa'');' ret ' eegh(''% no history yet for BESA dipole localization'');' e_try];
checkepochica     = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''ica'');' ret ' eegh(LASTCOM);' e_try];
checkplot         = ['[EEG LASTCOM] = eeg_checkset(EEG, ''chanloc'');' ret ' eegh(LASTCOM);' e_try];
checkicaplot      = ['[EEG LASTCOM] = eeg_checkset(EEG, ''ica'', ''chanloc'');' ret ' eegh(LASTCOM);' e_try];
checkepochplot    = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''chanloc'');' ret ' eegh(LASTCOM);' e_try];
checkepochicaplot = ['[EEG LASTCOM] = eeg_checkset(EEG, ''epoch'', ''ica'', ''chanloc'');' ret ' eegh(LASTCOM);' e_try];

% check string and backup old dataset
% -----------------------------------
backup =     [ 'if CURRENTSET ~= 0,' ...
               '    [ ALLEEG EEG ] = eeg_store(ALLEEG, EEG, CURRENTSET, ''savegui'');' ...
               '    eegh(''[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET, ''''savedata'''');'');' ...
               'end;' ];

storecall    = '[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); eegh(''[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);'');';
storenewcall = '[ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET, ''study'', ~isempty(STUDY)+0); eegh(LASTCOM);';
storeallcall = [ 'if ~isempty(ALLEEG) & ~isempty(ALLEEG(1).data), ALLEEG = eeg_checkset(ALLEEG);' ...
                 'EEG = eeg_retrieve(ALLEEG, CURRENTSET); eegh(''ALLEEG = eeg_checkset(ALLEEG); EEG = eeg_retrieve(ALLEEG, CURRENTSET);''); end;' ];

testeegtmp   =  'if exist(''EEGTMP'') == 1, EEG = EEGTMP; clear EEGTMP; end;'; % for backward compatibility
ifeeg        =  'if ~isempty(LASTCOM) & ~isempty(EEG),';
ifeegnh      =  'if ~isempty(LASTCOM) & ~isempty(EEG) & ~isempty(findstr(''='',LASTCOM)),';

% nh = no dataset history
% -----------------------
e_storeall_nh   = [e_catch 'eegh(LASTCOM);' ifeeg storeallcall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_hist_nh       = [e_catch 'eegh(LASTCOM);'];

% same as above but also save history in dataset
% ----------------------------------------------
e_newset        = [e_catch 'EEG = eegh(LASTCOM, EEG);' testeegtmp ifeeg   storenewcall 'disp(''Done.''); end; eeglab(''redraw'');'];
e_store         = [e_catch 'EEG = eegh(LASTCOM, EEG);' ifeegnh storecall    'disp(''Done.''); end; eeglab(''redraw'');'];
e_hist          = [e_catch 'EEG = eegh(LASTCOM, EEG);'];
e_histdone      = [e_catch 'EEG = eegh(LASTCOM, EEG); if ~isempty(LASTCOM), disp(''Done.''); end;' ];

% study checking
% --------------
e_load_study    = [e_catch 'eegh(LASTCOM); if ~isempty(LASTCOM), ALLEEG = ALLEEGTMP; EEG = ALLEEG; CURRENTSET = [1:length(EEG)]; eegh(''CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];''); STUDY = STUDYTMP; CURRENTSTUDY = 1; disp(''Done.''); end; clear ALLEEGTMP STUDYTMP; eeglab(''redraw'');'];

% build structures for plugins
% ----------------------------
trystrs.no_check                 = e_try;
trystrs.check_data               = check;
trystrs.check_ica                = checkica;
trystrs.check_cont               = checkcont;
trystrs.check_epoch              = checkepoch;
trystrs.check_event              = checkevent;
trystrs.check_epoch_ica          = checkepochica;
trystrs.check_chanlocs           = checkplot;
trystrs.check_epoch_chanlocs     = checkepochplot;
trystrs.check_epoch_ica_chanlocs = checkepochicaplot;
catchstrs.add_to_hist            = e_hist;
catchstrs.store_and_hist         = e_store;
catchstrs.new_and_hist           = e_newset;
catchstrs.new_non_empty          = e_newset;
catchstrs.update_study           = e_load_study;

end

