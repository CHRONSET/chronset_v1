% -cleans up following compilation (linux, OSX only)
%
% copyright 2009-2012 Blair Armstrong, Christine Watson, David Plaut
%
%    This file is part of SOS
%
%    SOS is free software: you can redistribute it and/or modify
%    it for academic and non-commercial purposes
%    under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.  For commercial or for-profit
%    uses, please contact the authors (sos@cnbc.cmu.edu).
%
%    SOS is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with SOS (see COPYING.txt).
%    If not, see <http://www.gnu.org/licenses/>.



% removes the additional files created when generating the executable
% similar functionality to 'make clean'

if exist('./bin/sos.prj','file'); delete('./bin/sos.prj'); end;
if exist('./bin/sos_delay_load.c','file'); delete('./bin/sos_delay_load.c'); end;
if exist('./bin/sos_main.c','file'); delete('./bin/sos_main.c'); end;
if exist('./bin/sos_mcc_component_data.c','file'); delete('./bin/sos_mcc_component_data.c'); end;
if exist('./bin/mccExcludedFiles.log','file'); delete('./bin/mccExcludedFiles.log'); end;
if exist('./bin/readme.txt','file'); delete('./bin/readme.txt'); end;
if exist('./bin/sos_linux_main.c','file'); delete('./bin/sos_linux_main.c'); end;
if exist('./bin/sos_linux_mcc_component_data.c','file'); delete('./bin/sos_linux_mcc_component_data.c'); end;
if exist('./bin/sos_linux.prj','file'); delete('./bin/sos_linux.prj'); end;