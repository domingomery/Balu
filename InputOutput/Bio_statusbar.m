%Display a status/progress bar and inform about the elapsed
%as well as the remaining time (linear estimation).
%
%Synopsis:
%
%  f = Bio_statusbar
%      Get all status/progress bar handles.
%
%  f = Bio_statusbar(title)
%      Create a new status/progress bar. If title is an empty
%      string, the default 'Progress ...' will be used.
%
%  f = Bio_statusbar(title,f)
%      Reset an existing status/progress bar or create a new
%      if the handle became invalid.
%
%  f = Bio_statusbar(done,f)
%      For 0 < done < 1, update the progress bar and the elap-
%      sed time. Estimate the remaining time until completion.
%      On user abort, return an empty handle.
%
%  v = Bio_statusbar('on')
%  v = Bio_statusbar('off')
%      Set default visibility for new statusbars and return
%      the previous setting.
%
%  v = Bio_statusbar('on',f)
%  v = Bio_statusbar('off',f)
%      Show or hide an existing statusbar and return the last
%      visibility setting.
%
%  delete(Bio_statusbar)
%      Remove all status/progress bars.
%
%  drawnow
%      Refresh all GUI windows.
%
% Example:
%     f=Bio_statusbar('Wait some seconds ...');
%     for p=0:0.01:1
%        pause(0.2);
%        if isempty(Bio_statusbar(p,f))
%           break;
%        end
%     end
%     if ishandle(f)
%        delete(f);
%     end
%
% Copyright (c) 2004, Marcel Leutenegger
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the Ecole Polytechnique F??d??rale de Lausanne,
%       Laboratoire d'Optique Biom??dicale nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function f=Bio_statusbar(p,f)
persistent visible;
%vr2014b = strcmp(version('-release'),'2014b');
vr_new  = versionyear()>=2014;

if nargin < nargout           % get handles
    o='ShowHiddenHandles';
    t=get(0,o);
    set(0,o,'on');
    f=findobj(get(0,'Children'),'flat','Tag','StatusBar');
    set(0,o,t);
    return;
end
curtime=86400*now;
if vr_new
    tt = nargin == 2 && isa(check(f),'matlab.ui.Figure');
else
    tt = nargin == 2 && check(f);
end
if nargin & ischar(p)
    if isequal(p,'on') | isequal(p,'off')
        if nargin == 2
            if check(f)             % show/hide
                v=get(f,'Visible');
                set(f,'Visible',p);
            end
        else
            v=visible;
            visible=p;              % default
            if ~strcmp(v,'off')
                v='on';
            end
        end
        if nargout
            f=v;
        end
    else
        if nargin == 2 & check(f)  % reset
            modify(f,'Line','XData',[4 4 4]);
            modify(f,'Rect','Position',[4 54 0.1 22]);
            modify(f,'Done','String','0');
            modify(f,'Time','String','0:00:00');
            modify(f,'Task','String','0:00:00');
        else
            f=create(visible);      % create
        end
        if p
            set(f,'Name',p);
        end
        set(f,'CloseRequestFcn','set(gcbo,''UserData'',-abs(get(gcbo,''UserData'')));','UserData',[curtime curtime 0]);
    end
    drawnow;
elseif tt % update
    t=get(f,'UserData');
    if any(t < 0)              % confirm
        if p >= 1 | isequal(questdlg({'Are you sure to stop the execution now?',''},'Abort requested','Stop','Resume','Resume'),'Stop')
            delete(f);
            f=[];                % interrupt
            return;
        end
        t=abs(t);
        set(f,'UserData',t);    % continue
    end
    p=min(1,max([0 p]));
    %
    % Refresh display if
    %
    %  1. still computing
    %  2. computation just finished
    %    or
    %     more than a second passed since last refresh
    %    or
    %     more than 0.4% computed since last refresh
    %
    if any(t) & (p >= 1 | curtime-t(2) > 1 | p-t(3) > 0.004)
        set(f,'UserData',[t(1) curtime p]);
        t=round(curtime-t(1));
        h=floor(t/60);
        modify(f,'Line','XData',[4 4+248*p 4+248*p]);
        modify(f,'Rect','Position',[4 54 max(0.1,248*p) 22]);
        modify(f,'Done','String',sprintf('%u',floor(p*100+0.5)));
        modify(f,'Time','String',sprintf('%u:%02u:%02u',[floor(h/60);mod(h,60);mod(t,60)]));
        if p > 0.05 | t > 60
            t=ceil(t/p-t);
            if t < 1e7
                h=floor(t/60);
                modify(f,'Task','String',sprintf('%u:%02u:%02u',[floor(h/60);mod(h,60);mod(t,60)]));
            end
        end
        if p == 1
            set(f,'CloseRequestFcn','delete(gcbo);','UserData',[]);
        end
        drawnow;
    end
end
if ~nargout
    clear;
end


%Check if a given handle is a progress bar.
%
function f=check(f)
%vr2014b = strcmp(version('-release'),'2014b');



if versionyear>=2014
    if isa(f,'matlab.ui.Figure') && ishandle(f(1)) && isequal(get(f(1),'Tag'),'StatusBar')
        f=f(1);
    else
        f=[];
    end
else
    if f && ishandle(f(1)) && isequal(get(f(1),'Tag'),'StatusBar')
        f=f(1);
    else
        f=[];
    end
end

%Create the progress bar.
%
function f=create(visible)
if ~nargin | isempty(visible)
    visible='on';
end
s=[256 80];
t=get(0,'ScreenSize');
f=figure('DoubleBuffer','on','HandleVisibility','off','MenuBar','none','Name','Progress ...','IntegerHandle','off','NumberTitle','off','Resize','off','Position',[floor((t(3:4)-s)/2) s],'Tag','StatusBar','ToolBar','none','Visible',visible);
a.Parent=axes('Parent',f,'Position',[0 0 1 1],'Visible','off','XLim',[0 256],'YLim',[0 80]);
%
%Horizontal bar
%
rectangle('Position',[4 54 248 22],'EdgeColor','white','FaceColor',[0.6 0.6 0.6],a);
line([4 4 252],[55 76 76],'Color',[0.5 0.5 0.5],a);
rectangle('Position',[4 54 0.1 22],'EdgeColor','white','FaceColor',[0.1 0.3 1],'Tag','Rect',a);
line([4 4 4],[54 54 77],'Color',[0.2 0.2 0.2],'Tag','Line',a);
%
%Description texts
%
a.FontWeight='bold';
a.Units='pixels';
a.VerticalAlignment='middle';
text(136,66,1,'%',a);
text(16,36,'Elapsed time:',a);
text(16,20,'Remaining:',a);
text(190,36,'hh:mm:ss',a);
text(190,20,'hh:mm:ss',a);
%
%Information texts
%
a.HorizontalAlignment='right';
text(136,66,1,'0',a,'Tag','Done');
text(178,36,'0:00:00',a,'Tag','Time');
text(178,20,'0:00:00',a,'Tag','Task');


%Modify an object property.
%
function modify(f,t,p,v)
set(findobj(f,'Tag',t),p,v);

function y = versionyear()
st = version('-release');
y  = str2num(st(1:4)); %#ok<ST2NM>

