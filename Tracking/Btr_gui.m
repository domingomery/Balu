function varargout = Btr_gui(varargin)
% BTR_GUI M-file for Btr_gui.fig
%      BTR_GUI, by itself, creates a new BTR_GUI or raises the existing
%      singleton*.
%
%      H = BTR_GUI returns the handle to a new BTR_GUI or the handle to
%      the existing singleton*.
%
%      BTR_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BTR_GUI.M with the given input arguments.
%
%      BTR_GUI('Property','Value',...) creates a new BTR_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Btr_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Btr_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Btr_gui

% Last Modified by GUIDE v2.5 13-Jun-2011 15:11:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Btr_gui_OpeningFcn, ...
    'gui_OutputFcn',  @Btr_gui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Btr_gui is made visible.
function Btr_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Btr_gui (see VARARGIN)

% Choose default command line output for Btr_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Btr_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global gdata
gdata.Iseq          = 178*ones(128,512);
gdata.img           = 0;
gdata.tok           = 0;
gdata.allobj        = 0;
gdata.monobj        = 0;
gdata.showdetails   = 0;
gdata.showstructure = 0;
gdata.linkframe     = 0;
gdata.tracks        = 1;
gdata.multilines    = 0;
gdata.nz            = [];
imshow(uint8(178*ones(256,256)))


% --- Outputs from this function are returned to the command line.
function varargout = Btr_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in NextImage.
function NextImage_Callback(hObject, eventdata, handles)
% hObject    handle to NextImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gdata
gdata.img = min([gdata.f.imgmax gdata.img+1]);
ShowImage(handles);

% --- Executes on button press in PreviousImage.
function PreviousImage_Callback(hObject, eventdata, handles)
% hObject    handle to PreviousImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gdata
gdata.img = max([gdata.f.imgmin gdata.img-1]);
ShowImage(handles);

% --- Executes on button press in LoadAll.
function LoadAll_Callback(hObject, eventdata, handles)
% hObject    handle to LoadAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gdata
d = dir('*.mat');
str = {d.name};
[ss,vv] = listdlg('PromptString','Select a file:',...
    'SelectionMode','single',...
    'ListString',str);
if vv
    fstr = d(ss,:).name;
    load(fstr)
    gdata.f = f;
    gdata.op1 = op1;
    gdata.op2 = op2;
    gdata.op3 = op3;
    % fprintf('Loading images %s%d...%d%s\n',f.prefix,f.imgmin,f.imgmax,f.extension)
end
set(handles.printfilename, 'String', fstr);
set(handles.trackstatus, 'String', 'ready');
gdata.img = gdata.f.imgmin;
gdata.tok = 0;
gdata.Iseq = 178*ones(128,512);
ShowImage(handles);

% --- Executes on button press in Track.
function Track_Callback(hObject, eventdata, handles)
% hObject    handle to Track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gdata
if gdata.tok == 0
    set(handles.trackstatus, 'String', 'searching...');
    drawnow
    tic
    [Tf,kp1,kp2,f,P,Iz,Wd,inx] = Btr_detection(gdata.f,gdata.op1,gdata.op2,gdata.op3);

    gdata.Tf  = Tf;
    gdata.kp1 = kp1;
    gdata.kp2 = kp2;
    gdata.f   = f;
    gdata.P   = P;
    gdata.Iz  = Iz;
    gdata.Wd  = Wd;
    gdata.inx = inx;
    t         = round(toc);
    gdata.tok = 1;
    gdata.nz = 1;
    gdata.nobj = size(Wd,3);
    gdata.Iseq = Iz(1:128,:);
    nf = f.imgmax-f.imgmin+1;
    set(handles.trackstatus, 'String', [num2str(t) 'sec']);
    drawnow
    op = gdata.op1;
    A = Btr_siftn(kp1,op);
    B = A;
    gdata.inx.S2 = A;
    op.join_iter = 1;
    for i=3:nf
        B = Btr_join(B,A,op);
        eval(['gdata.inx.S' num2str(i) '= B;']);
    end

    pause(1)
    set(handles.trackstatus, 'String', 'ok');
    drawnow
    gdata.img = gdata.f.imgmin;
    ShowImage(handles);
end
% --- Executes on button press in Select.
function Select_Callback(hObject, eventdata, handles)
% hObject    handle to Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gdata
if gdata.tok
    img0 = gdata.img-gdata.f.imgmin+1;
    I = Bio_loadimg(gdata.f,gdata.f.imgmin);
    [N,M] = size(I);
    dth = sqrt(N*M)/10;
    mi = gdata.mm(1);
    bi = gdata.mm(2);
    mj = gdata.mm(3);
    bj = gdata.mm(4);
    [x,y] = getpts;
    p = [mean(x) mean(y)];
    ip = (p(2)-bi)/mi;
    jp = (p(1)-bj)/mj;
    di = zeros(gdata.nobj,1);
    for i=1:gdata.nobj
        i1=gdata.Wd(img0,1,i);
        i2=gdata.Wd(img0,2,i);
        j1=gdata.Wd(img0,3,i);
        j2=gdata.Wd(img0,4,i);
        im = (i1+i2)/2;
        jm = (j1+j2)/2;
        di(i) = sqrt((im-ip)^2+(jm-jp)^2);
    end
    [dmin,imin] = min(di);
    if dmin<dth
        gdata.nz = imin;
        gdata.Iseq = gdata.Iz((gdata.nz-1)*128+1:128*gdata.nz,:);
        gdata.allobj = 0;
        gdata.monobj = 1;
        gdata.linkframe = 1;
        set(handles.All,'Value',0);
        set(handles.Mono,'Value',1);
        set(handles.Link,'Value',1);
        ShowImage(handles);
    end
end

% --- Executes on button press in NextObject.
function NextObject_Callback(hObject, eventdata, handles)
% hObject    handle to NextObject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gdata
if gdata.tok
    gdata.nz = min([gdata.nz+1 gdata.nobj]);
    gdata.Iseq = gdata.Iz((gdata.nz-1)*128+1:128*gdata.nz,:);
    ShowImage(handles)
end

% --- Executes on button press in PreviousObject.
function PreviousObject_Callback(hObject, eventdata, handles)
% hObject    handle to PreviousObject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gdata
if gdata.tok
    gdata.nz = max([gdata.nz-1 1]);
    gdata.Iseq = gdata.Iz((gdata.nz-1)*128+1:128*gdata.nz,:);
    ShowImage(handles)
end

% --- Executes on button press in All.
function All_Callback(hObject, eventdata, handles)
% hObject    handle to All (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of All
global gdata
if gdata.nz==0
    gdata.nz = 1;
end
gdata.allobj = get(hObject,'Value');
if gdata.allobj
    gdata.monobj=0;
    set(handles.Mono,'Value',0);
end
ShowImage(handles);

% --- Executes on button press in Mono.
function Mono_Callback(hObject, eventdata, handles)
% hObject    handle to Mono (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of Mono
global gdata
if gdata.nz==0
    gdata.nz = 1;
end
gdata.monobj = get(hObject,'Value');
if gdata.monobj
    gdata.allobj=0;
    set(handles.All,'Value',0);
end
ShowImage(handles);


% --- Executes on button press in Link.
function Link_Callback(hObject, eventdata, handles)
% hObject    handle to Link (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of Link
global gdata
gdata.linkframe = get(hObject,'Value');
ShowImage(handles)


% --- Executes on button press in Structure.
function Structure_Callback(hObject, eventdata, handles)
% hObject    handle to Structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of Structure
global gdata
gdata.showstructure = get(hObject,'Value');
if and(gdata.tracks>1,gdata.showstructure)
    gdata.showdetails = 0;
    set(handles.Details,'Value',0);
end
ShowImage(handles)


% --- Executes on button press in Details.
function Details_Callback(hObject, eventdata, handles)
% hObject    handle to Details (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of Details
global gdata
gdata.showdetails = get(hObject,'Value');
if and(gdata.tracks>1,gdata.showdetails)
    gdata.showstructure = 0;
    set(handles.Structure,'Value',0);
end
ShowImage(handles)


function ShowImage(handles)
global gdata
if gdata.img>0

    n1 = 512;          %image I1 is n1xn1
    n2 = round(n1/4);  %image I2 is n2xn1
    nf = gdata.f.imgmax-gdata.f.imgmin+1;
    img0 = gdata.img-gdata.f.imgmin+1;

    I  = Bio_loadimg(gdata.f,gdata.img);
    I1 = 178*ones(n1,n1);
    [h,w] = size(I);
    if w>h
        h1 = round(n1/w*h);
        h2 = round((n1-h1)/2)+1;
        I10 = imresize(I,[h1 n1]);
        I1 (h2:h2+h1-1,:) = I10;
        i1 = h2; i2 = h2+h1-1;
        j1 = 1;  j2 = n1;
    else
        w1 = round(n1/h*w);
        w2 = round((n1-w1)/2)+1;
        I10 = imresize(I,[n1 w1]);
        I1 (:,w2:w2+w1-1,:) = I10;
        i1 = 1;  i2 = n1;
        j1 = w2; j2 = w2+w1-1;
    end
    mi = (i2-i1)/(h-1);
    bi = i2-mi*h;
    mj = (j2-j1)/(w-1);
    bj = j2-mj*w;
    gdata.mm = [mi bi mj bj];
    si = [ num2str(img0) '/' num2str(gdata.f.imgmax-gdata.f.imgmin+1)];
    set(handles.imagenumber , 'String', si);
    set(handles.objectnumber, 'String', ' ');



    cstr = 'bgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykwbgrcmykw';

    objs = or(gdata.monobj,gdata.allobj);

    if and(and(gdata.linkframe,gdata.tok),and(gdata.nz>0,objs))
        I2 = imresize(gdata.Iseq,[n2 n1]);
    else
        I2 = imresize(178*ones(128,512),[n2 n1]);
    end

    imshow(uint8([I1;I2]));
    hold on

    if gdata.tok && (gdata.nz>0)
        if gdata.allobj
            nd = ones(gdata.nobj,1);
        end
        if gdata.monobj
            nd = zeros(gdata.nobj,1);
            nd(gdata.nz) = 1;
        end

        if and(gdata.showstructure,gdata.tracks==1) % structure from motion
            ii = find(gdata.kp1.img == img0);
            x  = gdata.kp1.fra(ii,1);
            y  = gdata.kp1.fra(ii,2);
            plot(mj*x+bj,mi*y+bi,'b.');
        end

        if and(gdata.showdetails,gdata.tracks==1) % fine detection
            ii = find(gdata.kp2.img == img0);
            x  = gdata.kp2.fra(ii,1);
            y  = gdata.kp2.fra(ii,2);
            plot(mj*x+bj,mi*y+bi,'r.');
        end
        if gdata.tracks>1
            ok = 0;
            if gdata.showstructure
                ct = 'b';
                kp = gdata.kp1;
                % op = gdata.op1;
                if gdata.tracks > nf
                    A = gdata.inx.H1;
                else
                    eval(['A = gdata.inx.S' num2str(gdata.tracks) ';']);
                end
                ok = ~isempty(A);

            end
            if and(gdata.showdetails,gdata.tracks<=6)
                ct = 'r';
                gdata.showstructure = 0;
                set(handles.Structure,'Value',0);
                kp = gdata.kp2;
                s  = ' 23rmf';
                eval(['A = gdata.inx.T' s(gdata.tracks) ';']);
                ok = ~isempty(A);
            end

            if ok
                if gdata.multilines
                    n = size(A,1);
                    if n>1
                        for i=1:n
                            Y = A(i,:);
                            Y = Y(Y>0);
                            x = kp.fra(Y,1:2);
                            plot(mj*x(:,1)+bj,mi*x(:,2)+bi,ct)
                            plot(mj*x(:,1)+bj,mi*x(:,2)+bi,[ct '.'])
                            %                            plot(mj*x(:,1)+bj,mi*x(:,2)+bi,[cstr(i) '.'])
                        end
                    end
                else
                    ii = find((kp.img(A(:,1))==img0));
                    n = length(ii);
                    if n>1
                        for i=1:n
                            Y = A(ii(i),:);
                            Y = Y(Y>0);
                            x = kp.fra(Y,1:2);
                            plot(mj*x(:,1)+bj,mi*x(:,2)+bi,ct)
                            %                            plot(mj*x(:,1)+bj,mi*x(:,2)+bi,cstr(i))
                            plot(mj*x(:,1)+bj,mi*x(:,2)+bi,[ct '.'])
                            %                            plot(mj*x(:,1)+bj,mi*x(:,2)+bi,[cstr(i) '.'])
                        end
                    end
                end
            end
        end
        if objs
            i10=n1+1;
            i20=n1+n2;
            j10=(n1/nf)*(img0-1)+1;
            j20=n1/nf*img0;
            if gdata.nz<=8
                lstr = cstr(gdata.nz);
            else
                if gdata.nz<=16
                    lstr = [cstr(gdata.nz) '--'];
                else
                    if gdata.nz<=24
                        lstr = [cstr(gdata.nz) ':'];
                    else
                        lstr = [cstr(gdata.nz) '-.'];
                    end
                end
            end
            if gdata.linkframe
                Bio_plotsquare([i10 i20 j10 j20],lstr)
            end
            sa = [num2str(gdata.nz) '/' num2str(gdata.nobj)];
            set(handles.objectnumber, 'String', sa);
            for i=1:gdata.nobj
                if nd(i)
                    i1=gdata.Wd(img0,1,i);
                    i2=gdata.Wd(img0,2,i);
                    j1=gdata.Wd(img0,3,i);
                    j2=gdata.Wd(img0,4,i);
                    if i<=8
                        lstr = cstr(i);
                    else
                        if i<=16
                            lstr = [cstr(i) '--'];
                        else
                            if i<=24
                                lstr = [cstr(i) ':'];
                            else
                                lstr = [cstr(i) '-.'];

                            end
                        end
                    end
                    i1 = mi*i1+bi;
                    i2 = mi*i2+bi;
                    j1 = mj*j1+bj;
                    j2 = mj*j2+bj;
                    Bio_plotsquare([i1 i2 j1 j2],lstr)

%                     Shadow regions                    
%                     [N1,M1] = size(I1);
%                     I1c = zeros(N1,M1,3);
%                     I1c(:,:,1) = I1;
%                     I1c(:,:,2) = I1;
%                     I1c(:,:,3) = I1;
%                     Is = I1(round(i1:i2),round(j1:j2));
%                     Rs = Bim_segbalu(Is);
%                     figure(1)
%                     imshow(Is,[])
%                     figure(2)
%                     imshow(Rs)
%                     ii = find(Rs==1);
%                     Is(ii) = 255;
%                     I1c(round(i1:i2),round(j1:j2),3) = Is;
%                     figure(3)
%                     imshow(double(I1c)/255)

                    if and(i==gdata.nz,gdata.linkframe)
                        ils = [i2 i10-10 i10-10 i10];
                        jls = [(j1+j2) (j1+j2) (j10+j20) (j10+j20)]/2;
                        plot(jls,ils,lstr)
                        plot(jls(4),ils(4),[lstr(1) 'o'])
                    end

                end
            end
        end
    end
    drawnow
    hold off
    pm = 0;

    if gdata.showstructure
        S = [' 1 view    ';' 2 views   '];
        for i=3:nf
            S  = [S;sprintf('%2d views   ',i)];
        end
        S = [S;'SfM: RANSAC'];
        pm = 1;
    end
    if gdata.showdetails
        S = [' 1 view    ';' 2 views   ';' 3 views   ';' 4 views   ';' Merged    ';' Classified'];
        pm = 1;
    end
    if pm

        g = get(handles.TrackPop,'Value');
        if g>size(S,1)
            set(handles.TrackPop,'Value',size(S,1));
            gdata.tracks = size(S,1);
            ShowImage(handles);
        end
        set(handles.TrackPop,'String',cellstr(S))
    end




end


% --- Executes on button press in Demo.
function Demo_Callback(hObject, eventdata, handles)
% hObject    handle to Demo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gdata


%gdata.img = gdata.f.imgmin;
for i=gdata.f.imgmin:gdata.f.imgmax
    gdata.img = i;
    ShowImage(handles)
    %NextImage_Callback(hObject, eventdata, handles)
    pause(0.5)
end

% if gdata.tok==0
%    Track_Callback(hObject, eventdata, handles)
% end
% gdata.img = gdata.f.imgmin;
% set(handles.All,'Value',0);
% set(handles.Mono,'Value',0);
% set(handles.Link,'Value',0);
% gdata.allobj = 0;
% gdata.monobj = 0;
% gdata.linkframe = 0;
%
% for i=gdata.f.imgmin:gdata.f.imgmax-1
%     NextImage_Callback(hObject, eventdata, handles)
%     pause(0.5)
% end
%
% set(handles.Mono,'Value',1);
% set(handles.Link,'Value',1);
% gdata.monobj = 1;
% gdata.linkframe = 1;
%
% gdata.nz = 1;
% gdata.img = gdata.f.imgmin;
% for i=gdata.f.imgmin:gdata.f.imgmax-1
%     NextImage_Callback(hObject, eventdata, handles)
%     pause(0.5)
% end
%
% gdata.img = gdata.f.imgmin;
% gdata.nz = 1;
% for i=1:gdata.nobj-1
%    NextObject_Callback(hObject, eventdata, handles)
%    pause(1)
% end
%
% set(handles.Mono,'Value',0);
% gdata.monobj = 0;
% set(handles.All,'Value',1);
% gdata.allobj = 1;
% ShowImage(handles)
%
% gdata.nz = 1;
% for i=gdata.f.imgmin:gdata.f.imgmax-1
%     NextImage_Callback(hObject, eventdata, handles)
%     pause(0.5)
% end
%
%


% --- Executes on selection change in TrackPop.
function TrackPop_Callback(hObject, eventdata, handles)
% hObject    handle to TrackPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns TrackPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TrackPop
global gdata
%contents = cellstr(get(hObject,'String'));
%gdata.tracks = str2num(contents{get(hObject,'Value')});
gdata.tracks = get(hObject,'Value');
ShowImage(handles)



% --- Executes during object creation, after setting all properties.
function TrackPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrackPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Multi.
function Multi_Callback(hObject, eventdata, handles)
% hObject    handle to Multi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Multi
global gdata
%contents = cellstr(get(hObject,'String'));
%gdata.tracks = str2num(contents{get(hObject,'Value')});
gdata.multilines = get(hObject,'Value');
ShowImage(handles)
