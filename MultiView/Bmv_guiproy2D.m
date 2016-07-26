function varargout = Bmv_guiproy2D(varargin)
% BMV_GUIPROY2D M-file for Bmv_guiproy2D.fig
%      BMV_GUIPROY2D, by itself, creates a new BMV_GUIPROY2D or raises the existing
%      singleton*.
%
%      H = BMV_GUIPROY2D returns the handle to a new BMV_GUIPROY2D or the handle to
%      the existing singleton*.
%
%      BMV_GUIPROY2D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BMV_GUIPROY2D.M with the given input arguments.
%
%      BMV_GUIPROY2D('Property','Value',...) creates a new BMV_GUIPROY2D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before matrixH_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Bmv_guiproy2D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help Bmv_guiproy2D

% Last Modified by GUIDE v2.5 23-Aug-2011 10:29:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Bmv_guiproy2D_OpeningFcn, ...
                   'gui_OutputFcn',  @Bmv_guiproy2D_OutputFcn, ...
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


% --- Executes just before Bmv_guiproy2D is made visible.
function Bmv_guiproy2D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Bmv_guiproy2D (see VARARGIN)

% Choose default command line output for Bmv_guiproy2D
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Bmv_guiproy2D wait for user response (see UIRESUME)
% uiwait(handles.figure1);
dx = 0;
dy = 0;
th = 0;
u0 = 0;
v0 = 0;
z  = 1;
x_p = [92   389   420    14]';
y_p = [57    71   305   278]';
save points dx dy th u0 v0 z x_p y_p





% --- Outputs from this function are returned to the command line.
function varargout = Bmv_guiproy2D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in left.
function left_Callback(hObject, eventdata, handles)
% hObject    handle to left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load points
dx = dx-[0 0 -20 20]';
x = x_p+dx;
y = y_p+dy;
save points u0 v0 z th dx dy x_p y_p
hshape(x,y)



% --- Executes on button press in right.
function right_Callback(hObject, eventdata, handles)
% hObject    handle to right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load points

dx = dx+[0 0 -20 20]';
x = x_p+dx;
y = y_p+dy;
save points u0 v0 z th dx dy x_p y_p
hshape(x,y)

% --- Executes on button press in A.
function A_Callback(hObject, eventdata, handles)
% hObject    handle to A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load points

dy = dy-[-20 0 0 20]';
x = x_p+dx;
y = y_p+dy;
save points u0 v0 z th dx dy x_p y_p
hshape(x,y)


% --- Executes on button press in V.
function V_Callback(hObject, eventdata, handles)
% hObject    handle to V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load points

dy = dy-[20 0 0 -20]';
x = x_p+dx;
y = y_p+dy;
save points u0 v0 z th dx dy x_p y_p
hshape(x,y)



% --- Executes on button press in points.
function points_Callback(hObject, eventdata, handles)
% hObject    handle to points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Select 4 points this way and press <Enter>:')
disp(' ');
disp('   1------------4');
disp('   |            |');
disp('   |            |');
disp('   |            |');
disp('   2------------3');
figure(1)
[y_p,x_p] = getpts;
dy = 0;
dx = 0;
th = 0;
z  = 1;
u0 = 0;
v0 = 0;
save points u0 v0 z th dx dy x_p y_p


% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img_name = input('File name: ');
I = imread(img_name);
[N,M,P] = size(I);
if P==3
    I = rgb2gray(I);
end
if (N>500)
    I = imresize(I,500/N);
end
figure(1)
imshow(I)
title('original image')
save I I



% --- Executes on button press in rotplus.
function rotplus_Callback(hObject, eventdata, handles)
% hObject    handle to rotplus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load points

th = th+pi/40;
x = x_p+dx;
y = y_p+dy;
save points u0 v0 z th dx dy x_p y_p
hshape(x,y)


function Hc = matrixHxy(x,y,x_p,y_p)
m1 = [x(1) x(2) x(3) x(4)
      y(1) y(2) y(3) y(4)
      1    1    1    1    ];

m2 = [x_p(1) x_p(2) x_p(3) x_p(4)
      y_p(1) y_p(2) y_p(3) y_p(4)
      1      1      1      1    ];
Hc = Bmv_homographySVD(m1,m2);



function hshape(x,y)
load I
load points
[N,M] = size(I);
% A1=[x(1) y(1) 1 0 0 0 -x(1)*x_p(1) -y(1)*x_p(1);
%     0 0 0 x(1) y(1) 1 -x(1)*y_p(1) -y(1)*y_p(1)];
% A2=[x(2) y(2) 1 0 0 0 -x(2)*x_p(2) -y(2)*x_p(2);
%     0 0 0 x(2) y(2) 1 -x(2)*y_p(2) -y(2)*y_p(2)];
% A3=[x(3) y(3) 1 0 0 0 -x(3)*x_p(3) -y(3)*x_p(3);
%     0 0 0 x(3) y(3) 1 -x(3)*y_p(3) -y(3)*y_p(3)];
% A4=[x(4) y(4) 1 0 0 0 -x(4)*x_p(4) -y(4)*x_p(4);
%     0 0 0 x(4) y(4) 1 -x(4)*y_p(4) -y(4)*y_p(4)];
% b=[x_p(1) y_p(1) x_p(2) y_p(2) x_p(3) y_p(3) x_p(4) y_p(4)];
% 
% h=inv([A1;A2;A3;A4])*b';
% 
% Hc =[h(1) h(2) h(3); h(4) h(5)  h(6);h(7) h(8) 1];

% disp('aqui')

% m1 = [x(1) x(2) x(3) x(4)
%       y(1) y(2) y(3) y(4)
%       1    1    1    1    ];
% 
% m2 = [x_p(1) x_p(2) x_p(3) x_p(4)
%       y_p(1) y_p(2) y_p(3) y_p(4)
%       1      1      1      1    ];
Hc = matrixHxy(x,y,x_p,y_p); 
%Hc = Bmv_homographySVD(m1,m2);

Ht = [1 0 -N/2;0 1 -M/2;0 0 1];
Hr  = [cos(th)      -sin(th)     0 
     sin(th)       cos(th)     0
     0             0           1];

Hz  = [z 0 u0;0 z v0; 0 0 1];


HH = Hz*Hc*inv(Ht)*Hr*Ht;





[N,M] = size(I);

J = Bmv_projective2D(I,HH,[N M],0);
save J J
figure(2)
imshow(uint8(J))
title('transformed image')

% --- Executes on button press in rotminus.
function rotminus_Callback(hObject, eventdata, handles)
% hObject    handle to rotminus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


load points

th = th-pi/40;
x = x_p+dx;
y = y_p+dy;
save points u0 v0 z th dx dy x_p y_p
hshape(x,y)


% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


img_format = input('Image format: ');
img_name   = input('Image name  : ');
load J
imwrite([img_name '.' img_format],img_format);


% --- Executes on button press in zoomplus.
function zoomplus_Callback(hObject, eventdata, handles)
% hObject    handle to zoomplus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


load points

z = z/1.1;
x = x_p+dx;
y = y_p+dy;
save points u0 v0 z th dx dy x_p y_p
hshape(x,y)


% --- Executes on button press in zoomminus.
function zoomminus_Callback(hObject, eventdata, handles)
% hObject    handle to zoomminus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


load points

z = z*1.1;
x = x_p+dx;
y = y_p+dy;
save points u0 v0 z th dx dy x_p y_p
hshape(x,y)


% --- Executes on button press in shiftleft.
function shiftleft_Callback(hObject, eventdata, handles)
% hObject    handle to shiftleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load points

v0 = v0+50*z;
x = x_p+dx;
y = y_p+dy;
save points u0 v0 z th dx dy x_p y_p
hshape(x,y)


% --- Executes on button press in shiftdown.
function shiftdown_Callback(hObject, eventdata, handles)
% hObject    handle to shiftdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load points

u0 = u0-50*z;
x = x_p+dx;
y = y_p+dy;
save points u0 v0 z th dx dy x_p y_p
hshape(x,y)


% --- Executes on button press in shiftright.
function shiftright_Callback(hObject, eventdata, handles)
% hObject    handle to shiftright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load points

v0 = v0-50*z;
x = x_p+dx;
y = y_p+dy;
save points u0 v0 z th dx dy x_p y_p
hshape(x,y)


% --- Executes on button press in shiftup.
function shiftup_Callback(hObject, eventdata, handles)
% hObject    handle to shiftup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load points

u0 = u0+50*z;
x = x_p+dx;
y = y_p+dy;
save points u0 v0 z th dx dy x_p y_p
hshape(x,y)




% --- Executes on button press in Ideal.
function Ideal_Callback(hObject, eventdata, handles)
% hObject    handle to Ideal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load I
load points
[N,M] = size(I);
x = [1 N N 1]';
y = [1 1 M M]';
% A1=[x(1) y(1) 1 0 0 0 -x(1)*x_p(1) -y(1)*x_p(1);
%     0 0 0 x(1) y(1) 1 -x(1)*y_p(1) -y(1)*y_p(1)];
% A2=[x(2) y(2) 1 0 0 0 -x(2)*x_p(2) -y(2)*x_p(2);
%     0 0 0 x(2) y(2) 1 -x(2)*y_p(2) -y(2)*y_p(2)];
% A3=[x(3) y(3) 1 0 0 0 -x(3)*x_p(3) -y(3)*x_p(3);
%     0 0 0 x(3) y(3) 1 -x(3)*y_p(3) -y(3)*y_p(3)];
% A4=[x(4) y(4) 1 0 0 0 -x(4)*x_p(4) -y(4)*x_p(4);
%     0 0 0 x(4) y(4) 1 -x(4)*y_p(4) -y(4)*y_p(4)];
% b=[x_p(1) y_p(1) x_p(2) y_p(2) x_p(3) y_p(3) x_p(4) y_p(4)];
% 
% h=inv([A1;A2;A3;A4])*b';
% 
% Hc =[h(1) h(2) h(3); h(4) h(5)  h(6);h(7) h(8) 1];

Hc = matrixHxy(x,y,x_p,y_p);


% [N,M] = size(I);

J = Bmv_projective2D(I,Hc,[N M],0);
save J J
figure(3)
imshow(uint8(J),[])
title('transformed image (ideal)')


% --- Executes on button press in Restore.
function Restore_Callback(hObject, eventdata, handles)
% hObject    handle to Restore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load points
dy = 0;
dx = 0;
th = 0;
z  = 1;
u0 = 0;
v0 = 0;
x = x_p;
y = y_p;
save points u0 v0 z th dx dy x_p y_p
hshape(x,y)


