% Bfx_gui
%
% Toolbox: Balu
%
%    Graphic User Interface for feature extraction.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl

function varargout = Bfx_gui(varargin)
% BFX_GUI M-file for Bfx_gui.fig
%      BFX_GUI, by itself, creates a new BFX_GUI or raises the existing
%      singleton*.
%
%      H = BFX_GUI returns the handle to a new BFX_GUI or the handle to
%      the existing singleton*.
%
%      BFX_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BFX_GUI.M with the given input arguments.
%
%      BFX_GUI('Property','Value',...) creates a new BFX_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Bfx_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Bfx_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Bfx_gui

% Last Modified by GUIDE v2.5 17-Nov-2011 10:21:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Bfx_gui_OpeningFcn, ...
    'gui_OutputFcn',  @Bfx_gui_OutputFcn, ...
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


% --- Executes just before Bfx_gui is made visible.
function Bfx_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Bfx_gui (see VARARGIN)

% Choose default command line output for Bfx_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Bfx_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
clc
disp('Balu 3.0: GUI for Feature Extraction.')
disp(' ')
disp('This GUI will extract features from current directory:')
cd
disp('Warning: not all features implemented in Balu are in this GUI.');
disp(' ')
disp('Please define the feature extraction process in the GUI window,');
disp('and press [Go] when you finish.')
disp(' ')

global Img Imgnow fextractor Nimg sfiles pseg mseg seg fsname

fextractor.formatimg     = 1;
fextractor.resizeimg     = 1;
fextractor.geostandard   = 0;
fextractor.geoinvariant  = 0;
fextractor.geofourierdes = 0;
fextractor.geoelliptical = 0;
fextractor.intstandard   = 0;
fextractor.intcontrast   = 0;
fextractor.intharalick   = 0;
fextractor.intfourierdct = 0;
fextractor.inthuint      = 0;
fextractor.intgabor      = 0;
fextractor.intlbp        = 0;
fextractor.inthog        = 0;
fextractor.colorgray     = 0;
fextractor.colorred      = 0;
fextractor.colorgreen    = 0;
fextractor.colorblue     = 0;
fextractor.colorhue      = 0;
fextractor.colorsat      = 0;
fextractor.colorval      = 0;
fextractor.colorl        = 0; %%%% Conversion RGB -> L*a*b*
fextractor.colora        = 0; % 0 = no, 1 = calibrated (Bim_rgb2lab), 
fextractor.colorb        = 0; % 2 = formulas (Bim_rgb2lab0)
fextractor.segmentation  = 0;
fextractor.partition     = 1;
fextractor.outmatlab     = 0;
fextractor.outascii      = 0;
fextractor.outexcel      = 0;
axis off

pseg = -0.05;
mseg = 0;
Nimg = 1;
seg  = 0;
fsname = 0;
sfiles = dir('*.jpg');
if not(isempty(sfiles))
    %    error('Balu Error: Current directory does not contain any image')
    %else
    Img = imread(sfiles(Nimg).name);
    Imgnow = Img;
    subplot(1,1,1); imshow(Img)
    title(sfiles(Nimg).name)

end


% --- Outputs from this function are returned to the command line.
function varargout = Bfx_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ButtonGo.
function ButtonGo_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonGo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fextractor
ok = 1;
fgeo = fextractor.geostandard+fextractor.geoinvariant+...
    fextractor.geofourierdes+fextractor.geoelliptical;

fint = fextractor.intstandard+fextractor.intcontrast+...
    fextractor.intharalick+fextractor.intfourierdct+...
    fextractor.inthuint+fextractor.intgabor+fextractor.intlbp+...
    fextractor.inthog;

fcol = fextractor.colorgray+fextractor.colorred+fextractor.colorgreen+...
    fextractor.colorblue+fextractor.colorhue+fextractor.colorsat+...
    fextractor.colorval+fextractor.colorl+fextractor.colora+...
    fextractor.colorb;
fout = fextractor.outmatlab+fextractor.outascii+fextractor.outexcel;


if fgeo+fint==0
    beep
    wwarndlg('You must set some features.','Error in Bfx_gui');
    ok = 0;
end
if and(fint>0,fcol==0)
    beep
    wwarndlg('If you select color features, you must set color component(s).','Error in Bfx_gui');
    ok = 0;
end
if and(fint==0,fcol>0)
    beep
    wwarndlg('If you select color component(s), you must set color features.','Error in Bfx_gui');
    ok = 0;
end
if (fout==0)
    beep
    wwarndlg('You must set output file format (Matlab, Text or Excel).','Error in Bfx_gui');
    ok = 0;
end

if ok
    yesnoans = questdlg('Bfx_gui will start to extract all features. This process could take several minutes. Are you sure?', ...
        'Bfx_gui Information', ...
        'Yes', 'No', 'No');
    if yesnoans(1)=='Y'
        save Bfx_guidata fextractor
        [f,fn] = Bfx_guifun(fextractor);
        questdlg(sprintf('Bfx_gui ended successfully: %d features extracted from %d images.',size(f,2),size(f,1)),'Bfx_gui Information','Ok', 'Ok');
    end
end
% --- Executes on button press in Gray.
function Gray_Callback(hObject, eventdata, handles)
% hObject    handle to Gray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Gray
global Imgnow fextractor
val = get(hObject,'Value');
if val
    if size(Imgnow,3)==1
        subplot(1,1,1); imshow(Imgnow)
    else
        subplot(1,1,1); imshow(rgb2gray(Imgnow))
    end
else
    subplot(1,1,1); imshow(Imgnow);
end
fextractor.colorgray = val;


% --- Executes on button press in Red.
function Red_Callback(hObject, eventdata, handles)
% hObject    handle to Red (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Red
global Imgnow fextractor
val = get(hObject,'Value');
if val
    subplot(1,1,1); imshow(uint8(Imgnow(:,:,1)),[(1:256)'/256 zeros(256,2)])
else
    subplot(1,1,1); imshow(Imgnow);
end
fextractor.colorred = val;


% --- Executes on button press in Green.
function Green_Callback(hObject, eventdata, handles)
% hObject    handle to Green (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Green
global Imgnow fextractor
val = get(hObject,'Value');
if val
    subplot(1,1,1); imshow(uint8(Imgnow(:,:,2)),[zeros(256,1) (1:256)'/256 zeros(256,1)])
else
    subplot(1,1,1); imshow(Imgnow);
end
fextractor.colorgreen = val;



% --- Executes on button press in Blue.
function Blue_Callback(hObject, eventdata, handles)
% hObject    handle to Blue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Blue
global Imgnow fextractor
val = get(hObject,'Value');
if val
    subplot(1,1,1); imshow(uint8(Imgnow(:,:,3)),[zeros(256,2) (1:256)'/256])
else
    subplot(1,1,1); imshow(Imgnow);
end
fextractor.colorblue = val;


% --- Executes on button press in Hue.
function Hue_Callback(hObject, eventdata, handles)
% hObject    handle to Hue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Hue
global Imgnow fextractor
val = get(hObject,'Value');
if val
    H = rgb2hsv(Imgnow);
    H(:,:,2) = 1;
    H(:,:,3) = 1;
    subplot(1,1,1); imshow(hsv2rgb(H))
else
    subplot(1,1,1); imshow(Imgnow);
end
fextractor.colorhue = val;


% --- Executes on button press in Sat.
function Sat_Callback(hObject, eventdata, handles)
% hObject    handle to Sat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Sat

global Imgnow fextractor
val = get(hObject,'Value');
if val
    H = rgb2hsv(Imgnow);
    subplot(1,1,1); imshow(H(:,:,2))
else
    subplot(1,1,1); imshow(Imgnow);
end
fextractor.colorsat = val;




% --- Executes on button press in Val.
function Val_Callback(hObject, eventdata, handles)
% hObject    handle to Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Val
global Imgnow fextractor
val = get(hObject,'Value');
if val
    H = rgb2hsv(Imgnow);
    subplot(1,1,1); imshow(H(:,:,3))
else
    subplot(1,1,1); imshow(Imgnow);
end
fextractor.colorval = val;



% --- Executes on button press in L.
function L_Callback(hObject, eventdata, handles)
% hObject    handle to L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of L
global Imgnow fextractor
val = get(hObject,'Value');
if val
    if (exist('LAB.mat','file'))
        load LAB
        L = Bim_rgb2lab(Imgnow,M);
        subplot(1,1,1); imshow(L(:,:,1),[])
    else
        beep        
        yesnoans = questdlg('Calibrated LAB conversion cannot be possible because LAB.mat file does not exist. Do you want to convert RGB to L*a*b* using CIE formulas?', ...
        'Bfx_gui RGB -> L*a*b conversion', ...
        'Yes', 'No', 'No');
        if yesnoans(1)=='Y'
           val = 2;
           L = Bim_rgb2lab0(Imgnow);
           subplot(1,1,1); imshow(L(:,:,1),[])
        else
           wwarndlg('L*a*b* conversion is not possible, because LAB.mat file does not exist.','Error in Bfx_gui');
           set(hObject,'Value',0);
           val = 0;
        end
    end
else
    subplot(1,1,1); imshow(Imgnow);
end
fextractor.colorl = val;





% --- Executes on button press in a.
function a_Callback(hObject, eventdata, handles)
% hObject    handle to a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of a
global Imgnow fextractor
val = get(hObject,'Value');
if val
    if (exist('LAB.mat','file'))
        load LAB
        L = rgb2lab(Imgnow,M);
        subplot(1,1,1); imshow(L(:,:,2),[])
    else
        beep        
        yesnoans = questdlg('Calibrated LAB conversion cannot be possible because LAB.mat file does not exist. Do you want to convert RGB to L*a*b* using CIE formulas?', ...
        'Bfx_gui RGB -> L*a*b conversion', ...
        'Yes', 'No', 'No');
        if yesnoans(1)=='Y'
           val = 2;
           L = Bim_rgb2lab0(Imgnow);
           subplot(1,1,1); imshow(L(:,:,2),[])
        else
           wwarndlg('L*a*b* conversion is not possible, because LAB.mat file does not exist.','Error in Bfx_gui');
           set(hObject,'Value',0);
           val = 0;
        end
    end
else
    subplot(1,1,1); imshow(Imgnow);
end
fextractor.colora = val;


% --- Executes on button press in b.
function b_Callback(hObject, eventdata, handles)
% hObject    handle to b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of b
global Imgnow fextractor
val = get(hObject,'Value');
if val
    if (exist('LAB.mat','file'))
        load LAB
        L = rgb2lab(Imgnow,M);
        subplot(1,1,1); imshow(L(:,:,3),[])
    else
        beep        
        yesnoans = questdlg('Calibrated LAB conversion cannot be possible because LAB.mat file does not exist. Do you want to convert RGB to L*a*b* using CIE formulas?', ...
        'Bfx_gui RGB -> L*a*b conversion', ...
        'Yes', 'No', 'No');
        if yesnoans(1)=='Y'
           val = 2;
           L = Bim_rgb2lab0(Imgnow);
           subplot(1,1,1); imshow(L(:,:,3),[])

        else
           wwarndlg('L*a*b* conversion is not possible, because LAB.mat file does not exist.','Error in Bfx_gui');
           set(hObject,'Value',0);
           val = 0;
        end
    end
else
    subplot(1,1,1); imshow(Imgnow);
end
fextractor.colorb = val;


% --- Executes on button press in Matlab.
function Matlab_Callback(hObject, eventdata, handles)
% hObject    handle to Matlab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Matlab
global fextractor
fextractor.outmatlab = get(hObject,'Value');


% --- Executes on button press in Ascii.
function Ascii_Callback(hObject, eventdata, handles)
% hObject    handle to Ascii (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Ascii
global fextractor
fextractor.outascii = get(hObject,'Value');


% --- Executes on button press in Excel.
function Excel_Callback(hObject, eventdata, handles)
% hObject    handle to Excel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Excel
global fextractor
fextractor.outexcel = get(hObject,'Value');


% --- Executes on button press in GeoStandard.
function GeoStandard_Callback(hObject, eventdata, handles)
% hObject    handle to GeoStandard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GeoStandard
global fextractor
fextractor.geostandard = get(hObject,'Value');
if (fextractor.segmentation == 0) && (fextractor.geostandard == 1)
    wwarndlg('Geometrical features with no segmentation? Bad idea :(','Warning in Bfx_gui');
end


% --- Executes on button press in IntMoments.
function IntMoments_Callback(hObject, eventdata, handles)
% hObject    handle to IntMoments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IntMoments
global fextractor
fextractor.geoinvariant = get(hObject,'Value');
if (fextractor.segmentation == 0) && (fextractor.geoinvariant == 1)
    wwarndlg('Geometrical features with no segmentation? Bad idea :(','Warning in Bfx_gui');
end


% --- Executes on button press in FourierDesc.
function FourierDesc_Callback(hObject, eventdata, handles)
% hObject    handle to FourierDesc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FourierDesc
global fextractor
fextractor.geofourierdes = get(hObject,'Value');
if (fextractor.segmentation == 0) && (fextractor.geofourierdes == 1)
    wwarndlg('Geometrical features with no segmentation? Bad idea :(','Warning in Bfx_gui');
end


% --- Executes on button press in Elliptical.
function Elliptical_Callback(hObject, eventdata, handles)
% hObject    handle to Elliptical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Elliptical
global fextractor
fextractor.geoelliptical = get(hObject,'Value');
if (fextractor.segmentation == 0) && (fextractor.geoelliptical == 1)
    wwarndlg('Geometrical features with no segmentation? Bad idea :(','Warning in Bfx_gui');
end


% --- Executes on button press in IntStandard.
function IntStandard_Callback(hObject, eventdata, handles)
% hObject    handle to IntStandard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IntStandard
global fextractor
fextractor.intstandard = get(hObject,'Value');


% --- Executes on button press in IntContrast.
function IntContrast_Callback(hObject, eventdata, handles)
% hObject    handle to IntContrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IntContrast
global fextractor
fextractor.intcontrast = get(hObject,'Value');


% --- Executes on button press in IntHaralick.
function IntHaralick_Callback(hObject, eventdata, handles)
% hObject    handle to IntHaralick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IntHaralick
global fextractor
fextractor.intharalick = get(hObject,'Value');


% --- Executes on button press in IntFourierDCT.
function IntFourierDCT_Callback(hObject, eventdata, handles)
% hObject    handle to IntFourierDCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IntFourierDCT
global fextractor
fextractor.intfourierdct = get(hObject,'Value');


% --- Executes on button press in IntHu.
function IntHu_Callback(hObject, eventdata, handles)
% hObject    handle to IntHu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IntHu
global fextractor
fextractor.inthuint = get(hObject,'Value');


% --- Executes on button press in IntGabor.
function IntGabor_Callback(hObject, eventdata, handles)
% hObject    handle to IntGabor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IntGabor
global fextractor
fextractor.intgabor = get(hObject,'Value');


% --- Executes on selection change in FormatImg.
function FormatImg_Callback(hObject, eventdata, handles)
% hObject    handle to FormatImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns FormatImg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FormatImg
global Img Imgnow fextractor Nimg sfiles
fextractor.formatimg = get(hObject,'Value');
s = dir('*.jpg');
if not(isempty(s))
    %    error('Balu Error: Current directory does not contain any image')
    %else
end
switch fextractor.formatimg
    case 2
        s = 'tif';
    case 3
        s = 'bmp';
    case 4
        s = 'png';
    case 5
        s = 'ppm';
    case 6
        s = 'gif';
    case 7
        s = 'pbm';
    otherwise
        s = 'jpg';
end
if isempty(s)
    beep
    wwarndlg('Current directory does not contain any image with this format.','Error in Bfx_gui');
else
    sfiles = [dir(['*.' lower(s)]); dir(['*.' upper(s)])];
    Nimg = 1;
    Img = imread(sfiles(1).name);
    Imgnow = Img;
    subplot(1,1,1); imshow(Img)
    title(sfiles(Nimg).name)
end

% --- Executes during object creation, after setting all properties.
function FormatImg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FormatImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function ImgResize_Callback(hObject, eventdata, handles)
% hObject    handle to ImgResize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImgResize as text
%        str2double(get(hObject,'String')) returns contents of ImgResize as a double
global fextractor Img Imgnow
fextractor.resizeimg = str2num(get(hObject,'String'));
warning off
Imgnow = imresize(Img,fextractor.resizeimg);
warning on
subplot(1,1,1); imshow(Imgnow)

% --- Executes during object creation, after setting all properties.
function ImgResize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgResize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ImgSegmentation.
function ImgSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to ImgSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ImgSegmentation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ImgSegmentation
global fextractor Imgnow pseg  seg mseg fsname
fextractor.segmentation = get(hObject,'Value')-1;
seg = fextractor.segmentation;
fsname = 0;
PlotImage(Imgnow,seg,pseg,mseg);

function [R,E,J] = SegImage(Imgnow,seg,pseg,mseg)
global fsname
[N,M,P] = size(Imgnow);
switch seg
    case 1 % balu
        [R,E,J] = Bim_segbalu(Imgnow,pseg);
    case 2 % maxvar
        [R,E,J] = Bim_segmaxvar(Imgnow,pseg);
    case 3 % maxfisher
        [R,E,J] = Bim_segmaxfisher(Imgnow,pseg);
    case 4 % otsu
        [R,E,J] = Bim_segotsu(Imgnow,pseg);
    case 5
        [R,E,J] = Bim_segpca(Imgnow,pseg);
    case 6
        if fsname == 0
        sname = input('Input segmentation program name: ');
        fsname = ['[R,E,J] = ' sname '(Imgnow,pseg);'];
        end
        eval(fsname);
    otherwise
        R = ones(N,M);
        if (P==3)
            J = rgb2gray(Imgnow);
        else
            J = Imgnow;
        end
        R = ones(size(J));
        E = zeros(size(J));
end
if mseg
    R = bwlabel(R);
end


function PlotImage(Imgnow,seg,pseg,mseg)
[R,E,J] = SegImage(Imgnow,seg,pseg,mseg);
[N,M,P] = size(Imgnow);
subplot(1,1,1); imshow(R,[])
pause(1)
ii = find(R==0);
if isempty('RR')
    Iw = Imgnow;
else
    if (P==3)
        Ir = Imgnow(:,:,1);Ir(ii)=0;
        Ig = Imgnow(:,:,2);Ig(ii)=0;
        Ib = Imgnow(:,:,3);Ib(ii)=0;
        Iw = zeros(size(Imgnow));
        Iw(:,:,1) = Ir;
        Iw(:,:,2) = Ig;
        Iw(:,:,3) = Ib;
    else
        Iw = Imgnow;
        Iw(ii) = 0;
    end
    subplot(1,1,1); imshow(uint8(Iw))
end


% --- Executes during object creation, after setting all properties.
function ImgSegmentation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in lbp.
function lbp_Callback(hObject, eventdata, handles)
% hObject    handle to lbp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lbp
global fextractor
fextractor.intlbp = get(hObject,'Value');

% --- Executes on selection change in Partition.
function Partition_Callback(hObject, eventdata, handles)
% hObject    handle to Partition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Partition contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Partition
global fextractor Imgnow
fextractor.partition = get(hObject,'Value');
[N,M,P] = size(Imgnow);
n = N/fextractor.partition;
m = M/fextractor.partition;
subplot(1,1,1); imshow(Imgnow)
hold on
for i=0:n:N
    plot([1 M],[i i])
end
for j=0:m:M
    plot([j j],[1 N])
end
hold off


% --- Executes during object creation, after setting all properties.
function Partition_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Partition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% Bfx_guifun
%
% Toolbox: Balu
%    Extract feactures of all images of courrent directory and store.
%    The extracted features are stored in excel, text or matlab files.
%    A list of the names of the images is stored in imgnames.txt
%    f = table of features, fn = name of features
%
% fetractor is a structure with following variables:
%
% fextractor.formatimg         1=jpg|2=tif|3=bmp|4=png|5=ppm|6=png|7=pbm
% fextractor.resizeimg         0=no|0.5, 0.25, [16 16], [32 32], [200 200]
% fextractor.geostandard       1=yes|0=no
% fextractor.geoinvariant           :
% fextractor.geofourierdes
% fextractor.geoelliptical
% fextractor.intstandard
% fextractor.intcontrast            :
% fextractor.intharalick
% fextractor.intfourierdct
% fextractor.inthuint
% fextractor.intgabor
% fextractor.intlbp
% fextractor.colorgray
% fextractor.colorred               :
% fextractor.colorgreen
% fextractor.colorblue
% fextractor.colorhue
% fextractor.colorsat
% fextractor.colorval               :
% fextractor.colorl
% fextractor.colora
% fextractor.colorb
% fextractor.segmentation  1=yes | 0 = no (whole image)
% fextractor.outmatlab
% fextractor.outascii
% fextractor.outexcel
%
% D.Mery, PUC-DCC, Apr. 2008
% http://dmery.ing.puc.cl
%
function [f,fn]=Bfx_guifun(fextractor)
fsb = Bio_statusbar('Extracting features');
global pseg mseg seg

disp('Balu Full Features Extractor');
disp(' ');
reset(gca);
imgfmt = fextractor.formatimg;

fg = zeros(10,1);
fi = zeros(10,1);
co = zeros(20,1);

fg(1) = fextractor.geostandard;
fg(2) = fextractor.geoinvariant;
fg(3) = fextractor.geofourierdes;
fg(4) = fextractor.geoelliptical;
rg = sum(fg)>0;

fi(1) = fextractor.intstandard;
fi(2) = fextractor.intcontrast;
fi(3) = fextractor.intharalick;
fi(4) = fextractor.intfourierdct;
fi(5) = fextractor.inthuint;
fi(6) = fextractor.intgabor;
fi(7) = fextractor.intlbp;
fi(8) = fextractor.inthog;
ri = sum(fi)>0;

co(1) = fextractor.colorgray;
co(2) = fextractor.colorred;
co(3) = fextractor.colorgreen;
co(4) = fextractor.colorblue;
co(5) = fextractor.colorhue;
co(6) = fextractor.colorsat;
co(7) = fextractor.colorval;
co(8) = fextractor.colorl;
co(9) = fextractor.colora;
co(10)= fextractor.colorb;
if and(sum(co)>0,sum(fi)==0)
    disp('Warning: No intensity feature will be extracted, because no intensity feature was marked. ')
    ri = 0;
end
if and(sum(co)==0,sum(fi)>0)
    disp('Warning: No intensity feature will be extracted, because no color component was marked. ')
    ri = 0;
end


if (sum(co(8:10)>0)>0)
    rgblabconv = (sum(co(8:10)==2)>0)+1;
else
    rgblabconv = 0;
end


if (sum([fg;fi])==0)
    error('No feature will be extracted, because there is no feature marked. ')
end
imgseg  = fextractor.segmentation;
imgpart = fextractor.partition;
if imgpart == 0
    imgpart = 1; % no partition means one partition!
end
oxls    = fextractor.outexcel;
otxt    = fextractor.outascii;
omat    = fextractor.outmatlab;
imgres  = fextractor.resizeimg;

if rg
    rgi = 0;
    % fg(1) = fextractor.geostandard;
    if fg(1)
        rgi = rgi+1;
        bg(rgi).name = 'basicgeo';
        bg(rgi).options.show = 0;
    end

    % fg(2) = fextractor.geoinvariant;
    if fg(2)
        rgi = rgi+1;
        bg(rgi).name = 'hugeo';
        bg(rgi).options.show = 0;
        rgi = rgi+1;
        bg(rgi).name = 'flusser';
        bg(rgi).options.show = 0;
    end

    % fg(3) = fextractor.geofourierdes;
    if fg(3)
        rgi = rgi+1;
        bg(rgi).name = 'fourierdes';
        bg(rgi).options.show = 0;
        bg(rgi).options.Nfourierdes = 8;
    end

    % fg(4) = fextractor.geoelliptical;
    if fg(4)
        rgi = rgi+1;
        bg(rgi).name = 'fitellipse';
        bg(rgi).options.show = 0;
    end

    opg.b = bg;

end

if ri
    rii = 0;
    % fi(1) = fextractor.intstandard;
    if fi(1)
        rii = rii+1;
        bi(rii).name = 'basicint';
        bi(rii).options.show = 0;
    end

    % fi(2) = fextractor.intcontrast;
    if fi(2)
        rii = rii+1;
        bi(rii).name = 'contrast';
        bi(rii).options.neighbor = 1;
        bi(rii).options.param = 1.5;
        bi(rii).options.show = 0;
    end

    % fi(3) = fextractor.intharalick;
    if fi(3)
        rii = rii+1;
        bi(rii).name = 'haralick';
        bi(rii).options.dharalick = 1;
        bi(rii).options.show = 0;
        rii = rii+1;
        bi(rii).name = 'haralick';
        bi(rii).options.dharalick = 2;
        bi(rii).options.show = 0;
        rii = rii+1;
        bi(rii).name = 'haralick';
        bi(rii).options.dharalick = 3;
        bi(rii).options.show = 0;
        rii = rii+1;
        bi(rii).name = 'haralick';
        bi(rii).options.dharalick = 4;
        bi(rii).options.show = 0;
        rii = rii+1;
        bi(rii).name = 'haralick';
        bi(rii).options.dharalick = 5;
        bi(rii).options.show = 0;
    end

    % fi(4) = fextractor.intfourierdct;
    if fi(4)
        rii = rii+1;
        bi(rii).name = 'fourier';
        bi(rii).options.Nfourier = 64;
        bi(rii).options.Mfourier = 64;
        bi(rii).options.nfourier = 4;
        bi(rii).options.mfourier = 4;
        bi(rii).options.show = 0;
        rii = rii+1;
        bi(rii).name = 'dct';
        bi(rii).options.Ndct = 64;
        bi(rii).options.Mdct = 64;
        bi(rii).options.ndct = 4;
        bi(rii).options.mdct = 4;
        bi(rii).options.show = 0;
    end

    % fi(5) = fextractor.inthuint;
    if fi(5)
        rii = rii+1;
        bi(rii).name = 'huint';
        bi(rii).options.show = 0;
    end

    % fi(6) = fextractor.intgabor;
    if fi(6)
        rii = rii+1;
        bi(rii).name = 'gabor';
        bi(rii).options.Lgabor = 8;
        bi(rii).options.Sgabor = 8;
        bi(rii).options.fhgabor = 2;
        bi(rii).options.flgabor = 0.1;
        bi(rii).options.Mgabor = 21;
        bi(rii).options.show = 0;
    end

    % fi(7) = fextractor.intlbp;
    if fi(7)
        rii = rii+1;
        bi(rii).name = 'lbp';
        bi(rii).options.vdiv = 1;
        bi(rii).options.hdiv = 1;
        bi(rii).options.semantic = 0;
        bi(rii).options.samples = 8;
        bi(rii).options.mappingtype = 'u2';
        bi(rii).options.show = 0;
    end
    if fi(8)
        rii = rii+1;
        bi(rii).name = 'hog';
        bi(rii).options.ni = 1;
        bi(rii).options.nj = 1;
        bi(rii).options.B = 9;
        bi(rii).options.show = 0;
    end
    opi.b = bi;
end

warning off
if ispc
    switch imgfmt
        case 2
            !dir *.tif/B > imgnames.txt
            sdir = dir('*.tif');
        case 3
            !dir *.bmp/B > imgnames.txt
            sdir = dir('*.bmp');
        case 4
            !dir *.png/B > imgnames.txt
            sdir = dir('*.png');
        case 5
            !dir *.ppm/B > imgnames.txt
            sdir = dir('*.ppm');
        case 6
            !dir *.gif/B > imgnames.txt
            sdir = dir('*.gif');
        case 7
            !dir *.pbm/B > imgnames.txt
            sdir = dir('*.pbm');

        otherwise
            !dir *.jpg/B > imgnames.txt
            sdir = dir('*.jpg');
    end
else % Asuming unix, ie isunix should be 1
    switch imgfmt
        case 2
            !ls *.tif > imgnames.txt
            sdir = dir('*.tif');
        case 3
            !ls *.bmp > imgnames.txt
            sdir = dir('*.bmp');
        case 4
            !ls *.png > imgnames.txt
            sdir = dir('*.png');
        case 5
            !ls *.ppm > imgnames.txt
            sdir = dir('*.ppm');
        case 6
            !ls *.png > imgnames.txt
            sdir = dir('*.png');
        case 7
            !ls *.pbm > imgnames.txt
            sdir = dir('*.pbm');
        otherwise
            !ls *.jpg > imgnames.txt
            sdir = dir('*.jpg');
    end
end
fid = fopen('imgnames.txt','rt');
sdirn = length(sdir);
ok = 1;
if rg
    FeatureGeo = [];
end

if ri
    neigh = 0;
    %costr = ['Gray '; 'Red  '; 'Green'; 'Blue '; 'Hue  '; 'Sat  '; 'Value'; 'L    '; 'a    '; 'b    ' ];
    costr = 'gRGBHSVLab';
    color_str = costr(co==1);
    %for i=1:10
    %    if co(i)
    %        s = sprintf('Feature%s = [];',costr(i,:));
    %        eval(s);
    %    end
    %end
end
fall  = []; % extracted features
falln = []; % feature names
imgnames = [];
sdiri = 0;

while(ok)
    fsb = Bio_statusbar(sdiri/sdirn,fsb);
    sdiri = sdiri+1;
    s = fscanf(fid,'%s\n',[1 1]);
    if isempty(s)
        ok = 0;
    else
        disp(sprintf('\n\nProcessing image %s...',s));
        Io = imread(s);
        ss = [s '                                      '];
        sname = ss(1:32);
        %[N,M,P] = size(Io);
        Io = double(Io);
        if imgres(1)>0
            %            k = 2^(imgres+3);
            %            I = imresize(Io,[k k]);
            mini = min(Io(:));
            maxi = max(Io(:));
            I = Bim_sat(imresize(Io,imgres),mini,maxi);
        else
            I = Io;
        end
        [N,M,P] = size(I);

        % segmentation yes/no
        [R,E,J] = SegImage(I,seg,pseg,mseg);
        if sum2(R)==0
            if P==3 % 3 channels
                subplot(2,2,1);imshow(I/256,[]); title(s)
            else
                subplot(2,2,1);imshow(Io,[]); title(s)
            end
            subplot(2,2,2);imshow(R,[]);title('Segmented Region')

            str = ['Image ' s ' has no segmentation. The feature extraction will fail. Do you want to segment the whole image?'];
            yesnoans = questdlg(str, ...
                'Bfx_gui: Segmentation Failure', ...
                'Yes', 'Abort', 'Abort');
            if yesnoans(1)=='Y'
                R = ones(N,M);
            end
        end


        %        switch imgseg
        %            case 1
        %                [R,E,J] = Bim_segbalu(I);
        %            case 2
        %                [R,E,J] = Bim_segbalu(I,0.20);
        %            case 3
        %                [R,E,J] = Bim_segbalu(I,-0.20);
        %            case 4
        %                [R,E,J] = Bim_segbalu(I);
        %                R = bwlabel(R);
        %            case 5
        %                [R,E,J] = Bim_segbalu(I,0.20);
        %                R = bwlabel(R);
        %            case 6
        %                [R,E,J] = Bim_segbalu(I,-0.20);
        %                R = bwlabel(R);
        %            otherwise
        %                R = ones(N,M);
        %                if (P==3)
        %                    J = rgb2gray(I);
        %                else
        %                    J = I;
        %                end
        %        end

        % I: imagen a procesar (color or grayvalues)
        % J: imagen en gray values (transformed or original)
        % R: segmented image (if no R = ones(N,M)


        % color vs. grayvalues

        if (P==3)
            Imgnow = uint8(round(I));
        else
            Imgnow = I;
        end
        ii = find(R==0);
        if P==3 % 3 channels
            subplot(2,2,1);imshow(I/256,[]); title(s)
            Ir = Imgnow(:,:,1);
            Ig = Imgnow(:,:,2);
            Ib = Imgnow(:,:,3);
        else
            subplot(2,2,1);imshow(Io,[]); title(s)
            Iw = Imgnow;
        end
        if not(isempty(ii))
            if (P==3)
                Ir(ii)=0;
                Ig(ii)=0;
                Ib(ii)=0;
            else
                Iw(ii) = 0;
            end
        end
        if P==3
            Iw = zeros(size(Imgnow));
            Iw(:,:,1) = Ir;
            Iw(:,:,2) = Ig;
            Iw(:,:,3) = Ib;
            imshow(uint8(Iw))
        else
            imshow(Iw,[])
        end
        title(s)
        subplot(2,2,2);imshow(R,[]);title('Segmented Region')
        pause(0.1)
        RRR = R;
        III = I;
        dpi = fix(N/imgpart);
        dpj = fix(M/imgpart);
        fpart = [];
        fnpart = [];
        fnimg = []; % feature names
        fimg = []; % fetaures for this image
        for parti=1:imgpart
            for partj=1:imgpart
                if imgpart > 1
                    spart = [num2str(parti) num2str(partj)];
                else
                    spart = '  ';
                end
                fn = []; % feature names
                fs = []; % extracted features for this partition
                R = RRR(dpi*(parti-1)+1:dpi*parti,dpj*(partj-1)+1:dpj*partj);
                I = III(dpi*(parti-1)+1:dpi*parti,dpj*(partj-1)+1:dpj*partj,:);
                %figure(10)
                %imshow(I/256,[])
                %figure(11)
                %imshow(III/256,[])
                if rg
                    %[NameGeo,Feature,UnitGeo] = geofeatures(R,fg);
                    [Feature,NameGeo] = Bfx_geo(R,opg);
                    %FeatureGeo = [FeatureGeo; Feature];
                    fn = [fn; NameGeo];
                    fs = [fs Feature];
                end
                if ri

                    if sum(co(2:4))>0
                        XRGB = I;
                    end
                    if sum(co(5:7))>0
                        XHSV = rgb2hsv(I);
                    end
                    if rgblabconv>0
                        if rgblabconv==1
                        load LAB
                        XLAB = Bim_rgb2lab(I,M);
                        else
                        XLAB = Bim_rgb2lab0(I);
                        end
                    end
                    color_str = '';
                    for i=1:10
                        if co(i)
                            switch i
                                case 1 % Gray
                                    if size(I,3)==3
                                        X = rgb2gray(I/256)*256;
                                        Xh = uint8(round(X));
                                    else
                                        X = double(I);
                                        if max(X(:))>256
                                            Xh = uint8(double(X)/256);
                                            X = X/256;
                                        else
                                            Xh = uint8(X);
                                        end
                                    end
                                case 2 % Red
                                    X = XRGB(:,:,1);
                                    Xh = uint8(X);
                                case 3 % Green
                                    X = XRGB(:,:,2);
                                    Xh = uint8(X);
                                case 4 % Blue
                                    X = XRGB(:,:,3);
                                    Xh = uint8(X);
                                case 5 % H
                                    X = XHSV(:,:,1);
                                    Xh = double(X);
                                case 6 % S
                                    X = XHSV(:,:,2);
                                    Xh = double(X);
                                case 7 % V
                                    X = XHSV(:,:,3);
                                    Xh = uint8(X);
                                case 8 % L*
                                    X = XLAB(:,:,1);
                                    Xh = double(X)/100;
                                case 9 % a*
                                    X = XLAB(:,:,2);
                                    Xh = (double(X)+120)/240;
                                case 10 % b*
                                    X = XLAB(:,:,3);
                                    Xh = (double(X)+120)/240;
                            end
                            opi.colstr = costr(i);

                            [Feature,FeatureName] = Bfx_int(X,R,opi);
                            fn = [fn; FeatureName];
                            fs = [fs Feature];

                            subplot(2,2,3);imshow(X,[]);title(['Channel-' costr(i)])
                            %if (P==3)
                            %    subplot(2,2,4);imhist(X);title([costr(i,:) ' Histogram'])
                            %else
                            subplot(2,2,4);imhist(Xh);title(['Histogram-' costr(i)])
                            %end

                            pause(0)
                        end
                    end
                end
                fn = [fn ones(size(fn,1),1)*spart];
                fnimg = [fnimg;fn]; % feature names
                fimg = [fimg fs]; % fetaures for this image
                for si = 1:size(fs,1)
                    imgnames = [imgnames; sname];
                end

            end
        end
        fall = [fall; fimg];
    end
    fnall = fnimg;
    save % ..backup
end
f = fall;
fn = fnall;
fclose(fid);
if omat
    disp('Bfx_gui: saving f (features) and fn (feature names) to Bfx_results.mat...')
    save Bfx_results f fn imgnames
end
if oxls
    % balu2xls
    disp('Bfx_gui: saving features to Bfx_results.csv (excel compatible)...')
    csvwrite('Bfx_results.csv',f);
    disp('Saving feature names in Bfx_guinames.txt...')
    save Bfx_imgnames.txt fn -ascii
end
if otxt
    disp('Bfx_gui: saving features to Bfx_results.txt...')
    save Bfx_results.txt f -ascii -double
    disp('Saving feature names in Bfx_guinames.txt...')
    save Bfx_imgnames.txt fn -ascii
end
warning on
disp(' ')
disp(sprintf('Bfx_gui ended successfully: %d features extracted from %d images.',size(f,2),size(f,1)))
disp('(all variables saved in matlab.mat as backup)')
delete(fsb);



% --- Executes on button press in PreviousImage.
function PreviousImage_Callback(hObject, eventdata, handles)
% hObject    handle to PreviousImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Nimg sfiles Imgnow
nf = length(sfiles);
if nf>0
    Nimg = Nimg-1;
    if Nimg<1
        Nimg = nf;
    end
    Img = imread(sfiles(Nimg).name);
    Imgnow = Img;
    subplot(1,1,1); imshow(Img)
    title(sfiles(Nimg).name)
end

% --- Executes on button press in NextImage.
function NextImage_Callback(hObject, eventdata, handles)
% hObject    handle to NextImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Nimg sfiles Imgnow
nf = length(sfiles);
if nf>0
    Nimg = Nimg+1;
    if Nimg>length(sfiles)
        Nimg = 1;
    end
    Img = imread(sfiles(Nimg).name);
    Imgnow = Img;
    subplot(1,1,1); imshow(Img)
    title(sfiles(Nimg).name)
end


% --- Executes on button press in MultipleSegmentation.
function MultipleSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to MultipleSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MultipleSegmentation
global mseg
mseg = 1-mseg;


% --- Executes on button press in MinusSegmentation.
function MinusSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to MinusSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pseg mseg seg Imgnow
pseg = pseg+0.05;
PlotImage(Imgnow,seg,pseg,mseg)

% --- Executes on button press in PlusSegmentation.
function PlusSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to PlusSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pseg mseg seg Imgnow
pseg = pseg-0.05;
PlotImage(Imgnow,seg,pseg,mseg)


% --- Executes on button press in DoSegmentation.
function DoSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to DoSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Imgnow pseg  seg mseg
PlotImage(Imgnow,seg,pseg,mseg);


% --- Executes on button press in IntHOG.
function IntHOG_Callback(hObject, eventdata, handles)
% hObject    handle to IntHOG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IntHOG

global fextractor
fextractor.inthog = get(hObject,'Value');

