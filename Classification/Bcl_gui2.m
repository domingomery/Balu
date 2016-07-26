% Bcl_gui2
%
% Toolbox: Balu
%
%    Graphic User Interface for feature extraction.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl

function varargout = Bcl_gui2(varargin)
% BCL_GUI2 M-file for Bcl_gui2.fig
%      BCL_GUI2, by itself, creates a new BCL_GUI2 or raises the existing
%      singleton*.
%
%      H = BCL_GUI2 returns the handle to a new BCL_GUI2 or the handle to
%      the existing singleton*.
%
%      BCL_GUI2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BCL_GUI2.M with the given input arguments.
%
%      BCL_GUI2('Property','Value',...) creates a new BCL_GUI2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Bcl_gui2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Bcl_gui2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Bcl_gui2

% Last Modified by GUIDE v2.5 06-Mar-2015 18:11:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Bcl_gui2_OpeningFcn, ...
    'gui_OutputFcn',  @Bcl_gui2_OutputFcn, ...
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


% --- Executes just before Bcl_gui2 is made visible.
function Bcl_gui2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Bcl_gui2 (see VARARGIN)

% Choose default command line output for Bcl_gui2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Bcl_gui2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

clc
disp('Balu 3.0: GUI for feature and classifier selection.')
disp(' ')
disp('This GUI will load data stored in current directory:')
cd
disp('Warning: not all features selection and classifier algorithm implemented');
disp('in Balu are in this GUI.');
disp(' ')
disp('Please define the feature selection and classifier processes in the GUI window,');
disp('and press [Go] when you finish.')
disp(' ')
xlabel('Number of features')
ylabel('Performance [%]');
s.dmin            = 0;
s.maha            = 0;
s.lda             = 0;
s.qda             = 0;
s.knn4            = 0;
s.knn6            = 0;
s.knn8            = 0;
s.knn10           = 0;
s.knn15           = 0;
s.knn20           = 0;
s.libsvma         = 0;
s.svma            = 0;
s.svmb            = 0;
s.svmc            = 0;
s.svmd            = 0;
s.svme            = 0;
s.pnn             = 0;
s.nnglma          = 0;
s.nnglmb          = 0;
s.nnglmc          = 0;
s.adaboost5       = 0;
s.adaboost10      = 0;
s.boosting025     = 0;
s.boosting050     = 0;
s.sfsfisher       = 0;
s.fosmod          = 0;
s.lsef            = 0;
s.frank           = 0;
s.sfsknn5         = 0;
s.sfslda          = 0;
s.sfsqda          = 0;
s.sfssvmlin       = 0;
s.sfslibsvmlin    = 0;
s.sfsglma         = 0;
s.fsrandom        = 0;
s.allfeatures     = 0;
s.mmax            = 20;
s.crossvalfolders = 10;
s.filename        = '';
s.advanced        = 0;
axis([0 s.mmax 50 100])

save Bcl_guidata s


% --- Outputs from this function are returned to the command line.
function varargout = Bcl_gui2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in dmin.
function dmin_Callback(hObject, eventdata, handles)
% hObject    handle to dmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dmin
load Bcl_guidata
s.dmin = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Minimal Distance Classifier = %s',offon(s.dmin+1,:)))


% --- Executes on button press in maha.
function maha_Callback(hObject, eventdata, handles)
% hObject    handle to maha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of maha
load Bcl_guidata
s.maha = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Mahalanobis Distance Classifier = %s',offon(s.maha+1,:)))


% --- Executes on button press in lda.
function lda_Callback(hObject, eventdata, handles)
% hObject    handle to lda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lda
load Bcl_guidata
s.lda = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Linear Discriminant Analysis Classifier = %s',offon(s.lda+1,:)))


% --- Executes on button press in qda.
function qda_Callback(hObject, eventdata, handles)
% hObject    handle to qda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of qda
load Bcl_guidata
s.qda = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Quadratic Discriminant Analysis Classifier = %s',offon(s.qda+1,:)))


% --- Executes on button press in knn4.
function knn4_Callback(hObject, eventdata, handles)
% hObject    handle to knn4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of knn4
load Bcl_guidata
s.knn4 = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Nearest Neighbours Classifier with 4 neighbours = %s',offon(s.knn4+1,:)))

% --- Executes on button press in knn6.
function knn6_Callback(hObject, eventdata, handles)
% hObject    handle to knn6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of knn6
load Bcl_guidata
s.knn6 = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Nearest Neighbours Classifier with 6 neighbours = %s',offon(s.knn6+1,:)))


% --- Executes on button press in knn8.
function knn8_Callback(hObject, eventdata, handles)
% hObject    handle to knn8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of knn8
load Bcl_guidata
s.knn8 = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Nearest Neighbours Classifier with 8 neighbours = %s',offon(s.knn8+1,:)))


% --- Executes on button press in knn10.
function knn10_Callback(hObject, eventdata, handles)
% hObject    handle to knn10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of knn10
load Bcl_guidata
s.knn10 = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Nearest Neighbours Classifier with 10 neighbours = %s',offon(s.knn10+1,:)))


% --- Executes on button press in knn15.
function knn15_Callback(hObject, eventdata, handles)
% hObject    handle to knn15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of knn15
load Bcl_guidata
s.knn15 = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Nearest Neighbours Classifier with 15 neighbours = %s',offon(s.knn15+1,:)))


% --- Executes on button press in knn20.
function knn20_Callback(hObject, eventdata, handles)
% hObject    handle to knn20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of knn20
load Bcl_guidata
s.knn20 = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Nearest Neighbours Classifier with 20 neighbours = %s',offon(s.knn20+1,:)))


% --- Executes on button press in svma.
function svma_Callback(hObject, eventdata, handles)
% hObject    handle to svma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of svma
load Bcl_guidata
s.svma = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Support Vector Machine with Linear kernel or dot product = %s',offon(s.svma+1,:)))



% --- Executes on button press in svmb.
function svmb_Callback(hObject, eventdata, handles)
% hObject    handle to svmb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of svmb
load Bcl_guidata
s.svmb = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Support Vector Machine with Quadratic kernel = %s',offon(s.svmb+1,:)))



% --- Executes on button press in svmc.
function svmc_Callback(hObject, eventdata, handles)
% hObject    handle to svmc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of svmc
load Bcl_guidata
s.svmc = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Support Vector Machine with Polynomial kernel (order 3) = %s',offon(s.svmc+1,:)))


% --- Executes on button press in svmd.
function svmd_Callback(hObject, eventdata, handles)
% hObject    handle to svmd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of svmd
load Bcl_guidata
s.svmd = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Support Vector Machine with Gaussian Radial Basis Function kernel = %s',offon(s.svmd+1,:)))


% --- Executes on button press in svme.
function svme_Callback(hObject, eventdata, handles)
% hObject    handle to svme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of svme
load Bcl_guidata
s.svme = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Support Vector Machine with Multilayer Perceptron kernel = %s',offon(s.svme+1,:)))


% --- Executes on button press in pnn.
function pnn_Callback(hObject, eventdata, handles)
% hObject    handle to pnn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pnn
load Bcl_guidata
s.pnn = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Probabilistic Neural Network = %s',offon(s.pnn+1,:)))





% --- Executes on button press in nnglma.
function nnglma_Callback(hObject, eventdata, handles)
% hObject    handle to nnglma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nnglma
load Bcl_guidata
s.nnglma = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Neural Network using a Generalized Linear Model (linear) = %s',offon(s.nnglma+1,:)))



% --- Executes on button press in nnglmb.
function nnglmb_Callback(hObject, eventdata, handles)
% hObject    handle to nnglmb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nnglmb
load Bcl_guidata
s.nnglmb = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Neural Network using a Generalized Linear Model (logistic) = %s',offon(s.nnglmb+1,:)))


% --- Executes on button press in nnglmc.
function nnglmc_Callback(hObject, eventdata, handles)
% hObject    handle to nnglmc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nnglmc
load Bcl_guidata
s.nnglmc = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Neural Network using a Generalized Linear Model (softmax) = %s',offon(s.nnglmc+1,:)))


% --- Executes on button press in adaboost5.
function adaboost5_Callback(hObject, eventdata, handles)
% hObject    handle to adaboost5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of adaboost5
load Bcl_guidata
s.adaboost5 = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Adaboost Classifier with 5 iterations = %s',offon(s.adaboost5+1,:)))


% --- Executes on button press in adaboost10.
function adaboost10_Callback(hObject, eventdata, handles)
% hObject    handle to adaboost10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of adaboost10
load Bcl_guidata
s.adaboost10 = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Adaboost Classifier with 10 iterations = %s',offon(s.adaboost10+1,:)))


% --- Executes on button press in boosting025.
function boosting025_Callback(hObject, eventdata, handles)
% hObject    handle to boosting025 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of boosting025
load Bcl_guidata
s.boosting025 = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Boosting with 25%% of samples for first Classifier = %s',offon(s.boosting025+1,:)))

% --- Executes on button press in boosting050.
function boosting050_Callback(hObject, eventdata, handles)
% hObject    handle to boosting050 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of boosting050
load Bcl_guidata
s.boosting050 = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Boosting with 50%% of samples for first Classifier = %s',offon(s.boosting050+1,:)))


% --- Executes on button press in sfsfisher.
function sfsfisher_Callback(hObject, eventdata, handles)
% hObject    handle to sfsfisher (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sfsfisher
load Bcl_guidata
s.sfsfisher = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Feature Selection using SFS with Fisher Criteria = %s',offon(s.sfsfisher+1,:)))


% --- Executes on button press in fosmod.
function fosmod_Callback(hObject, eventdata, handles)
% hObject    handle to fosmod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fosmod
load Bcl_guidata
s.fosmod = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Feature Selection using FOSMOD = %s',offon(s.fosmod+1,:)))


% --- Executes on button press in lsef.
function lsef_Callback(hObject, eventdata, handles)
% hObject    handle to lsef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lsef
load Bcl_guidata
s.lsef = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Feature Selection using LSEF = %s',offon(s.lsef+1,:)))


% --- Executes on button press in frank.
function frank_Callback(hObject, eventdata, handles)
% hObject    handle to frank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of frank
load Bcl_guidata
s.frank = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Feature Selection using Rank Criteria = %s',offon(s.frank+1,:)))


% % --- Executes on button press in sfsfisherpca.
% function sfsfisherpca_Callback(hObject, eventdata, handles)
% % hObject    handle to sfsfisherpca (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%
% % Hint: get(hObject,'Value') returns toggle state of sfsfisherpca
% load Bcl_guidata
% s.sfsfisherpca = get(hObject,'Value');
% save Bcl_guidata s
% offon = ['off';'on '];
% disp(sprintf('Feature Selection from all features and PCA features using SFS with Fisher Criteria = %s',offon(s.sfsfisherpca+1,:)))


% --- Executes on button press in sfsknn5.
function sfsknn5_Callback(hObject, eventdata, handles)
% hObject    handle to sfsknn5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sfsknn5
load Bcl_guidata
s.sfsknn5 = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Feature Selection using SFS with KNN-5 Classifier = %s',offon(s.sfsknn5+1,:)))


% --- Executes on button press in sfslda.
function sfslda_Callback(hObject, eventdata, handles)
% hObject    handle to sfslda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sfslda
load Bcl_guidata
s.sfslda = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Feature Selection using SFS with LDA Classifier = %s',offon(s.sfslda+1,:)))


% --- Executes on button press in sfsqda.
function sfsqda_Callback(hObject, eventdata, handles)
% hObject    handle to sfsqda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sfsqda
load Bcl_guidata
s.sfsqda = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Feature Selection using SFS with QDA Classifier = %s',offon(s.sfsqda+1,:)))


% --- Executes on button press in sfssvmlin.
function sfssvmlin_Callback(hObject, eventdata, handles)
% hObject    handle to sfssvmlin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sfssvmlin
load Bcl_guidata
s.sfssvmlin = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Feature Selection using SFS with SVM-lin Classifier = %s',offon(s.sfssvmlin+1,:)))


% --- Executes on button press in go.
function go_Callback(hObject, eventdata, handles)
% hObject    handle to go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load Bcl_guidata
ok = 1;

if s.advanced == 0

    fclas = s.dmin         +...
        s.maha         +...
        s.lda          +...
        s.qda          +...
        s.knn4         +...
        s.knn6         +...
        s.knn8         +...
        s.knn10        +...
        s.knn15        +...
        s.knn20        +...
        s.libsvma      +...
        s.svma         +...
        s.svmb         +...
        s.svmc         +...
        s.svmd         +...
        s.svme         +...
        s.pnn          +...
        s.nnglma       +...
        s.nnglmb       +...
        s.nnglmc       +...
        s.adaboost5    +...
        s.adaboost10   +...
        s.boosting025  +...
        s.boosting050;

    ffsel = s.sfsfisher    +...
        s.fosmod       +...
        s.lsef         +...
        s.frank        +...
        s.sfsknn5      +...
        s.sfslda       +...
        s.sfsqda       +...
        s.sfslibsvmlin +...
        s.sfssvmlin    +...
        s.sfsglma      +...
        s.fsrandom     +...
        s.allfeatures;




    if ffsel==0
        beep
        wwarndlg('You must set at least one feature selector.','Error in Bcl_gui');
        ok = 0;
    end
    if fclas==0
        beep
        wwarndlg('You must set at least one classifier.','Error in Bcl_gui');
        ok = 0;
    end
end
if s.mmax<1
    beep
    wwarndlg('Invalid maximal number of features to be selected.','Error in Bcl_gui');
    ok = 0;
end
if s.crossvalfolders<1
    beep
    wwarndlg('Invalid cros validations folders.','Error in Bcl_gui');
    ok = 0;
end

if exist(s.filename,'file')==0
    beep
    wwarndlg(sprintf('.mat file %s is not defined.',s.filename),'Error in Bcl_gui');
    ok = 0;
end



if ok
    yesnoans = questdlg('Bcl_gui will start to find the best classification. This process could take several minutes. Are you sure?', ...
        'Bfex Information', ...
        'Yes', 'No', 'No');
    if yesnoans(1)=='Y'
        eval(['load ' s.filename]);
        set(handles.Result1, 'String', ' ');
        set(handles.Result2, 'String', ' ');
        pause(0)
        if not(exist('fn','var'))
            fn = [];
        end
        Bcl_guifun(s,f,fn,d,handles);
        questdlg(sprintf('Bcl_gui ended successfully.'),'Bcl_gui Information','Ok', 'Ok');
    end
end


function mmax_Callback(hObject, eventdata, handles)
% hObject    handle to mmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mmax as text
%        str2double(get(hObject,'String')) returns contents of mmax as a double
load Bcl_guidata
s.mmax = str2num(get(hObject,'String'));
save Bcl_guidata s
disp(sprintf('Maximal Number of Features to be selected = %d ',s.mmax))


% --- Executes during object creation, after setting all properties.
function mmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fsrandom.
function fsrandom_Callback(hObject, eventdata, handles)
% hObject    handle to fsrandom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fsrandom
load Bcl_guidata
s.fsrandom = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Random Feature Selection = %s',offon(s.fsrandom+1,:)))


% --- Executes on button press in allfeatures.
function allfeatures_Callback(hObject, eventdata, handles)
% hObject    handle to allfeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allfeatures
load Bcl_guidata
s.allfeatures = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('All features = %s',offon(s.allfeatures+1,:)))


% --- Executes on button press in ensvote.
function ensvote_Callback(hObject, eventdata, handles)
% hObject    handle to ensvote (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ensvote
load Bcl_guidata
s.ensvote = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Ensamble using majority vote = %s',offon(s.ensvote+1,:)))


% % --- Executes on button press in ensindividual.
% function ensindividual_Callback(hObject, eventdata, handles)
% % hObject    handle to ensindividual (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%
% % Hint: get(hObject,'Value') returns toggle state of ensindividual
% load Bcl_guidata
% s.ensindividual = get(hObject,'Value');
% save Bcl_guidata s
% offon = ['off';'on '];
% disp(sprintf('Individual Classifiers = %s',offon(s.ensindividual+1,:)))
%


function filename_Callback(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename as text
%        str2double(get(hObject,'String')) returns contents of filename as a double
load Bcl_guidata
s.filename = get(hObject,'String');
save Bcl_guidata s
disp(sprintf('MAT filename = %s',s.filename))



% --- Executes during object creation, after setting all properties.
function filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in select_filename.
function select_filename_Callback(hObject, eventdata, handles)
% hObject    handle to select_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
d = dir('*.mat');
str = {d.name};
[ss,vv] = listdlg('PromptString','Select a file:',...
    'SelectionMode','single',...
    'ListString',str);
if vv
    load Bcl_guidata
    fstr = d(ss,:).name;
    s.filename = fstr;
    save Bcl_guidata s
    disp(sprintf('MAT filename = %s',s.filename))
    load(s.filename)
    disp('Features:')
    howis(f)
    disp('Classes:')
    howis(d)
end
set(handles.printfilename, 'String', fstr);




function CrossValFolders_Callback(hObject, eventdata, handles)
% hObject    handle to CrossValFolders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CrossValFolders as text
%        str2double(get(hObject,'String')) returns contents of CrossValFolders as a double
load Bcl_guidata
s.crossvalfolders = str2num(get(hObject,'String'));
save Bcl_guidata s
disp(sprintf('Number of folders in Cross Validation = %d ',s.crossvalfolders))


% --- Executes during object creation, after setting all properties.
function CrossValFolders_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CrossValFolders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sfsglma.
function sfsglma_Callback(hObject, eventdata, handles)
% hObject    handle to sfsglma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sfsglma
load Bcl_guidata
s.sfsglma = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Feature selection using SFS with GLM-a = %s',offon(s.sfsglma+1,:)))


function [bcs,selec] = Bcl_guifun(s,f,fn,d,handles)
% s structur with classifier and feature selector definition (see Bcl_gui2.m)
% f features
% d ideal classification

if s.advanced == 0

    % Classifier defintion

    k = 0;
    if s.dmin
        k=k+1;
        bcl(k).name = 'dmin';
        bcl(k).options = [];
    end

    if s.maha
        k=k+1;
        bcl(k).name = 'maha';
        bcl(k).options = [];
    end

    if s.lda
        k=k+1;
        bcl(k).name = 'lda';
        bcl(k).options.p = [];
    end

    if s.qda
        k=k+1;
        bcl(k).name = 'qda';
        bcl(k).options.p = [];
    end

    if s.knn4
        k=k+1;
        bcl(k).name = 'knn';
        bcl(k).options.k = 4;
    end

    if s.knn6
        k=k+1;
        bcl(k).name = 'knn';
        bcl(k).options.k = 6;
    end

    if s.knn8
        k=k+1;
        bcl(k).name = 'knn';
        bcl(k).options.k = 8;
    end

    if s.knn10
        k=k+1;
        bcl(k).name = 'knn';
        bcl(k).options.k = 10;
    end

    if s.knn15
        k=k+1;
        bcl(k).name = 'knn';
        bcl(k).options.k = 15;
    end

    if s.knn20
        k=k+1;
        bcl(k).name = 'knn';
        bcl(k).options.k = 20;
    end

    if s.libsvma
        k=k+1;
        bcl(k).name = 'libsvm';
        bcl(k).options.kernel = '-t 0 -h 0';
    end
    
    if s.svma
        k=k+1;
        bcl(k).name = 'svmplus';
        bcl(k).options.kernel = 1;
    end

    if s.svmb
        k=k+1;
        bcl(k).name = 'svmplus';
        bcl(k).options.kernel = 2;
    end

    if s.svmc
        k=k+1;
        bcl(k).name = 'svmplus';
        bcl(k).options.kernel = 3;
    end

    if s.svmd
        k=k+1;
        bcl(k).name = 'svmplus';
        bcl(k).options.kernel = 4;
    end

    if s.svme
        k=k+1;
        bcl(k).name = 'svmplus';
        bcl(k).options.kernel = 5;
    end

    if s.pnn
        k=k+1;
        bcl(k).name = 'pnn';
        bcl(k).options = [];
    end

    if s.nnglma
        k=k+1;
        bcl(k).name = 'nnglm';
        bcl(k).options.method = 1;
        bcl(k).options.iter   = 12;
    end

    if s.nnglmb
        k=k+1;
        bcl(k).name = 'nnglm';
        bcl(k).options.method = 2;
        bcl(k).options.iter   = 12;
    end

    if s.nnglmc
        k=k+1;
        bcl(k).name = 'nnglm';
        bcl(k).options.method = 3;
        bcl(k).options.iter   = 12;
    end

    if s.adaboost5
        k=k+1;
        bcl(k).name = 'adaboost';
        bcl(k).options.iter   =5;
    end

    if s.adaboost10
        k=k+1;
        bcl(k).name = 'adaboost';
        bcl(k).options.iter   =10;
    end

    if s.boosting025
        k=k+1;
        bcl(k).name = 'boosting';
        bcl(k).options.s  =0.25;
    end

    if s.boosting050
        k=k+1;
        bcl(k).name = 'boosting';
        bcl(k).options.s  =0.5;
    end

    % Feature selection definition
    k = 0;
    if s.sfsfisher
        k=k+1;
        bfs(k).name = 'sfs';
        bfs(k).options.b.name = 'fisher';
    end

    if s.fosmod
        k=k+1;
        bfs(k).name = 'fosmod';
    end

    if s.lsef
        k=k+1;
        bfs(k).name = 'lsef';
    end

    if s.frank
        k=k+1;
        bfs(k).name = 'rank';
        bfs(k).options.criterion = 'ttest';
    end


    if s.sfslda
        k=k+1;
        bfs(k).name = 'sfs';
        bfs(k).options.b.name = 'lda';
        bfs(k).options.b.options.p = [];
    end

    if s.sfsqda
        k=k+1;
        bfs(k).name = 'sfs';
        bfs(k).options.b.name = 'qda';
        bfs(k).options.b.options.p = [];
    end

    if s.sfsknn5
        k=k+1;
        bfs(k).name = 'sfs';
        bfs(k).options.b.name = 'knn';
        bfs(k).options.b.options.k = 8;

    end

    if s.sfslibsvmlin
        k=k+1;
        bfs(k).name = 'sfs';
        bfs(k).options.b.name = 'libsvm';
        bfs(k).options.b.options.kernel = '-t 0 -h 0';
    end
    
    if s.sfssvmlin
        k=k+1;
        bfs(k).name = 'sfs';
        bfs(k).options.b.name = 'svmplus';
        bfs(k).options.b.options.kernel = 1;
    end

    if s.sfsglma
        k=k+1;
        bfs(k).name = 'sfs';
        bfs(k).options.b.name = 'nnglm';
        bfs(k).options.b.options.method = 1;
        bfs(k).options.b.options.iter   = 10;
    end

    if s.fsrandom
        k=k+1;
        bfs(k).name = 'random';
        bfs(k).options.b.name = 'fisher';
        bfs(k).options.M = 200;
    end

    if s.allfeatures
        k=k+1;
        bfs(k).name = 'all';
    end
else
    [bfs,bcl] = feval(s.advanced);
end



m = min(s.mmax,size(f,2));

options.Xn  = fn;
options.bcl = bcl;
options.bfs = bfs;
options.m   = m;
options.v   = s.crossvalfolders;

[bcs,selec] = Bcl_balu(f,d,bcl,bfs,options);
op.b = bcs; op.v = s.crossvalfolders; op.show = 1; op.c = 0.95;
[p,ci] = Bev_crossval(f(:,selec),d,op);

axis([0 length(selec)+1 50 100])
strper1 = sprintf('%s: %5.2f%%',bcs.name,p*100);
set(handles.Result1, 'String', strper1);

strper2 = sprintf('with %d features.',length(selec));
set(handles.Result2, 'String', strper2);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d = dir('B*.m');
str = {d.name};
[ss,vv] = listdlg('PromptString','Select a file:',...
    'SelectionMode','single',...
    'ListString',str);
if vv
    load Bcl_guidata
    fstr = d(ss,:).name;
    s.advanced = fstr(1:end-2);
    save Bcl_guidata s
    disp(sprintf('m filename = %s.m',s.advanced))
    [bfs,bcl] = feval(s.advanced);
    disp('Feature Extraction Models:')
    for i=1:length(bfs)
        fprintf(' %2d) %s-type\n',i,bfs(i).name);
    end
    disp(' ');
    disp('Classifier Models:')
    for i=1:length(bcl)
        fprintf(' %2d) %s-type\n',i,bcl(i).name);
    end
end
set(handles.text8, 'String', fstr);


% --- Executes on button press in libsvm_lin.
function libsvm_lin_Callback(hObject, eventdata, handles)
% hObject    handle to libsvm_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of libsvm_lin

load Bcl_guidata
s.libsvma = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('LIBSVM: Support Vector Machine with Linear kernel or dot product = %s',offon(s.libsvma+1,:)))





% --- Executes on button press in sfs_libsvm_lin.
function sfs_libsvm_lin_Callback(hObject, eventdata, handles)
% hObject    handle to sfs_libsvm_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sfs_libsvm_lin
load Bcl_guidata
s.sfslibsvmlin = get(hObject,'Value');
save Bcl_guidata s
offon = ['off';'on '];
disp(sprintf('Feature Selection using SFS with SVM-lin Classifier = %s',offon(s.sfslibsvmlin+1,:)))
