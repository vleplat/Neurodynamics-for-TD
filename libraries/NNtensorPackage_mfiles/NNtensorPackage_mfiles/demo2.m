function varargout = demo2(varargin)
% DEMO2 M-file for demo2.fig
%      DEMO2, by itself, creates a new DEMO2 or raises the existing
%      singleton*.
%
%      H = DEMO2 returns the handle to a new DEMO2 or the handle to
%      the existing singleton*.
%
%      DEMO2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEMO2.M with the given input arguments.
%
%      DEMO2('Property','Value',...) creates a new DEMO2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before demo2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to demo2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Edit the above text to modify the response to help demo2
%
% Last Modified by GUIDE v2.5 19-Jul-2011 17:00:25
%
% HELP:
% 1) generate the 3-rd order tensor first, with a desired size. 
% 2) compute the decomposition pushing the Run button.
% 3) Check the backtracking box to use an approximated line search,
%  which alternates with exact line search, and increases the number of 
% iterations to reach the convergence, but actually boosts the computation
% time for large tensors.
% 4) Add penalization L1/L2 with small alpha value, to take into account
% the sparsity of data, or smooth the results. It is often relevant when 
% real data are dealt with.
%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @demo2_OpeningFcn, ...
                   'gui_OutputFcn',  @demo2_OutputFcn, ...
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


% --- Executes just before demo2 is made visible.
function demo2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to demo2 (see VARARGIN)

% Choose default command line output for demo2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

h = findobj(gcbf, 'style', 'radio');
set(h, 'value', 0);
set(handles.ID_RadioNoPen, 'value', 1);
set(handles.ID_EditPen, 'Enable', 'Off');

rand('twister',sum(100*clock))

% UIWAIT makes demo2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = demo2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3



function ID_EditPen_Callback(hObject, eventdata, handles)
% hObject    handle to ID_EditPen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ID_EditPen as text
%        str2double(get(hObject,'String')) returns contents of ID_EditPen as a double


% --- Executes during object creation, after setting all properties.
function ID_EditPen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ID_EditPen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ID_RadioNoPen.
function ID_RadioNoPen_Callback(hObject, eventdata, handles)
% hObject    handle to ID_RadioNoPen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (get(hObject,'Value') == get(hObject,'Max'))
  % Radio button is selected-take approriate action
  h = findobj(gcbf, 'style', 'radio');
  set(h, 'value', 0);
  set(hObject, 'value', 1);
  
  set(handles.ID_EditPen,'Enable','off');
else
  % Radio button is not selected-take approriate action
  set(handles.ID_EditPen,'Enable','on');
end

% Hint: get(hObject,'Value') returns toggle state of ID_RadioNoPen


% --- Executes on button press in ID_RadioL1.
function ID_RadioL1_Callback(hObject, eventdata, handles)
% hObject    handle to ID_RadioL1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (get(hObject,'Value') == get(hObject,'Max'))
  % Radio button is selected-take approriate action
  h = findobj(gcbf, 'style', 'radio');
  set(h, 'value', 0);
  set(hObject, 'value', 1);
  
  set(handles.ID_EditPen,'Enable','on');
else
  % Radio button is not selected-take approriate action
end

% Hint: get(hObject,'Value') returns toggle state of ID_RadioL1


% --- Executes on button press in ID_RadioL2.
function ID_RadioL2_Callback(hObject, eventdata, handles)
% hObject    handle to ID_RadioL2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (get(hObject,'Value') == get(hObject,'Max'))
  % Radio button is selected-take approriate action
  h = findobj(gcbf, 'style', 'radio');
  set(h, 'value', 0);
  set(hObject, 'value', 1);
  
  set(handles.ID_EditPen,'Enable','on');
else
  % Radio button is not selected-take approriate action
end
% Hint: get(hObject,'Value') returns toggle state of ID_RadioL2


% --- Executes on button press in ID_ButtonRun.
function ID_ButtonRun_Callback(hObject, eventdata, handles)
% hObject    handle to ID_ButtonRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

T = evalin('base', 'T');

backtracking = get(handles.checkbox1,'Value');
pen = 0;
alpha = 0;
if(get(handles.ID_RadioL1,'Value') == 1)
    pen = 1;
    alpha = str2double(get(handles.ID_EditPen,'String'));
elseif(get(handles.ID_RadioL2,'Value'))
    pen = 2;
    alpha = str2double(get(handles.ID_EditPen,'String'));
end

holdPlot = get(handles.checkbox3,'Value');

nbIterations = 10000;

nbFactors = str2double(get(handles.ID_EditNbFactors,'String'));

AInit = rand(size(T, 1), nbFactors);
BInit = rand(size(T, 2), nbFactors);
CInit = rand(size(T, 3), nbFactors);

axes(handles.axes1);

%if(strcmp(get(handles.ID_ComboAlgo,'String'), 'Gradient'))
[A B C error] = cgP(T, nbFactors, createOptions(pen, alpha, backtracking, nbIterations), AInit, BInit, CInit);

if(holdPlot)
    hold on;
else
    hold off;
end
plot(10 * log10(error));

hold on

[A B C error] = gradP(T, nbFactors, createOptions(pen, alpha, backtracking, nbIterations), AInit, BInit, CInit);
plot(10 * log10(error),'r');

[A B C error] = bfgsP(T, nbFactors, createOptions(pen, alpha, backtracking, nbIterations), AInit, BInit, CInit);
plot(10 * log10(error), 'g');

legend('Conjugate gradient', 'Gradient', 'BFGS');
xlabel('Iterations')
ylabel('Reconstruction error')
    %h=gca;
%     set(handles.nc,'String',nc);
%     set(handles.nl,'String',nl);
    %imagesc(em, ex,log10(img.'),'parent',handles.axes1);
%end

function T = constructTensor(s1, s2, s3, nbFactors)
Af = rand(s1, nbFactors);
Bf = rand(s2, nbFactors);
Cf = rand(s3, nbFactors);

%Construct a positive 3rd tensor with these loadings
T = Af * khatriRao(Cf, Bf)';
T = reshape(T, s1, s2, s3);

assignin('base', 'T', T);

% --- Executes on selection change in ID_ComboAlgo.
function ID_ComboAlgo_Callback(hObject, eventdata, handles)
% hObject    handle to ID_ComboAlgo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ID_ComboAlgo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ID_ComboAlgo


% --- Executes during object creation, after setting all properties.
function ID_ComboAlgo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ID_ComboAlgo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ID_ButtonGenerate.
function ID_ButtonGenerate_Callback(hObject, eventdata, handles)
% hObject    handle to ID_ButtonGenerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dim1 = str2double(get(handles.edit1,'String'));
if(dim1 < 2), dim1 = 2; end

dim2 = str2double(get(handles.edit2,'String'));
if(dim2 < 2), dim2 = 2; end

dim3 = str2double(get(handles.edit3,'String'));
if(dim3 < 2), dim3 = 2; end

nbFactors = str2double(get(handles.ID_EditNbFactors,'String'));

constructTensor(dim1, dim2, dim3, nbFactors);



function ID_EditNbFactors_Callback(hObject, eventdata, handles)
% hObject    handle to ID_EditNbFactors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ID_EditNbFactors as text
%        str2double(get(hObject,'String')) returns contents of ID_EditNbFactors as a double


% --- Executes during object creation, after setting all properties.
function ID_EditNbFactors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ID_EditNbFactors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


