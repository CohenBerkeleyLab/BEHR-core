function varargout = no2_column_map_GUI(varargin)
% NO2_COLUMN_MAP_GUI MATLAB code for no2_column_map_GUI.fig
%      NO2_COLUMN_MAP_GUI, by itself, creates a new NO2_COLUMN_MAP_GUI or raises the existing
%      singleton*.
%
%      H = NO2_COLUMN_MAP_GUI returns the handle to a new NO2_COLUMN_MAP_GUI or the handle to
%      the existing singleton*.
%
%      NO2_COLUMN_MAP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NO2_COLUMN_MAP_GUI.M with the given input arguments.
%
%      NO2_COLUMN_MAP_GUI('Property','Value',...) creates a new NO2_COLUMN_MAP_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before no2_column_map_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to no2_column_map_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help no2_column_map_GUI

% Last Modified by GUIDE v2.5 01-Apr-2016 13:43:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @no2_column_map_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @no2_column_map_GUI_OutputFcn, ...
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


% --- Executes just before no2_column_map_GUI is made visible.
function no2_column_map_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to no2_column_map_GUI (see VARARGIN)

% Set the figure to use the correct close request functions
%set(handles.no2_column_map_GUI,'CloseRequestFcn',@no2_column_map_CloseRequestFcn);

% Choose default command line output for no2_column_map_GUI
handles.output = hObject;

% Draw a continental US map
usa = shaperead('usastatehi.shp');
for b=1:50
    if b~=2 && b~=11 %Don't draw AK and HI
        line(usa(b).X, usa(b).Y, 'color', 'k', 'parent', handles.map_axes);
    end
end

% Draw the default lat/lon boundaries (25 N, 125 W to 50 N, 65 W).  Save
% the line handle, it will be needed for deletion later.
latlon_line = line([-125, -125, -65, -65, -125], [25, 50, 50, 25, 25], 'color', 'r', 'linewidth',1, 'parent', handles.map_axes);
set(handles.map_axes,'YLim',[20 55]);
handles.latlon_line = latlon_line;
handles.latbdy = [25 50];
handles.lonbdy = [-125 -65];

% Set the default border color to white
handles.border_color = [1,1,1];

% Set a flag to know if the cloud fraction was set manually; if not, always
% set it to default whenever the cloud product is changed
handles.manual_clouds = false;

% Set the default value for state plotting and holidays
handles.states = 1;
handles.holidays = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes no2_column_map_GUI wait for user response (see UIRESUME)
uiwait(handles.no2_column_map_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = no2_column_map_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isfield(handles,'return')
    varargout{1} = handles.return;
    varargout{2} = handles.output;
else
    varargout{1} = handles.output;
end
delete(hObject);


% --- Executes on button press in make_map.
function make_map_Callback(hObject, eventdata, handles)
% hObject    handle to make_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.start_date) || isempty(handles.end_date)
    msgbox('Start and end dates must be filled in!','Date error','warn')
else
    try
        sdn = datenum(handles.start_date); edn = datenum(handles.end_date);
    if edn < sdn
        msgbox('Ending date is before start date.','Date error','warn')
    else
        output_struct.start_date = handles.start_date; output_struct.end_date = handles.end_date;
        output_struct.latbdy = handles.latbdy; output_struct.lonbdy = handles.lonbdy;
        output_struct.border_color = handles.border_color;
        output_struct.state = handles.states;
        output_struct.behr_dir = handles.behr_dir;
        output_struct.behr_prefix = handles.behr_prefix;
        output_struct.rowanomaly = handles.rowanomaly;
        output_struct.cloud_type = handles.cloudtype;
        output_struct.cloud_max = handles.cloud_max;
        output_struct.data_field = handles.datafield;
        output_struct.projection = handles.projection;
        output_struct.coast = handles.coast;
        output_struct.grid_resolution = handles.gridres;
        
        if handles.dayofweek==1
            output_struct.flags = {};
        else
            dayofweek = {'weekday', 'r_weekday', 'weekend', 'r_weekend'};
            output_struct.flags = dayofweek(handles.dayofweek-1);
        end
        if handles.holidays
            output_struct.flags{end+1} = 'US_holidays';
        end
        
        handles.return = output_struct;
        guidata(hObject,handles);
        no2_column_map_CloseRequestFcn(handles.no2_column_map_GUI,eventdata,handles);
    end
    catch err
        if strcmp(err.identifier,'MATLAB:datenum:ConvertDateString')
            msgbox('One of the input dates is not a valid format.','Date error','Error')
        else
            rethrow(err)
        end
    end
end

% --- Executes on button press in cancel_button.
function cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.return = 0;
guidata(hObject,handles);
no2_column_map_CloseRequestFcn(handles.no2_column_map_GUI,eventdata,handles);


function dir_text_Callback(hObject, eventdata, handles)
% hObject    handle to dir_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dir_text as text
%        str2double(get(hObject,'String')) returns contents of dir_text as a double
handles.behr_dir = get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function dir_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dir_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
behr_dir = '/Volumes/share-sat/SAT/BEHR/BEHR_Files_2014';
if ~exist(behr_dir,'dir');
    behr_dir = '';
end
set(hObject,'String',behr_dir);
handles.behr_dir = behr_dir;
guidata(hObject,handles);


function prefix_text_Callback(hObject, eventdata, handles)
% hObject    handle to prefix_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prefix_text as text
%        str2double(get(hObject,'String')) returns contents of prefix_text as a double
handles.behr_prefix = get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function prefix_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prefix_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
behr_prefix = 'OMI_BEHR_';
set(hObject,'String',behr_prefix);
handles.behr_prefix = behr_prefix;
guidata(hObject,handles)

% --- Executes on button press in getfile_button.
function getfile_button_Callback(hObject, eventdata, handles)
% hObject    handle to getfile_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.mat','Choose a file to set the BEHR directory and prefix');
p = regexp(filename,'_\d\d\d\d\d\d\d\d'); % Find the index of the filename where the prefix ends and the date starts
prefix = filename(1:p);
if filename ~= 0
    set(handles.prefix_text,'String',prefix);
    handles.behr_prefix = prefix;
    set(handles.dir_text,'String',pathname);
    handles.behr_dir = pathname;
    guidata(hObject,handles);
end

% --- Executes on selection change in dayofweek_popup.
function dayofweek_popup_Callback(hObject, eventdata, handles)
% hObject    handle to dayofweek_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dayofweek_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dayofweek_popup
handles.dayofweek = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function dayofweek_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dayofweek_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.dayofweek = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in holidays_check.
function holidays_check_Callback(hObject, eventdata, handles)
% hObject    handle to holidays_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of holidays_check
handles.holidays = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on selection change in rowanom_popup.
function rowanom_popup_Callback(hObject, eventdata, handles)
% hObject    handle to rowanom_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns rowanom_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from rowanom_popup
contents = cellstr(get(hObject,'String'));
handles.rowanomaly = contents{get(hObject,'Value')};
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function rowanom_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rowanom_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
contents = cellstr(get(hObject,'String'));
handles.rowanomaly = contents{get(hObject,'Value')};
guidata(hObject,handles);

% --- Executes on selection change in clouds_popup.
function clouds_popup_Callback(hObject, eventdata, handles)
% hObject    handle to clouds_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns clouds_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clouds_popup
contents = cellstr(get(hObject,'String'));
clouds = contents{get(hObject,'Value')};
if strcmp(clouds,'OMI') && ~handles.manual_clouds
    set(handles.cloudmax_text,'String','0.2');
    handles.cloud_max = 0.2;
elseif strcmp(clouds,'MODIS') && ~handles.manual_clouds
    set(handles.cloudmax_text,'String','0.0');
    handles.cloud_max = 0;
end
handles.cloudtype = clouds;
guidata(hObject,handles);
    

% --- Executes during object creation, after setting all properties.
function clouds_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clouds_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
contents = cellstr(get(hObject,'String'));
handles.cloudtype = contents{get(hObject,'Value')};
guidata(hObject,handles);


function cloudmax_text_Callback(hObject, eventdata, handles)
% hObject    handle to cloudmax_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cloudmax_text as text
%        str2double(get(hObject,'String')) returns contents of cloudmax_text as a double
handles.manual_clouds = true;
handles.cloud_max = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function cloudmax_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cloudmax_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
cloud_max = '0.2';
set(hObject,'String',cloud_max);
handles.cloud_max = str2double(cloud_max);
guidata(hObject,handles);

function data_text_Callback(hObject, eventdata, handles)
% hObject    handle to data_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_text as text
%        str2double(get(hObject,'String')) returns contents of data_text as a double
handles.datafield = get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function data_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
datafield = 'BEHRColumnAmountNO2Trop';
set(hObject,'String',datafield);
handles.datafield = datafield;
guidata(hObject,handles);


% --- Executes on selection change in projection_popup.
function projection_popup_Callback(hObject, eventdata, handles)
% hObject    handle to projection_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns projection_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from projection_popup
contents = cellstr(get(hObject,'String'));
handles.projection = contents{get(hObject,'Value')};
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function projection_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to projection_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
contents = cellstr(get(hObject,'String'));
handles.projection = contents{get(hObject,'Value')};
guidata(hObject,handles);

% --- Executes on button press in states_checkbox.
function states_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to states_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of states_checkbox
handles.states = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on selection change in coast_popup.
function coast_popup_Callback(hObject, eventdata, handles)
% hObject    handle to coast_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns coast_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from coast_popup
contents = cellstr(get(hObject,'String'));
handles.coast = contents{get(hObject,'Value')};
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function coast_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coast_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
contents = cellstr(get(hObject,'String'));
handles.coast = contents{get(hObject,'Value')};
guidata(hObject,handles);


function grid_text_Callback(hObject, eventdata, handles)
% hObject    handle to grid_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grid_text as text
%        str2double(get(hObject,'String')) returns contents of grid_text as a double
res = str2double(get(hObject,'String'));
handles.gridres = res;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function grid_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grid_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
res = 0.05;
set(hObject,'String',num2str(res));
handles.gridres = res;
guidata(hObject,handles);


function red_text_Callback(hObject, eventdata, handles)
% hObject    handle to red_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of red_text as text
%        str2double(get(hObject,'String')) returns contents of red_text as a double
c = getbordercolor(handles);
handles = setbordercolor(handles,c);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function red_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end



function green_edit_Callback(hObject, eventdata, handles)
% hObject    handle to green_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of green_edit as text
%        str2double(get(hObject,'String')) returns contents of green_edit as a double
c = getbordercolor(handles);
handles = setbordercolor(handles,c);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function green_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end



function blue_edit_Callback(hObject, eventdata, handles)
% hObject    handle to blue_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blue_edit as text
%        str2double(get(hObject,'String')) returns contents of blue_edit as a double
c = getbordercolor(handles);
handles = setbordercolor(handles,c);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function blue_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in colorpick_button.
function colorpick_button_Callback(hObject, eventdata, handles)
% hObject    handle to colorpick_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
color = handles.border_color;
c = uisetcolor(color,'Choose coast and state border color');
if ~isscalar(c) % Check that the user did not press cancel, returning a 0
    handles = setbordercolor(handles,c);
    guidata(hObject,handles);
end


function start_date_text_Callback(hObject, eventdata, handles)
% hObject    handle to start_date_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_date_text as text
%        str2double(get(hObject,'String')) returns contents of start_date_text as a double
handles.start_date = get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function start_date_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_date_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
start_date = '';
set(hObject,'String',start_date)
handles.start_date = start_date;
guidata(hObject,handles);


function end_date_text_Callback(hObject, eventdata, handles)
% hObject    handle to end_date_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_date_text as text
%        str2double(get(hObject,'String')) returns contents of end_date_text as a double
handles.end_date = get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function end_date_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_date_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end_date = '';
set(hObject,'String',end_date)
handles.end_date = end_date;
guidata(hObject,handles);


function lat_min_text_Callback(hObject, eventdata, handles)
% hObject    handle to lat_min_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lat_min_text as text
%        str2double(get(hObject,'String')) returns contents of lat_min_text as a double
handles = setlatlon(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function lat_min_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lat_min_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','25');



function lat_max_text_Callback(hObject, eventdata, handles)
% hObject    handle to lat_max_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lat_max_text as text
%        str2double(get(hObject,'String')) returns contents of lat_max_text as a double
handles = setlatlon(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function lat_max_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lat_max_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',50);


function lon_min_text_Callback(hObject, eventdata, handles)
% hObject    handle to lon_min_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lon_min_text as text
%        str2double(get(hObject,'String')) returns contents of lon_min_text as a double
handles = setlatlon(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function lon_min_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lon_min_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',-125);


function lon_max_text_Callback(hObject, eventdata, handles)
% hObject    handle to lon_max_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lon_max_text as text
%        str2double(get(hObject,'String')) returns contents of lon_max_text as a double
handles = setlatlon(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function lon_max_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lon_max_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',-65);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  ADDITIONAL FUNCTIONs  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = setlatlon(S)
latmin = str2double(get(S.lat_min_text,'String'));
latmax = str2double(get(S.lat_max_text,'String'));
lonmin = str2double(get(S.lon_min_text,'String'));
lonmax = str2double(get(S.lon_max_text,'String'));

if latmin > latmax;
    tmp = latmin; latmin = latmax; latmax = tmp; clear('tmp')
    set(S.lat_min_text,'String',num2str(latmin));
    set(S.lat_max_text,'String',num2str(latmax));
end

if lonmin > lonmax;
    tmp = lonmin; lonmin = lonmax; lonmax = tmp; clear('tmp')
    set(S.lon_min_text,'String',num2str(lonmin));
    set(S.lon_max_text,'String',num2str(lonmax));
end

delete(S.latlon_line);
S.latlon_line = line([lonmin, lonmin, lonmax, lonmax, lonmin], [latmin, latmax, latmax, latmin, latmin], 'color', 'r', 'linewidth', 2, 'parent', S.map_axes);
S.latbdy = [latmin, latmax];
S.lonbdy = [lonmin, lonmax];

function c = getbordercolor(S)
    c = [0,0,0];
    c(1) = str2double(get(S.red_text,'String'));
    c(2) = str2double(get(S.green_edit,'String'));
    c(3) = str2double(get(S.blue_edit,'String'));

function S = setbordercolor(S,c)
    S.border_color = c;
    set(S.red_text,'String',num2str(c(1)));
    set(S.red_text,'BackgroundColor',c)
    set(S.green_edit,'String',num2str(c(2)));
    set(S.green_edit,'BackgroundColor',c)
    set(S.blue_edit,'String',num2str(c(3)));
    set(S.blue_edit,'BackgroundColor',c)


% --- Executes when user attempts to close no2_column_map_GUI.
function no2_column_map_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to no2_column_map_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if exist('handles','var') && isstruct(handles) && isfield(handles,'no2_column_map_GUI') && isequal(get(handles.no2_column_map_GUI,'waitstatus'),'waiting')
uiresume(handles.no2_column_map_GUI);
else
delete(hObject);
end
