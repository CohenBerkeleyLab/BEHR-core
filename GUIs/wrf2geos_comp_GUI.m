function varargout = wrf2geos_comp_GUI(varargin)
% WRF2GEOS_COMP_GUI MATLAB code for wrf2geos_comp_GUI.fig
%      WRF2GEOS_COMP_GUI, by itself, creates a new WRF2GEOS_COMP_GUI or raises the existing
%      singleton*.
%
%      H = WRF2GEOS_COMP_GUI returns the handle to a new WRF2GEOS_COMP_GUI or the handle to
%      the existing singleton*.
%
%      WRF2GEOS_COMP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WRF2GEOS_COMP_GUI.M with the given input arguments.
%
%      WRF2GEOS_COMP_GUI('Property','Value',...) creates a new WRF2GEOS_COMP_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before wrf2geos_comp_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to wrf2geos_comp_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help wrf2geos_comp_GUI

% Last Modified by GUIDE v2.5 23-Jun-2016 13:51:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wrf2geos_comp_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @wrf2geos_comp_GUI_OutputFcn, ...
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


% --- Executes just before wrf2geos_comp_GUI is made visible.
function wrf2geos_comp_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wrf2geos_comp_GUI (see VARARGIN)

% Choose default command line output for wrf2geos_comp_GUI
handles.output = hObject;

% varargin should contain 6 inputs, in order:
%   gc_no2 profiles
%   wrf_no2_profiles
%   lon coordinates
%   lat coordinates
%   joint pressure coordinate
%   dates of each slice
handles.data.gc_no2 = varargin{1};
handles.data.wrf_no2 = varargin{2};
handles.data.lon = varargin{3};
handles.data.lat = varargin{4};
handles.data.pres = varargin{5};
handles.data.datevec = varargin{6};

% Also initialize the indicies
handles.indicies.x = 1;
handles.indicies.t = 1;

% And the ranges each index may take on
handles.indicies.xrange = [1 length(handles.data.lon)];
handles.indicies.trange = [1 length(handles.data.datevec)];

% Plus the boundaries we want for the map
handles.lonlim = [min(handles.data.lon(:))-5, max(handles.data.lon(:))+5];
handles.latlim = [min(handles.data.lat(:))-5, max(handles.data.lat(:))+5];
handles.gcline = [];
handles.makemap = true;

handles.gcprof = [];
handles.wrfprof = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wrf2geos_comp_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = wrf2geos_comp_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_xup.
function button_xup_Callback(hObject, eventdata, handles)
% hObject    handle to button_xup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.indicies.x < handles.indicies.xrange(2)
    handles.indicies.x = handles.indicies.x + 1;
    handles = update_plots(handles);
    guidata(hObject, handles);
end


% --- Executes on button press in button_xdown.
function button_xdown_Callback(hObject, eventdata, handles)
% hObject    handle to button_xdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.indicies.x > handles.indicies.xrange(1)
    handles.indicies.x = handles.indicies.x - 1;
    handles = update_plots(handles);
    guidata(hObject, handles);
end


% --- Executes on button press in button_time_forward.
function button_time_forward_Callback(hObject, eventdata, handles)
% hObject    handle to button_time_forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.indicies.t < handles.indicies.trange(2)
    handles.indicies.t = handles.indicies.t + 1;
    handles = update_plots(handles);
    guidata(hObject, handles);
end

% --- Executes on button press in button_time_backward.
function button_time_backward_Callback(hObject, eventdata, handles)
% hObject    handle to button_time_backward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.indicies.t > handles.indicies.trange(1)
    handles.indicies.t = handles.indicies.t - 1;
    handles = update_plots(handles);
    guidata(hObject, handles);
end

function handles = update_plots(handles)
% This function will update all the plots and descriptive text labels
% when called. Should be called by any of the button callback functions.
% It will only plot the state outlines on the map once, after that it will
% just redraw the GEOS-Chem grid cell.

%%% MAP %%%

if handles.makemap
    states(handles.axes_map, 'k')
    set(handles.axes_map,'XLim', handles.lonlim);
    set(handles.axes_map,'YLim', handles.latlim);
    handles.makemap = false;
end

if ~isempty(handles.gcline)
    delete(handles.gcline);
end

[gloncorn, glatcorn] = geos_chem_corners;
gloncornvec = gloncorn(1,:)';
glatcornvec = glatcorn(:,1)';
[glon, glat] = geos_chem_centers('2x25');

xx = find(glon == handles.data.lon(handles.indicies.x));
yy = find(glat == handles.data.lat(handles.indicies.x));

linex = [gloncornvec(xx), gloncornvec(xx+1), gloncornvec(xx+1), gloncornvec(xx), gloncornvec(xx)];
liney = [glatcornvec(yy), glatcornvec(yy), glatcornvec(yy+1), glatcornvec(yy+1), glatcornvec(yy)];

handles.gcline = line(linex, liney, 'color', 'r', 'linewidth', 2);
handles.text_pix_x.String = sprintf('%d of %d', handles.indicies.x, handles.indicies.xrange(2));
if glon(xx) < 0
    lonstr = '%.1f W';
else
    lonstr = '%.1f E';
end
handles.text_lon.String = sprintf(lonstr, abs(glon(xx)));
if glat(yy) < 0
    latstr = '%.1f S';
else
    latstr = '%.1f N';
end
handles.text_lat.String = sprintf(latstr, abs(glat(yy)));

handles.text_timestep.String = sprintf('%d of %d', handles.indicies.t, handles.indicies.trange(2));
handles.text_date.String = datestr(handles.data.datevec(handles.indicies.t),'yyyy-mm-dd');

%%% PROFILES %%%
if ~isempty(handles.gcprof); delete(handles.gcprof); end
if ~isempty(handles.wrfprof); delete(handles.wrfprof); end

set(handles.axes_profiles,'ydir','reverse');

handles.gcprof = line(handles.data.gc_no2(:,handles.indicies.x, handles.indicies.t), handles.data.pres(:,handles.indicies.x, handles.indicies.t), 'color', 'r', 'linewidth', 2, 'parent', handles.axes_profiles);
handles.wrfprof = line(handles.data.wrf_no2(:, handles.indicies.x, handles.indicies.t)*1e3, handles.data.pres(:,handles.indicies.x, handles.indicies.t), 'color', 'b', 'linewidth', 2, 'parent', handles.axes_profiles);

legend(handles.axes_profiles,[handles.gcprof; handles.wrfprof],{'GEOS-Chem','WRF-Chem'});
