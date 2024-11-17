function varargout = all(varargin)
% all MATLAB code for all.fig
%      all, by itself, creates a new all or raises the existing
%      singleton*.
%
%      H = all returns the handle to a new all or the handle to
%      the existing singleton*.
%
%      all('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in all.M with the given input arguments.
%
%      all('Property','Value',...) creates a new all or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the all before all_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to all_OpeningFcn via varargin.
%
%      *See all Options on GUIDE's Tools menu.  Choose "all allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help all

% Last Modified by GUIDE v2.5 20-Mar-2024 10:52:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @all_OpeningFcn, ...
                   'gui_OutputFcn',  @all_OutputFcn, ...
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


% --- Executes just before all is made visible.
function all_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to all (see VARARGIN)

% Choose default command line output for all
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes all wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = all_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edtGridDim_Callback(hObject, eventdata, handles)
% hObject    handle to edtGridDim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtGridDim as text
%        str2double(get(hObject,'String')) returns contents of edtGridDim as a double


% --- Executes during object creation, after setting all properties.
function edtGridDim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtGridDim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtD0_Callback(hObject, eventdata, handles)
% hObject    handle to edtD0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtD0 as text
%        str2double(get(hObject,'String')) returns contents of edtD0 as a double


% --- Executes during object creation, after setting all properties.
function edtD0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtD0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtAlfa_Callback(hObject, eventdata, handles)
% hObject    handle to edtAlfa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtAlfa as text
%        str2double(get(hObject,'String')) returns contents of edtAlfa as a double


% --- Executes during object creation, after setting all properties.
function edtAlfa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtAlfa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtMmin_Callback(hObject, eventdata, handles)
% hObject    handle to edtMmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtMmin as text
%        str2double(get(hObject,'String')) returns contents of edtMmin as a double


% --- Executes during object creation, after setting all properties.
function edtMmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtMmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in drawgrid.
function drawgrid_Callback(hObject, eventdata, handles)
axes(handles.axes1);
set(handles.axes1, 'Visible', 'on');
cla;
n = get(handles.edtGridDim,'String');
n = str2num(n);
d0 = get(handles.edtD0,'String');
d0 = str2double(d0);

x = 1;y=1;
while numel(x) < n
    x(end+1) = x(end) + d0;
end
while numel(y) < n
    y(end+1) = y(end) + d0;
end

[x, y] = meshgrid(x, y);
grid_points = [x(:), y(:)];

scatter(grid_points(:, 1), grid_points(:, 2), 130, 'filled','b');
hold on;
for i = 1:size(grid_points, 1)
    text(grid_points(i, 1), grid_points(i, 2), num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'yellow');
end
% xlabel('Column');
% ylabel('Row');
title('Grid Visualization Before The Algo');
axis([min(grid_points(:,1))-1 max(grid_points(:,1))+1 min(grid_points(:,2))-1 max(grid_points(:,2))+1]);
grid on;
% hObject    handle to drawgrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
axes(handles.axes3);
set(handles.axes3, 'Visible', 'on');
cla;
n = get(handles.edtGridDim,'String');
n = str2num(n);
d0 = get(handles.edtD0,'String');
d0 = str2double(d0);
alfa = get(handles.edtAlfa,'String');
alfa = str2double(alfa);
Mmin = get(handles.edtMmin,'String');
Mmin = str2double(Mmin);
M = ones(n^2, 1);
nb_sensors = 1;

%Define Grid points
x = 1;y=1;
while numel(x) < n
    x(end+1) = x(end) + d0;
end
while numel(y) < n
    y(end+1) = y(end) + d0;
end

[x, y] = meshgrid(x, y);
grid_points = [x(:), y(:)];

scatter(grid_points(:, 1), grid_points(:, 2), 100, 'filled','b');
hold on;
for i = 1:size(grid_points, 1)
    text(grid_points(i, 1), grid_points(i, 2), num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'yellow');
end

title('Grid Visualization After Applying MAX-AVG-COV');
axis([min(grid_points(:,1))-1 max(grid_points(:,1))+1 min(grid_points(:,2))-1 max(grid_points(:,2))+1]);
grid on;
% Distance matrix 
Distance = zeros(n^2, n^2);

% Calculate distances between each grid point and all other points
for i = 1:size(grid_points, 1)
    point = grid_points(i, :);
    distances_to_point = zeros(1, n^2);
    for j = 1:size(grid_points, 1)
        distances_to_point(j) = sqrt((point(1) - grid_points(j, 1))^2 + (point(2) - grid_points(j, 2))^2);
    end
    %using condititon on distance for non grid points
    distances_to_point = distances_to_point + d0/sqrt(2);
    Distance(i, :) = distances_to_point;
end

% Display matrix Distance
disp('Distance Matrix:');
disp(Distance);

% Apply formula p = e^(-d) for every number in the matrix Distance
Detection_Matrix = exp(-alfa*Distance);

% Display the new matrix
disp('Detection Matrix:');
disp(Detection_Matrix);


% Calculate Miss Probability Matrix
Miss = 1 - Detection_Matrix;

% Display the new matrix MISS
disp('Matrix MISS:');
disp(Miss);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   MAX_AVG_COV   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('starting MAX_AVG_COV:');
step=1;
indexes = [];
grid_coordinates = [];


while true
    disp('STEP');disp(step);
    step = step + 1;
   
    % Calculate ?i for every grid point
    disp('somme i');
    sum_M = sum(Miss, 2)
   

    %place the sensor at k when sum_k is min
    indices = find(sum_M > 0)
    values_greater_than_zero = sum_M(indices)
    min_sigma_value = min(values_greater_than_zero)
    min_sigma_index = find(sum_M == min_sigma_value,1)
    %[min_sigma_value, min_sigma_index] = min(sum_M);
    scatter(grid_points(min_sigma_index, 1), grid_points(min_sigma_index, 2), 230, 'r', 'filled');

    % Update vector M
    M = M .* Miss(min_sigma_index, :)';

    % Display the updated vector M
    disp('Updated vector M:');
    disp(M);

    % Delete k-th row and column of the Miss
    Miss(min_sigma_index, :) = 0; % Set the k-th row to zero
    Miss(:, min_sigma_index) = 0; % Set the k-th column to zero

    % Display the updated Miss Matrix
    disp('Updated Miss Matrix:');
    disp(Miss);
    
    % Save the index
    indexes(end + 1) = min_sigma_index
    
     % Store the coordinates of the selected grid point
    grid_coordinates(end+1, :) = grid_points(min_sigma_index, :)
    
    % Increment the variable nb_sensors by 1
    nb_sensors = nb_sensors + 1;
    
    % Check if any condition for exiting the loop is met
    if all(M <= Mmin) || nb_sensors > n^2
        break; % Exit the loop
    end
end
indexesString = num2str(indexes); % Convert parameter j to a string
set(handles.senspos, 'String', ['Sensors positions are : ', indexesString]);
set(handles.senspos, 'Visible', 'on');
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function senspos_Callback(hObject, eventdata, handles)
% hObject    handle to senspos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of senspos as text
%        str2double(get(hObject,'String')) returns contents of senspos as a double


% --- Executes during object creation, after setting all properties.
function senspos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to senspos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in resett.
function resett_Callback(hObject, eventdata, handles)
axes(handles.axes1);
cla(handles.axes1);
set(handles.axes1, 'Visible', 'off');
axes(handles.axes3);
cla(handles.axes3);
set(handles.axes3, 'Visible', 'off');
axes(handles.axes4);
cla(handles.axes4);
set(handles.axes4, 'Visible', 'off');

set(handles.edtGridDim, 'String', '');
set(handles.edtD0, 'String', '');
set(handles.edtAlfa, 'String', '');
set(handles.edtMmin, 'String', '');
set(handles.avgobst, 'String', '');
set(handles.minobst, 'String', '');

set(handles.senspos, 'String', '');
set(handles.senspos, 'Visible', 'off');
set(handles.senspostwo, 'String', '');
set(handles.senspostwo, 'Visible', 'off');

clear all;
clc;
% hObject    handle to resett (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
axes(handles.axes4);
set(handles.axes4, 'Visible', 'on');
cla;
n = get(handles.edtGridDim,'String');
n = str2num(n);
d0 = get(handles.edtD0,'String');
d0 = str2double(d0);
alfa = get(handles.edtAlfa,'String');
alfa = str2double(alfa);
Mmin = get(handles.edtMmin,'String');
Mmin = str2double(Mmin);
M = ones(n^2, 1);
nb_sensors = 1;
%Define Grid points
x = 1;y=1;
while numel(x) < n
    x(end+1) = x(end) + d0;
end
while numel(y) < n
    y(end+1) = y(end) + d0;
end

[x, y] = meshgrid(x, y);
grid_points = [x(:), y(:)];

scatter(grid_points(:, 1), grid_points(:, 2), 100, 'filled','b');
hold on;
for i = 1:size(grid_points, 1)
    text(grid_points(i, 1), grid_points(i, 2), num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'yellow');
end

title('Grid Visualization After Applying MAX-MIN-COV');
axis([min(grid_points(:,1))-1 max(grid_points(:,1))+1 min(grid_points(:,2))-1 max(grid_points(:,2))+1]);
grid on;
% Distance matrix 
Distance = zeros(n^2, n^2);

% Calculate distances between each grid point and all other points
for i = 1:size(grid_points, 1)
    point = grid_points(i, :);
    distances_to_point = zeros(1, n^2);
    for j = 1:size(grid_points, 1)
        distances_to_point(j) = sqrt((point(1) - grid_points(j, 1))^2 + (point(2) - grid_points(j, 2))^2);
    end
    %using condititon on distance for non grid points
    distances_to_point = distances_to_point + d0/sqrt(2);
    Distance(i, :) = distances_to_point;
end

% Display matrix Distance
disp('Distance Matrix:');
disp(Distance);

% Apply formula p = e^(-d) for every number in the matrix Distance
Detection_Matrix = exp(-alfa*Distance);

% Display the new matrix
disp('Detection Matrix:');
disp(Detection_Matrix);


% Calculate Miss Probability Matrix
Miss = 1 - Detection_Matrix;

% Display the new matrix MISS
disp('Matrix MISS:');
disp(Miss);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     MAX_MIN_COV   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('starting MAX_MIN_COV:');
step=1;
indexes = [];
grid_coordinates = [];

k = randi(n^2);  
str = sprintf('k = %d', k);
disp(str);

%place at first point k
scatter(grid_points(k, 1), grid_points(k, 2), 200, 'r', 'filled');
% text(grid_points(k, 1), grid_points(k, 2), num2str(nb_sensors), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'yellow');

indexes(end + 1) = k;
grid_coordinates(end+1, :) = grid_points(k, :);


while true
    disp('STEP');disp(step);
    step = step + 1;
   
    % Update vector M
    M = M .* Miss(k, :)';
    % Display the updated vector M
    disp('Updated vector M:');
    disp(M);
    
    %Check if any condition for exiting the loop is met
    if all(M <= Mmin) || nb_sensors >= n^2
        break; 
    end
    
    % Increment the variable nb_sensors by 1
    nb_sensors = nb_sensors + 1;

    %Place sensor when Mk is max
    [max_M, max_index] = max(M);
    scatter(grid_points(max_index, 1), grid_points(max_index, 2), 200, 'r', 'filled');
%     text(grid_points(max_index, 1), grid_points(max_index, 2), num2str(nb_sensors), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'yellow');
    
    % Delete k-th row and column of the Miss
    Miss(k, :) = 0; % Set the k-th row to zero
    Miss(:, k) = 0; % Set the k-th column to zero
    
    % Display the updated Miss Matrix
    disp('Updated Miss Matrix:');
    disp(Miss);
    
    % Save the index
    indexes(end + 1) = max_index
    
     % Store the coordinates of the selected grid point
    grid_coordinates(end+1, :) = grid_points(max_index, :)
    
    %Update the value of k for the next iteration
    k = max_index;
    
end
indexesString = num2str(indexes); % Convert parameter j to a string
set(handles.senspostwo, 'String', ['Sensors positions are : ', indexesString]);
set(handles.senspostwo, 'Visible', 'on');

% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
n = get(handles.edtGridDim,'String');
n = str2num(n);
d0 = get(handles.edtD0,'String');
d0 = str2double(d0);
alfa = get(handles.edtAlfa,'String');
alfa = str2double(alfa);
Mmin = get(handles.edtMmin,'String');
Mmin = str2double(Mmin);
M = ones(n^2, 1);
nb_sensors = 1;
%%%%%%%%%%%%%%%%%%%%%%%%
obstacle_points = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create grid points
x = 1;y=1;
while numel(x) < n
    x(end+1) = x(end) + d0;
end
while numel(y) < n
    y(end+1) = y(end) + d0;
end

[x, y] = meshgrid(x, y);
grid_points = [x(:), y(:)];
axes(handles.axes3);
set(handles.axes3, 'Visible', 'on');
cla;
scatter(grid_points(:, 1), grid_points(:, 2), 130, 'filled', 'b');
hold on
for i = 1:size(grid_points, 1)
    text(grid_points(i, 1), grid_points(i, 2), num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'yellow');
end
axis([min(grid_points(:,1))-1 max(grid_points(:,1))+1 min(grid_points(:,2))-1 max(grid_points(:,2))+1]);
% xlabel('Column');
% ylabel('Row');
title('Grid Visualization After Applying MAX-AVG-COV with obstacles');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%    Introduce OBSTACLES        %%%%%%%%%%%%%%%%%%%
input_string = get(handles.avgobst, 'String');
sets = strsplit(input_string, ';');
% Loop through each set of coordinates
for i = 1:numel(sets)
% Split the current set of coordinates based on the delimiter (comma in this case)
    values = strsplit(sets{i}, ',');
        
    % Convert the strings to numbers and create the coordinate pair
    value1 = str2double(values{1});
    value2 = str2double(values{2});
    coordinate = [value1, value2];
    
    point1 = grid_points(value1,:);
    point2 = grid_points(value2,:);
    
    obstacle_x = (point1(1) + point2(1)) / 2;
    obstacle_y = (point1(2) + point2(2)) / 2;
    obstacle_coordinates = [obstacle_x, obstacle_y];
    
    obstacle_points=[obstacle_points;obstacle_coordinates];
end
% figure for algo and grid with obstacles
axes(handles.axes3);
cla(handles.axes3);
scatter(grid_points(:, 1), grid_points(:, 2), 130, 'filled', 'b'); % Plot grid points in blue
hold on
for i = 1:size(grid_points, 1)
    text(grid_points(i, 1), grid_points(i, 2), num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'yellow');
end
% Plot obstacles as rectangles
for i = 1:size(obstacle_points, 1)
    width = 0.15; % Width of the rectangle
    height = 0.25; % Height of the rectangle
    x = obstacle_points(i, 1)- width/2; % Adjust x-coordinate for rectangle center
    y = obstacle_points(i, 2) - height/2; % Adjust y-coordinate for rectangle center
    rectangle('Position', [x, y, width, height], 'FaceColor', 'k');
end
% Set axis limits to include the entire grid and obstacles
axis([min(grid_points(:,1))-1 max(grid_points(:,1))+1 min(grid_points(:,2))-1 max(grid_points(:,2))+1]);
title('Grid Visualization with Obstacles');
grid on;
set(handles.axes3, 'Visible', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    Variables and Parameters    %%%%%%%%%%%%%%%%

% Distance matrix 
Distance = zeros(n^2, n^2);

% Calculate distances between each grid point and all other points
for i = 1:size(grid_points, 1)
    point = grid_points(i, :);
    distances_to_point = zeros(1, n^2);
    for j = 1:size(grid_points, 1)
        distances_to_point(j) = sqrt((point(1) - grid_points(j, 1))^2 + (point(2) - grid_points(j, 2))^2);
    end
    %using condititon on distance for non grid points
    distances_to_point = distances_to_point + d0/sqrt(2);
    Distance(i, :) = distances_to_point;
end

% Display matrix Distance
disp('Distance Matrix:');
disp(Distance);

% Apply formula p = e^(-d) for every number in the matrix Distance
Detection_Matrix = exp(-alfa*Distance);

% Display the new matrix
disp('Detection Matrix:');
disp(Detection_Matrix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   Set Pij=0      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through all obstacles
for k = 1:size(obstacle_points,1)
    obstacle_x = obstacle_points(k, 1);
    obstacle_y = obstacle_points(k, 2);
    % Loop through all pairs of grid points
    for i = 1:size(grid_points, 1)
        for j = 1:size(grid_points, 1)
            if i ~= j % Exclude the same grid point
                % Calculate the equation of the line connecting grid point i and j
                % The equation of a line passing through points (x1, y1) and (x2, y2) is given by:
                % y - y1 = (y2 - y1) / (x2 - x1) * (x - x1)
                % In our case, x1 and y1 are the coordinates of grid point i, and x2 and y2 are the coordinates of grid point j
            
                x1 = grid_points(i, 1);
                y1 = grid_points(i, 2);
                x2 = grid_points(j, 1);
                y2 = grid_points(j, 2);
            
                % Check if the line is vertical (x1 = x2)
                if x1 == x2
                    % Check if the obstacle lies between the grid points i and j on the same vertical line
                    if obstacle_x == x1 && obstacle_y >= min(y1, y2) && obstacle_y <= max(y1, y2)
                        Detection_Matrix(i, j) = 0; % Set the detection probability to zero
                    end
                else
                    % Calculate the slope (m) of the line
                    m = (y2 - y1) / (x2 - x1);
                
                    % Calculate the y-intercept (b) of the line
                    b = y1 - m * x1;
                
                    % Check if the obstacle lies on the line between grid points i and j
                    % If the coordinates of the obstacle satisfy the equation of the line, and it lies between the grid points, set the detection probability to zero
                    if abs(obstacle_y - (m * obstacle_x + b)) < eps && obstacle_x >= min(x1, x2) && obstacle_x <= max(x1, x2)
                        Detection_Matrix(i, j) = 0; % Set the detection probability to zero
                    end
                end
            end
        end
    end
end

disp('New Detection_Matrix After Obstacles:');
disp(Detection_Matrix);

% Calculate Miss Probability Matrix
Miss = 1 - Detection_Matrix;

% Display the new matrix MISS
disp('Matrix MISS:');
disp(Miss);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     MAX_AVG_COV   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('starting MAX_AVG_COV:');
step=1;
indexes = [];
grid_coordinates = [];


while true
    disp('STEP');disp(step);
    step = step + 1;
   
    % Calculate ?i for every grid point
    disp('somme i');
    sum_M = sum(Miss, 2)
   

    %place the sensor at k when sum_k is min
    indices = find(sum_M > 0)
    values_greater_than_zero = sum_M(indices)
    min_sigma_value = min(values_greater_than_zero)
    min_sigma_index = find(sum_M == min_sigma_value,1)
    %[min_sigma_value, min_sigma_index] = min(sum_M);
    scatter(grid_points(min_sigma_index, 1), grid_points(min_sigma_index, 2), 230, 'r', 'filled');
    text(grid_points(min_sigma_index, 1), grid_points(min_sigma_index, 2), num2str(nb_sensors), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'yellow');

    % Update vector M
    M = M .* Miss(min_sigma_index, :)';

    % Display the updated vector M
    disp('Updated vector M:');
    disp(M);

    % Delete k-th row and column of the Miss
    Miss(min_sigma_index, :) = 0; % Set the k-th row to zero
    Miss(:, min_sigma_index) = 0; % Set the k-th column to zero

    % Display the updated Miss Matrix
    disp('Updated Miss Matrix:');
    disp(Miss);
    
    % Save the index
    indexes(end + 1) = min_sigma_index
    
     % Store the coordinates of the selected grid point
    grid_coordinates(end+1, :) = grid_points(min_sigma_index, :)
    
    % Increment the variable nb_sensors by 1
    nb_sensors = nb_sensors + 1;
    
    % Check if any condition for exiting the loop is met
    if all(M <= Mmin) || nb_sensors > n^2
        break; % Exit the loop
    end
end
indexesString = num2str(indexes); % Convert parameter j to a string
set(handles.senspos, 'String', ['Sensors positions are : ', indexesString]);
set(handles.senspos, 'Visible', 'on');



% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function avgobst_Callback(hObject, eventdata, handles)
% hObject    handle to avgobst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of avgobst as text
%        str2double(get(hObject,'String')) returns contents of avgobst as a double


% --- Executes during object creation, after setting all properties.
function avgobst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avgobst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
n = get(handles.edtGridDim,'String');
n = str2num(n);
d0 = get(handles.edtD0,'String');
d0 = str2double(d0);
alfa = get(handles.edtAlfa,'String');
alfa = str2double(alfa);
Mmin = get(handles.edtMmin,'String');
Mmin = str2double(Mmin);
M = ones(n^2, 1);
nb_sensors = 1;
%%%%%%%%%%%%%%%%%%%%%%%%
obstacle_points = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create grid points
x = 1;y=1;
while numel(x) < n
    x(end+1) = x(end) + d0;
end
while numel(y) < n
    y(end+1) = y(end) + d0;
end

[x, y] = meshgrid(x, y);
grid_points = [x(:), y(:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%    Introduce OBSTACLES        %%%%%%%%%%%%%%%%%%%
input_string = get(handles.minobst, 'String');
sets = strsplit(input_string, ';');

% Loop through each set of coordinates
for i = 1:numel(sets)
% Split the current set of coordinates based on the delimiter (comma in this case)
    values = strsplit(sets{i}, ',');
        
    % Convert the strings to numbers and create the coordinate pair
    value1 = str2double(values{1});
    value2 = str2double(values{2});
    coordinate = [value1, value2];
    
    point1 = grid_points(value1,:);
    point2 = grid_points(value2,:);
    
    obstacle_x = (point1(1) + point2(1)) / 2;
    obstacle_y = (point1(2) + point2(2)) / 2;
    obstacle_coordinates = [obstacle_x, obstacle_y];
    
    obstacle_points=[obstacle_points;obstacle_coordinates];
end

% figure for algo and grid with obstacles
axes(handles.axes4);
cla(handles.axes4);
scatter(grid_points(:, 1), grid_points(:, 2), 130, 'filled', 'b'); % Plot grid points in blue
hold on
for i = 1:size(grid_points, 1)
    text(grid_points(i, 1), grid_points(i, 2), num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'yellow');
end
% Plot obstacles as rectangles
for i = 1:size(obstacle_points, 1)
    width = 0.15;
    height = 0.25;
    x = obstacle_points(i, 1)- width/2;
    y = obstacle_points(i, 2) - height/2;
    rectangle('Position', [x, y, width, height], 'FaceColor', 'k');
end
axis([min(grid_points(:,1))-1 max(grid_points(:,1))+1 min(grid_points(:,2))-1 max(grid_points(:,2))+1]);
title('Grid Visualization After Applying MAX-MIN-COV with obstacles');
grid on;
set(handles.axes3, 'Visible', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    Variables and Parameters    %%%%%%%%%%%%%%%%

% Distance matrix 
Distance = zeros(n^2, n^2);

% Calculate distances between each grid point and all other points
for i = 1:size(grid_points, 1)
    point = grid_points(i, :);
    distances_to_point = zeros(1, n^2);
    for j = 1:size(grid_points, 1)
        distances_to_point(j) = sqrt((point(1) - grid_points(j, 1))^2 + (point(2) - grid_points(j, 2))^2);
    end
    %using condititon on distance for non grid points
    distances_to_point = distances_to_point + d0/sqrt(2);
    Distance(i, :) = distances_to_point;
end

% Display matrix Distance
disp('Distance Matrix:');
disp(Distance);

% Apply formula p = e^(-d) for every number in the matrix Distance
Detection_Matrix = exp(-alfa*Distance);

% Display the new matrix
disp('Detection Matrix:');
disp(Detection_Matrix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   Set Pij=0      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through all obstacles
for k = 1:size(obstacle_points,1)
    obstacle_x = obstacle_points(k, 1);
    obstacle_y = obstacle_points(k, 2);
    % Loop through all pairs of grid points
    for i = 1:size(grid_points, 1)
        for j = 1:size(grid_points, 1)
            if i ~= j % Exclude the same grid point
                % Calculate the equation of the line connecting grid point i and j
                % The equation of a line passing through points (x1, y1) and (x2, y2) is given by:
                % y - y1 = (y2 - y1) / (x2 - x1) * (x - x1)
                % In our case, x1 and y1 are the coordinates of grid point i, and x2 and y2 are the coordinates of grid point j
            
                x1 = grid_points(i, 1);
                y1 = grid_points(i, 2);
                x2 = grid_points(j, 1);
                y2 = grid_points(j, 2);
            
                % Check if the line is vertical (x1 = x2)
                if x1 == x2
                    % Check if the obstacle lies between the grid points i and j on the same vertical line
                    if obstacle_x == x1 && obstacle_y >= min(y1, y2) && obstacle_y <= max(y1, y2)
                        Detection_Matrix(i, j) = 0; % Set the detection probability to zero
                    end
                else
                    % Calculate the slope (m) of the line
                    m = (y2 - y1) / (x2 - x1);
                
                    % Calculate the y-intercept (b) of the line
                    b = y1 - m * x1;
                
                    % Check if the obstacle lies on the line between grid points i and j
                    % If the coordinates of the obstacle satisfy the equation of the line, and it lies between the grid points, set the detection probability to zero
                    if abs(obstacle_y - (m * obstacle_x + b)) < eps && obstacle_x >= min(x1, x2) && obstacle_x <= max(x1, x2)
                        Detection_Matrix(i, j) = 0; % Set the detection probability to zero
                    end
                end
            end
        end
    end
end

disp('New Detection_Matrix After Obstacles:');
disp(Detection_Matrix);

% Calculate Miss Probability Matrix
Miss = 1 - Detection_Matrix;

% Display the new matrix MISS
disp('Matrix MISS:');
disp(Miss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     MAX_MIN_COV   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('starting MAX_MIN_COV:');
step=1;
indexes = [];
grid_coordinates = [];

k = randi(n^2);  
str = sprintf('k = %d', k);
disp(str);

%place at first point k
scatter(grid_points(k, 1), grid_points(k, 2), 200, 'r', 'filled');
text(grid_points(k, 1), grid_points(k, 2), num2str(nb_sensors), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'yellow');

indexes(end + 1) = k;
grid_coordinates(end+1, :) = grid_points(k, :);


while true
    disp('STEP');disp(step);
    step = step + 1;
   
    % Update vector M
    M = M .* Miss(k, :)';
    % Display the updated vector M
    disp('Updated vector M:');
    disp(M);
    
    %Check if any condition for exiting the loop is met
    if all(M <= Mmin) || nb_sensors >= n^2
        break; 
    end
    
    % Increment the variable nb_sensors by 1
    nb_sensors = nb_sensors + 1;

    %Place sensor when Mk is max
    [max_M, max_index] = max(M);
    scatter(grid_points(max_index, 1), grid_points(max_index, 2), 200, 'r', 'filled');
    text(grid_points(max_index, 1), grid_points(max_index, 2), num2str(nb_sensors), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'yellow');
    
    % Delete k-th row and column of the Miss
    Miss(k, :) = 0; % Set the k-th row to zero
    Miss(:, k) = 0; % Set the k-th column to zero
    
    % Display the updated Miss Matrix
    disp('Updated Miss Matrix:');
    disp(Miss);
    
    % Save the index
    indexes(end + 1) = max_index
    
     % Store the coordinates of the selected grid point
    grid_coordinates(end+1, :) = grid_points(max_index, :)
    
    %Update the value of k for the next iteration
    k = max_index;
    
end
indexesString = num2str(indexes); % Convert parameter j to a string
set(handles.senspostwo, 'String', ['Sensors positions are : ', indexesString]);
set(handles.senspostwo, 'Visible', 'on');

% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function minobst_Callback(hObject, eventdata, handles)
% hObject    handle to minobst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minobst as text
%        str2double(get(hObject,'String')) returns contents of minobst as a double


% --- Executes during object creation, after setting all properties.
function minobst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minobst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
