function varargout = FeroResyanto_1301154318(varargin)
% UNTITLED MATLAB code for untitled.fig
%      UNTITLED, by itself, creates a new UNTITLED or raises the existing
%      singleton*.
%
%      H = UNTITLED returns the handle to a new UNTITLED or the handle to
%      the existing singleton*.
%
%      UNTITLED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED.M with the given input arguments.
%
%      UNTITLED('Property','Value',...) creates a new UNTITLED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled

% Last Modified by GUIDE v2.5 09-May-2018 08:22:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled_OutputFcn, ...
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


% --- Executes just before untitled is made visible.
function untitled_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled (see VARARGIN)

% Choose default command line output for untitled
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = untitled_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- pushbutton2 Picture ---
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[nama_file,nama_path] = uigetfile({'*.jpg';'*.jpeg';'*.bmp';'*.png';'*.tif';},...
    'Pilih Gambar');
if ~isequal (nama_file,0)
    handles.data1 = imread(fullfile(nama_path,nama_file));
    [filepath,name,ext] = fileparts(nama_file) 
    gambar =  strcat(nama_path,nama_file);
    infoGambar = imfinfo(gambar);
    guidata(hObject,handles);
    gambar =  handles.data1;

    axes(handles.axes1);
    imshow(handles.data1);
    set(handles.name_image,'string',name);
    set(handles.height,'string',infoGambar.Height);
    set(handles.width,'string',infoGambar.Width);
    
else
    return
end


function name_image_Callback(hObject, eventdata, handles)
% hObject    handle to name_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of name_image as text
%        str2double(get(hObject,'String')) returns contents of name_image as a double


% --- Executes during object creation, after setting all properties.
function name_image_CreateFcn(hObject, eventdata, handles)
% hObject    handle to name_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function height_Callback(hObject, eventdata, handles)
% hObject    handle to height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of height as text
%        str2double(get(hObject,'String')) returns contents of height as a double


% --- Executes during object creation, after setting all properties.
function height_CreateFcn(hObject, eventdata, handles)
% hObject    handle to height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function width_Callback(hObject, eventdata, handles)
% hObject    handle to width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of width as text
%        str2double(get(hObject,'String')) returns contents of width as a double


% --- Executes during object creation, after setting all properties.
function width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- GrayScale or RGB Picture
% --- Executes on selection change in grayscale.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gambar= handles.data1;
idMode = get(hObject,'value');
if idMode == 1

elseif idMode == 2 
    axes(handles.axes2);
    gray = .300*gambar(:,:,1) + .400*gambar(:,:,2) + .300*gambar(:,:,3);
    imshow(gray);
elseif idMode == 3
    axes(handles.axes2);
    imshow(gambar);
end

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% 
% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Inverse Picture ---
% --- Executes on button press in inverse.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gambar= handles.data1;
axes(handles.axes2);
hasil = 255-gambar;
imshow(hasil);


% --- Zoom In Picture ---
% --- Executes on button press in zoomin.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles.data1;
size_x = size(a,1)*2;
size_y = size(a,2)*2;
d = uint8(zeros(size_x,size_y,3));

%Proses ZoomIn sebanyak 1 kali
for b_asli=1:size(a,1)
    for k_asli=1:size(a,2)
        temp = a(b_asli,k_asli,:);
        for b_baru=1:2
            for k_baru=1:2
                dummyb=((b_asli-1)*2 + b_baru);
                dummyk= ((k_asli-1)*2 + k_baru);
                d(dummyb,dummyk,:)=temp;
            end
        end
        d(dummyb,dummyk);
    end
end

figure
imshow(a); title('Citra Sebelum ZoomIn');
figure
imshow(d); title('Citra Setelah ZoomIn');

% --- zoomout Picture ---
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles.data1;
size_x = round(size(a,1)/2);
size_y = round(size(a,2)/2);
d = uint8(zeros(size_x,size_y,3));

% Proses zoomout Citra Sebanyak 1 kali
m = 1; n = 1; x=2;
for i = 1:size(d,1)
    for j = 1:size(d,2)
        d(i,j,:) = a(m,n,:);
        n = round(n+x);
    end
    m = round(m+x);
    n = 1;
end

% figure
% imshow(a); title('Citra Sebelum ZoomOut');
figure
imshow(d); title('Citra Setelah ZoomOut');

% --- Executes on button press in Crop.
function Crop_Callback(hObject, eventdata, handles)
% hObject    handle to Crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gambar = handles.data1;
size_x = size(gambar,1);
size_y = size(gambar,2);
crop = uint8(zeros(size_x,size_y,3));
height_first = str2double(get(handles.height_crop,'String'));
height_last = str2double(get(handles.height_crop_to,'String'));
width_first = str2double(get(handles.width_crop,'String'));
width_last = str2double(get(handles.width_crop_to,'String'));

for i = height_first:height_last
    for j=width_first:width_last
        crop(i,j,:) = gambar(i,j,:);
    end
end

figure
imshow(crop); title('Citra Setelah Crop');

function height_crop_Callback(hObject, eventdata, handles)
% hObject    handle to height_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of height_crop as text
%        str2double(get(hObject,'String')) returns contents of height_crop as a double


% --- Executes during object creation, after setting all properties.
function height_crop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to height_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function width_crop_Callback(hObject, eventdata, handles)
% hObject    handle to width_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of width_crop as text
%        str2double(get(hObject,'String')) returns contents of width_crop as a double


% --- Executes during object creation, after setting all properties.
function width_crop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to width_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function height_crop_to_Callback(hObject, eventdata, handles)
% hObject    handle to height_crop_to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of height_crop_to as text
%        str2double(get(hObject,'String')) returns contents of height_crop_to as a double


% --- Executes during object creation, after setting all properties.
function height_crop_to_CreateFcn(hObject, eventdata, handles)
% hObject    handle to height_crop_to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function width_crop_to_Callback(hObject, eventdata, handles)
% hObject    handle to width_crop_to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of width_crop_to as text
%        str2double(get(hObject,'String')) returns contents of width_crop_to as a double


% --- Executes during object creation, after setting all properties.
function width_crop_to_CreateFcn(hObject, eventdata, handles)
% hObject    handle to width_crop_to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in flip_vertical.
function flip_vertical_Callback(hObject, eventdata, handles)
% hObject    handle to flip_vertical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gambar = handles.data1;
x = size(gambar);
flip = gambar;
k = x(1);
for i = 1:x(1)
    for j = 1:x(2)
        flip(k,j,:) = gambar(i,j,:);
    end
    k = k-1;
end
axes(handles.axes2);
imshow(flip);


% --- Executes on button press in flip_horizontal.
function flip_horizontal_Callback(hObject, eventdata, handles)
% hObject    handle to flip_horizontal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gambar = handles.data1;
x = size(gambar);
flip = gambar;
for i = 1:x(1)
    k = x(2);
    for j = 1:x(2)
        flip(i,k,:) = gambar(i,j,:);
        k = k-1;
    end
end
axes(handles.axes2);
imshow(flip);

function angka_Callback(hObject, eventdata, handles)
% hObject    handle to angka (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angka as text
%        str2double(get(hObject,'String')) returns contents of angka as a double


% --- Executes during object creation, after setting all properties.
function angka_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angka (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in add_brightness.
function add_brightness_Callback(hObject, eventdata, handles)
% hObject    handle to add_brightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


gambar = handles.data1;
x = str2double(get(handles.angka,'String'));
add = gambar+x;
axes(handles.axes2);
imshow(add);


% --- Executes on button press in minus_brightness.
function minus_brightness_Callback(hObject, eventdata, handles)
% hObject    handle to minus_brightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gambar = handles.data1;
x = str2double(get(handles.angka,'String'));
minus = gambar-x;
axes(handles.axes2);
imshow(minus);


% --- Executes on button press in cross_brightness.
function cross_brightness_Callback(hObject, eventdata, handles)
% hObject    handle to cross_brightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gambar = handles.data1;
x = str2double(get(handles.angka,'String'));
add = gambar*x;
axes(handles.axes2);
imshow(add);

% --- Executes on button press in divide_brightness.
function divide_brightness_Callback(hObject, eventdata, handles)
% hObject    handle to divide_brightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gambar = handles.data1;
x = str2double(get(handles.angka,'String'));
minus = gambar/x;
axes(handles.axes2);
imshow(minus);


% --- Executes on button press in histogram_rgb.
function histogram_rgb_Callback(hObject, eventdata, handles)
% hObject    handle to histogram_rgb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gambar = handles.data1;
red = gambar(:,:,1);
[M,N] = size(red);

t=1:256;
n=0:255;
count=0;

for z=1:256
    for i=1:M
        for j=1:N
            
            if red(i,j)==z-1
                count=count+1;
            end
        end
    end
    t(z)=count;
    count=0;
end
% axes(handles.axes2);
% disp(t')
figure
stem(n,t); 
grid on;
ylabel('number of pixels with such intensity levels---->');
xlabel('intensity value---->'); title('Image Histogram (RED)')


green = gambar(:,:,2);
[M,N] = size(green);

t=1:256;
n=0:255;
count=0;

for z=1:256
    for i=1:M
        for j=1:N
            
            if green(i,j)==z-1
                count=count+1;
            end
        end
    end
    t(z)=count;
    count=0;
end

figure
stem(n,t); 
grid on;
ylabel('number of pixels with such intensity levels---->');
xlabel('intensity value---->'); title('Image Histogram (green)')


blue = gambar(:,:,3);
[M,N] = size(blue);

t=1:256;
n=0:255;
count=0;

for z=1:256
    for i=1:M
        for j=1:N
            
            if blue(i,j)==z-1
                count=count+1;
            end
        end
    end
    t(z)=count;
    count=0;
end

figure
stem(n,t); 
grid on;
ylabel('number of pixels with such intensity levels---->');
xlabel('intensity value---->'); title('Image Histogram (blue)')





% --- Executes on button press in histogram_gray.
function histogram_gray_Callback(hObject, eventdata, handles)
% hObject    handle to histogram_gray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gambar = handles.data1;
h=rgb2gray(gambar);
[M,N]=size(h);

t=1:256;
n=0:255;
count=0;

for z=1:256
    for i=1:M
        for j=1:N
            
            if h(i,j)==z-1
                count=count+1;
            end
        end
    end
    t(z)=count;
    count=0;
end
% axes(handles.axes2);
% disp(t')
figure
stem(n,t); 
grid on;
ylabel('number of pixels with such intensity levels---->');
xlabel('intensity value---->'); title('Image Histogram')


% --- Executes on button press in rotate.
function rotate_Callback(hObject, eventdata, handles)
% hObject    handle to rotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gambar = handles.data1;
x = str2double(get(handles.rotate_value,'String'));

[rowsi,colsi,z]= size(gambar); 

angle=x;

rads=2*pi*angle/360;  

%calculating array dimesions such that  rotated image gets fit in it exactly.
% we are using absolute so that we get  positve value in any case ie.,any quadrant.

rowsf=ceil(rowsi*abs(cos(rads))+colsi*abs(sin(rads)));                      
colsf=ceil(rowsi*abs(sin(rads))+colsi*abs(cos(rads)));                     

% define an array withcalculated dimensionsand fill the array  with zeros ie.,black
rotate_img=uint8(zeros([rowsf colsf 3 ]));

%calculating center of original and final image
xo=ceil(rowsi/2);                                                            
yo=ceil(colsi/2);

midx=ceil((size(rotate_img,1))/2);
midy=ceil((size(rotate_img,2))/2);

% in this loop we calculate corresponding coordinates of pixel of A 
% for each pixel of C, and its intensity will be  assigned after checking
% weather it lie in the bound of A (original image)
for i=1:size(rotate_img,1)
    for j=1:size(rotate_img,2)                                                       

         x= (i-midx)*cos(rads)+(j-midy)*sin(rads);                                       
         y= -(i-midx)*sin(rads)+(j-midy)*cos(rads);                             
         x=round(x)+xo;
         y=round(y)+yo;

         if (x>=1 && y>=1 && x<=size(gambar,1) &&  y<=size(gambar,2) ) 
              rotate_img(i,j,:)=gambar(x,y,:);  
         end

    end
end

axes(handles.axes2);
imshow(rotate_img);



function rotate_value_Callback(hObject, eventdata, handles)
% hObject    handle to rotate_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rotate_value as text
%        str2double(get(hObject,'String')) returns contents of rotate_value as a double


% --- Executes during object creation, after setting all properties.
function rotate_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rotate_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sharp.
function sharp_Callback(hObject, eventdata, handles)
% hObject    handle to sharp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x = handles.data1;
[tinggi, lebar, panjang]=size(x);
image_baru=double(x);
y = [0 -1 0
    -1 5 -1
    0 -1 0];
for p=1 : 3
    for i=2 : tinggi-2
        for j=2 : lebar-2
            jum=image_baru(i-1,j-1,p)*y(1,1)+image_baru(i,j-1,p)*y(2,1)...
                + image_baru(i+1,j-1,p)*y(3,1)+image_baru(i-1,j,p)*y(1,2)...
                + image_baru(i,j,p)*y(2,2)+image_baru(i+1,j,p)*y(3,2)...
                +image_baru(i-1,j+1,p)*y(1,3)+image_baru(i,j+1,p)*y(2,3)...
                +image_baru(i+1,j+1,p)*y(3,3);
            img(i-1,j-1,p)=jum;
        end
    end
end
axes(handles.axes2);
imshow(uint8(img));


% --- Executes on button press in blur.
function blur_Callback(hObject, eventdata, handles)
% hObject    handle to blur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x = handles.data1;
[tinggi, lebar, panjang]=size(x);
m = 1/9;
image_baru=double(x);
y = [m m m
    m m m
    m m m];
for p=1 : 3
    for i=2 : tinggi-2
        for j=2 : lebar-2
            jum=image_baru(i-1,j-1,p)*y(1,1)+image_baru(i,j-1,p)*y(2,1)...
                + image_baru(i+1,j-1,p)*y(3,1)+image_baru(i-1,j,p)*y(1,2)...
                + image_baru(i,j,p)*y(2,2)+image_baru(i+1,j,p)*y(3,2)...
                +image_baru(i-1,j+1,p)*y(1,3)+image_baru(i,j+1,p)*y(2,3)...
                +image_baru(i+1,j+1,p)*y(3,3);
            img(i-1,j-1,p)=jum;
        end
    end
end
axes(handles.axes2);
imshow(uint8(img));
% figure
% imshow(uint8(img)); title('Blur Image');


% --- Executes on button press in edge.
function edge_Callback(hObject, eventdata, handles)
% hObject    handle to edge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x = handles.data1;
[tinggi, lebar, panjang]=size(x);
image_baru=double(x);
y = [1 1 1
    1 -8 1
    1 1 1];
for i=2 : tinggi-2
    for j=2 : lebar-2
        jum=image_baru(i-1,j-1)*y(1,1)+image_baru(i,j-1)*y(2,1)...
            + image_baru(i+1,j-1)*y(3,1)+image_baru(i-1,j)*y(1,2)...
            + image_baru(i,j)*y(2,2)+image_baru(i+1,j)*y(3,2)...
            +image_baru(i-1,j+1)*y(1,3)+image_baru(i,j+1)*y(2,3)...
            +image_baru(i+1,j+1)*y(3,3);
        img(i-1,j-1)=jum;
    end
end
axes(handles.axes2);
imshow(uint8(img));



function x1_Callback(hObject, eventdata, handles)
% hObject    handle to x1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x1 as text
%        str2double(get(hObject,'String')) returns contents of x1 as a double


% --- Executes during object creation, after setting all properties.
function x1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x2_Callback(hObject, eventdata, handles)
% hObject    handle to x2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x2 as text
%        str2double(get(hObject,'String')) returns contents of x2 as a double


% --- Executes during object creation, after setting all properties.
function x2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x3_Callback(hObject, eventdata, handles)
% hObject    handle to x3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x3 as text
%        str2double(get(hObject,'String')) returns contents of x3 as a double


% --- Executes during object creation, after setting all properties.
function x3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x4_Callback(hObject, eventdata, handles)
% hObject    handle to x4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x4 as text
%        str2double(get(hObject,'String')) returns contents of x4 as a double


% --- Executes during object creation, after setting all properties.
function x4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x5_Callback(hObject, eventdata, handles)
% hObject    handle to x5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x5 as text
%        str2double(get(hObject,'String')) returns contents of x5 as a double


% --- Executes during object creation, after setting all properties.
function x5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x6_Callback(hObject, eventdata, handles)
% hObject    handle to x6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x6 as text
%        str2double(get(hObject,'String')) returns contents of x6 as a double


% --- Executes during object creation, after setting all properties.
function x6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x7_Callback(hObject, eventdata, handles)
% hObject    handle to x7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x7 as text
%        str2double(get(hObject,'String')) returns contents of x7 as a double


% --- Executes during object creation, after setting all properties.
function x7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x8_Callback(hObject, eventdata, handles)
% hObject    handle to x8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x8 as text
%        str2double(get(hObject,'String')) returns contents of x8 as a double


% --- Executes during object creation, after setting all properties.
function x8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x9_Callback(hObject, eventdata, handles)
% hObject    handle to x9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x9 as text
%        str2double(get(hObject,'String')) returns contents of x9 as a double


% --- Executes during object creation, after setting all properties.
function x9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in custom_convolution.
function custom_convolution_Callback(hObject, eventdata, handles)
% hObject    handle to custom_convolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x = handles.data1;
[tinggi, lebar, panjang]=size(x);
image_baru=double(x);
x1 = str2double(get(handles.x1,'String'));
x2 = str2double(get(handles.x2,'String'));
x3 = str2double(get(handles.x3,'String'));
x4 = str2double(get(handles.x4,'String'));
x5 = str2double(get(handles.x5,'String'));
x6 = str2double(get(handles.x6,'String'));
x7 = str2double(get(handles.x7,'String'));
x8 = str2double(get(handles.x8,'String'));
x9 = str2double(get(handles.x9,'String'));

y = [x1 x2 x3
    x4 x5 x6
    x7 x8 x9];

for p=1 : 3
    for i=2 : tinggi-2
        for j=2 : lebar-2
            jum=image_baru(i-1,j-1,p)*y(1,1)+image_baru(i,j-1,p)*y(2,1)...
                + image_baru(i+1,j-1,p)*y(3,1)+image_baru(i-1,j,p)*y(1,2)...
                + image_baru(i,j,p)*y(2,2)+image_baru(i+1,j,p)*y(3,2)...
                +image_baru(i-1,j+1,p)*y(1,3)+image_baru(i,j+1,p)*y(2,3)...
                +image_baru(i+1,j+1,p)*y(3,3);
            img(i-1,j-1,p)=jum;
        end
    end
end
axes(handles.axes2);
imshow(uint8(img));


% --- Executes on button press in median_filter.
function median_filter_Callback(hObject, eventdata, handles)
% hObject    handle to median_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x = handles.data1;
[tinggi, lebar, panjang]=size(x);
img=double(x);
array = {}

for p=1 : panjang
   for i=2 : tinggi
       for j=2 : lebar
            array{1} = img(i-1,j-1,p);
            array{2} = img(i-1,j,p);
            array{3} = img(i,j-1,p);
            array{4} = img(i,j,p);
            for loop = 1 : 3
                for y=1 : 3
                    if ( array{y} >= array{y+1} )
                      con = array{y};
                      array{y} = array{y+1};
                      array{y+1} = con;
                    end
                end
            end
            median = round(( array{2}+array{3} ) / 2)
            image(i-1,j-1,p)=median;
       end
       array{1}= img(i-1,lebar,p);
       array{2} = 0; 
       array{3}= img(i,lebar,p);
       array{4} = 0;
       if ( array{1} < array{3} )
          con = array{1}; 
       else
           con = array{3};
       end   
%        for y=1 : 3
%             if ( array{y} >= array{y+1} )
%               con = array{y};
%               array{y} = array{y+1};
%               array{y+1} = con;
%             end
%        end
       median = round(( con + 0 ) / 2);
       image(i-1,lebar,p) = median;
   end
   for k=2 : lebar
        array{1} = img(tinggi,k-1,p);
        array{2} = img(tinggi,k,p);
        if ( array{1} < array{2} )
           con = array{1}; 
        else
           con = array{2};
        end
        median = round(( con + 0 ) / 2);
        image(tinggi,k-1,p) = median;
   end
end
axes(handles.axes2);
imshow(uint8(image));
            


% --- Executes on button press in mean_filter.
function mean_filter_Callback(hObject, eventdata, handles)
% hObject    handle to mean_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x = handles.data1;
[tinggi, lebar, panjang]=size(x);
img=double(x);

for p=1 : panjang
    for i=2 : tinggi
        for j=2 : lebar
           jum = img(i-1,j-1,p)+img(i-1,j,p)+img(i,j-1,p)+img(i,j,p);
           mean = jum/4;
           image(i-1,j-1,p)=mean;
        end
        jum = img(i-1,lebar,p)+img(i,lebar,p);
        mean = jum/4;
        image(i-1,lebar,p) = mean
    end
    for k=2 : lebar
        jum = img(tinggi,k-1,p)+img(tinggi,k,p);
        mean = jum/4;
        image(tinggi,k-1,p) = mean;
    end
end

axes(handles.axes2);
imshow(uint8(image));

% --- Executes on button press in modus_filter.
function modus_filter_Callback(hObject, eventdata, handles)
% hObject    handle to modus_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x = handles.data1;
[tinggi, lebar, panjang]=size(x);
img=double(x);
array = {};
banyak = {};


for p=1 : panjang
   for i=2 : tinggi
       for j=2 : lebar
            a=0;
            array{1} = img(i-1,j-1,p);
            array{2} = img(i-1,j,p);
            array{3} = img(i,j-1,p);
            array{4} = img(i,j,p);
%             sorting
            for loop=1 : 3 
                for y=1 : 3
                    if ( array{y} >= array{y+1} )
                      con = array{y};
                      array{y} = array{y+1};
                      array{y+1} = con;
                    end
                end
            end
%             hitung jumlah 
            for b=1 : 4
                banyak{b}=0;
                for c=1 : 4
                   if ( array{b} == array{c} )
                      banyak{b} = banyak{b}+1;
                   end
                end
            end
%             cari modus
            for d=1 : 3
                if ( banyak{d} >= a )
                    a = banyak{d};
                    modus = array{d};
                end
            end
            if ( a == 1)
                image(i-1,j-1,p)=array{1};
            else
                image(i-1,j-1,p)=modus;
            end
       end
       image(i-1,lebar,p) = 0;
   end
   for k=1 : lebar
        image(tinggi,k,p) = 0;
   end
end
axes(handles.axes2);
imshow(uint8(image));


% --- Executes on selection change in pilihan_segmentasi.
function pilihan_segmentasi_Callback(hObject, eventdata, handles)
% hObject    handle to pilihan_segmentasi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pilihan_segmentasi contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pilihan_segmentasi


% --- Executes during object creation, after setting all properties.
function pilihan_segmentasi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pilihan_segmentasi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Segmentasi.
function Segmentasi_Callback(hObject, eventdata, handles)
% hObject    handle to Segmentasi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

img = handles.data1;
x = get(handles.pilihan_segmentasi,'Value');

% Bikin matriks zero
white = uint8(zeros(size(img)));
blue = uint8(zeros(size(img)));
yellow = uint8(zeros(size(img)));
red = uint8(zeros(size(img)));

% Ambil warnanya
for i = 1:size(img,1)
    for j = 1:size(img,2)
        if ((img(i,j,1) <= 75) && (img(i,j,2) <= 115) && (img(i,j,3) >= 140))
            blue(i,j,:) = img(i,j,:);
        end
        if ((img(i,j,1) >= 130) && (img(i,j,2) <= 80) && (img(i,j,3) <= 80))
            red(i,j,:) = img(i,j,:);
        end
        if ((img(i,j,1) >= 200) && (img(i,j,2) >= 150) && (img(i,j,2) <= 190) && (img(i,j,3) >= 50) && (img(i,j,3) <= 90))
            yellow(i,j,:) = img(i,j,:);
        end
        if ((img(i,j,1) >= 210) && (img(i,j,2) >= 210) && (img(i,j,3) >= 210))
            white(i,j,:) = img(i,j,:);
        end
    end
end

if x == 1
    
elseif x == 2
    axes(handles.axes2);
    imshow(blue);
elseif x == 3
    axes(handles.axes2);
    imshow(red);
elseif x == 4
    axes(handles.axes2);
    imshow(yellow);
elseif x == 5
    axes(handles.axes2);
    imshow(white);
end

function green_treshold_Callback(hObject, eventdata, handles)
% hObject    handle to green_treshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of green_treshold as text
%        str2double(get(hObject,'String')) returns contents of green_treshold as a double


% --- Executes during object creation, after setting all properties.
function green_treshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_treshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function red_treshold_Callback(hObject, eventdata, handles)
% hObject    handle to red_treshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of red_treshold as text
%        str2double(get(hObject,'String')) returns contents of red_treshold as a double


% --- Executes during object creation, after setting all properties.
function red_treshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_treshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function blue_treshold_Callback(hObject, eventdata, handles)
% hObject    handle to blue_treshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blue_treshold as text
%        str2double(get(hObject,'String')) returns contents of blue_treshold as a double


% --- Executes during object creation, after setting all properties.
function blue_treshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_treshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function blue_to_treshold_Callback(hObject, eventdata, handles)
% hObject    handle to blue_to_treshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blue_to_treshold as text
%        str2double(get(hObject,'String')) returns contents of blue_to_treshold as a double


% --- Executes during object creation, after setting all properties.
function blue_to_treshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_to_treshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function green_to_treshold_Callback(hObject, eventdata, handles)
% hObject    handle to green_to_treshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of green_to_treshold as text
%        str2double(get(hObject,'String')) returns contents of green_to_treshold as a double


% --- Executes during object creation, after setting all properties.
function green_to_treshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_to_treshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function red_to_treshold_Callback(hObject, eventdata, handles)
% hObject    handle to red_to_treshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of red_to_treshold as text
%        str2double(get(hObject,'String')) returns contents of red_to_treshold as a double


% --- Executes during object creation, after setting all properties.
function red_to_treshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_to_treshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in treshold.
function treshold_Callback(hObject, eventdata, handles)
% hObject    handle to treshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

img = handles.data1;

% Ambil nilai input RGB
r1 = str2double(get(handles.red_treshold,'String'));
r2 = str2double(get(handles.red_to_treshold,'String'));
g1 = str2double(get(handles.green_treshold,'String'));
g2 = str2double(get(handles.green_to_treshold,'String'));
b1 = str2double(get(handles.blue_treshold,'String'));
b2 = str2double(get(handles.blue_to_treshold,'String'));


% Bikin matriks zero
new_img = uint8(zeros(size(img)));

% Ambil warnanya
for i = 1:size(img,1)
    for j = 1:size(img,2)
        if ((img(i,j,1) >= r1) && (img(i,j,1) <= r2) && (img(i,j,2) >= g1) && (img(i,j,2) <= g2) && (img(i,j,3) >= b1) && (img(i,j,3) <= b2))
            new_img(i,j,:) = 255;
        end
    end
end

axes(handles.axes2);
imshow(new_img);


% --- Executes on button press in dilation.
function dilation_Callback(hObject, eventdata, handles)
% hObject    handle to dilation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Buka Gambar 
img = handles.data1;
bw_img=im2bw(img);

%Structuring element
dil=[1 0 0; 0 1 0; 0 0 1];

pad=padarray(bw_img,[1 1]);
dil_img=false(size(bw_img));
for i=1:size(pad,1)-2
    for j=1:size(pad,2)-2
        dil_img(i,j)=sum(sum(dil&pad(i,j:j+2)));
    end
end

axes(handles.axes2);
imshow(dil_img);

% --- Executes on button press in erosion.
function erosion_Callback(hObject, eventdata, handles)
% hObject    handle to erosion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Buka Gambar 
img = handles.data1
bw_img=im2bw(img);

%Structuring element

dil=[1 1 0];

%Pad array with ones on both sides

pad=padarray(bw_img,[0 1],1);

%Intialize the matrix D of size A with zeros

dil_img=false(size(bw_img));

for i=1:size(pad,1)
    for j=1:size(pad,2)-2
        l=pad(i,j:j+2);
        %Find the position of ones in the structuring element
        k=find(dil==1);
       if(l(k)==1)
        dil_img(i,j)=1;
        end
    end
end

axes(handles.axes2);
imshow(dil_img);


% --- Executes on button press in compress_rle.
function compress_rle_Callback(hObject, eventdata, handles)
% hObject    handle to compress_rle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Buka Gambar 
img = handles.data1

% Transform matrices
dct_matrix = dctmtx(8);
dct = @(X)dct_matrix * X * dct_matrix';
idct = @(X)dct_matrix' * X * dct_matrix;

% Quantization tables
q_max = 255;
q_y = ...
    [16 11 10 16 24 40 51 61;
     12 12 14 19 26 58 60 55;
     14 13 16 24 40 57 69 56;
     14 17 22 29 51 87 80 62;
     18 22 37 56 68 109 103 77;
     24 35 55 64 81 104 113 92;
     49 64 78 87 103 121 120 101;
     72 92 95 98 112 100 103 99];
 q_c = ...
     [17 18 24 47 99 99 99 99;
      18 21 26 66 99 99 99 99;
      24 26 56 99 99 99 99 99;
      47 66 99 99 99 99 99 99;
      99 99 99 99 99 99 99 99;
      99 99 99 99 99 99 99 99;
      99 99 99 99 99 99 99 99;
      99 99 99 99 99 99 99 99];
  
%   Scale quantization matrices based on quality factor
qf = 30;

q_scale = floor(5000/qf);

q_y = round(q_y*q_scale/100);
q_c = round(q_c*q_scale/100);

% RGB to YCbCr
ycc = rgb2ycbcr(im2double(img));

% Down-sample and decimate chroma
cb = conv2(ycc(:,:,2),[1 1;1 1])./4.0;
cr = conv2(ycc(:,:,3),[1 1;1 1])./4.0;
cb = cb(2:2:size(cb,1),2:2:size(cb,2));
cr = cr(2:2:size(cr,1),2:2:size(cr,2));
y = ycc(:,:,1);

% DCT, with scaling before quantization
y = blkproc(y,[8 8],dct).*q_max;
cb = blkproc(cb,[8 8],dct).*q_max;
cr = blkproc(cr,[8 8],dct).*q_max;

% Quantize DCT coefficients
y = blkproc(y,[8 8],@(block_struct)round(round(block_struct)./q_y));
cb = blkproc(cb,[8 8],@(block_struct)round(round(block_struct)./q_c));
cr = blkproc(cr,[8 8],@(block_struct)round(round(block_struct)./q_c));

% Dequantize DCT coefficients
y = blkproc(y,[8 8],@(block_struct)round(round(block_struct).*q_y));
cb = blkproc(cb,[8 8],@(block_struct)round(round(block_struct).*q_c));
cr = blkproc(cr,[8 8],@(block_struct)round(round(block_struct).*q_c));

% Inverse DCT
y = blkproc(y./q_max,[8 8],idct);
cb = blkproc(cb./q_max,[8 8],idct);
cr = blkproc(cr./q_max,[8 8],idct);

% Up-sample chroma
upsample_filter_1d = [1 3 3 1]/4;
upsample_filter = upsample_filter_1d'* upsample_filter_1d;

cb = conv2(upsample_filter,upsample(upsample(padarray(cb,[1 1],'replicate'),2)',2)');
cb = cb(4:size(cb,1)-4,4:size(cb,2)-4);

cr = conv2(upsample_filter,upsample(upsample(padarray(cr,[1 1],'replicate'),2)',2)');
cr = cr(4:size(cr,1)-4,4:size(cr,2)-4);

% Concatenate the channels
jpeg_result = ycbcr2rgb(cat(3,y,cb,cr));

axes(handles.axes2);
imshow(jpeg_result);
