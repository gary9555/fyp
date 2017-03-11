%addpath('../optic_flow/') % load plot_flowc
addpath('./aedat/') 

%% Load dataset
% 1) dataset from Walter from moving camera
%load 'data/spikedata_0710.mat'
%x = double(dataM{1}); y = double(dataM{2}); pol = double(dataM{3}); ts = double(dataM{4});
%xmax = 128;
%ymax = 128;

% 2) Load x,y,ts,pol arrays:
%load 'data/seq_0001_e.mat'
%load 'data/seq_0020_e.mat'

% 3) Load data from spinning ball
% load 'selectedData.mat'
% x = x_selected_(15000:end);
% y = y_selected_(15000:end);
% ts = ts_selected_(15000:end);
% pol = pol_selected_(15000:end);

% 3) Load data from the 
[x,y,pol,ts] = getDVSeventsDavis('RotatingBar.aedat');
%[t,ax,ay,az,temperature,gx,gy,gz,data] = getIMUSamplesDavis;%('IMU_rotBars.aedat',10000);  
%ts = ts/30;

% Define specific operations for different datasets
translSquare = 1; % syn sample
rotBar = 2; % syn sample
translSin = 3; % real sample
rotDisk = 4;  % real sample
translBoxes = 5; % real sample
% here we assign the data set we are using 
data_name = rotBar;

%%%%%
% load 'selectedData.mat'
% x = x_selected_(15000:15200);
% y = y_selected_(15000:15200);
% ts = ts_selected_(15000:15200);
% pol = pol_selected_(15000:15200);
% [x1,y1,pol1,ts1] = getDVSeventsDavis('TranslatingSquare.aedat');
% x1 = x1(15000:15200);
% y1 = y1(15000:15200);
% ts1 = ts1(15000:15200);
% pol1 = pol1(15000:15200);
% figure(1);hold on;
% plot(x);
% plot(x1);hold off;
% figure(2);hold on;
% plot(y);
% plot(y1);hold off;
% figure(3);hold on;
% plot(ts);
% plot(ts1);hold off;
% figure(4);hold on;
% plot(pol);
% plot(pol1);hold off;


%%%%%
xmax = 240;
ymax = 180;
speed_thres = 0.01;

%% Intialize the data
% Ensure coordinates are 1-indexed:
xs = x+1;
ys = ymax-y; % y+1;
ts = double(ts);
pol =double(pol);


%num spikes
nts = length(ts);

%organize data in chunks
group_time = true;%true;
ms_per_frame = 5; %5 %ms
id_tsb = 1; % previous group of spike index
spikes_per_frame = 400;%5000; % visualize the flow every spikes_per_frame spikes

% size image
%xmax = 240;%max(xs);
%ymax = 180;%max(ys);

%of = SpikeOf(xmax,ymax);
of = SpikeNormalFlow(xmax,ymax);  %%%%%%%%% initialize the optical flow object and receive the size of the image %%%%%%
of.setPars(10000,200000,20, 1e4, 0.5);  %of.setPars(500,11000,30, 0.004*5e5, 0.5);  % default resolution is 20
% 1e6./of.t_val*duration % to compute max and min velocity in pixels per frame 

% variable to visualize the flow
[X,Y] = meshgrid(1:xmax,1:ymax);

U =zeros(ymax,xmax);
V =zeros(ymax,xmax);
Img = zeros(ymax,xmax);
T = zeros(ymax,xmax);

% average end point error
RAEE = 0;

if data_name ==rotBar
   load('gtRotatingBar.mat','vxGT','vyGT');
end

tic
for i = 1:nts,
    % Current spike:
    xi = xs(i);
    yi = ys(i);
    tsi = ts(i);
    poli =  sign(pol(i) - .5); % +1/-1
    
    % update the belief every spike
    if  (group_time && (tsi-ts(id_tsb) >= ms_per_frame*1E3))
        [u,v] = of.updateFlow(xi,yi,tsi,poli,true);
    else
        [u,v] = of.updateFlow(xi,yi,tsi,poli,false);
    end
 
%%%%%%%%%%%%%%%%%%% translating square %%%%%%%%%%%%%%%%%%    
    % only reasonable speed is considered
if data_name == translSquare
    if u < speed_thres
        u=0;
    end
    if v < speed_thres
        v=0;
    end
    
    tmp = 0;
    if u<10 
        tmp = u^2;
    else
        tmp = (20-u)^2;
    end
    if v<10
        tmp = tmp+v^2;
    else
        tmp = tmp+ (20-v)^2; 
    end
    RAEE = RAEE + sqrt(tmp)/20;

elseif data_name == rotBar
    gtx = vxGT(yi,xi);
    gty = vyGT(yi,xi);
    if abs(gtx)>0 && abs(gty)>0
        RAEE = RAEE + sqrt((u-gtx)^2+(v-gty)^2)/sqrt(gtx^2+gty^2);
    end
    
elseif data_name == translSin
    

elseif data_name == rotDisk
    

elseif data_name == translBoxes
    
end
%%%%%%%%%%%%%%%%%%% /translating square %%%%%%%%%%%%%%%%%%    
        
%% show flow
 
    U(yi,xi)= u;    
    V(yi,xi)= v;
    T(yi,xi)= tsi;
    
    % but only update the plot every group of spike
    if  (group_time && (tsi-ts(id_tsb) >= ms_per_frame*1E3)) || ...  % group spikes in terms of time
        (~group_time && mod(i,spikes_per_frame)==0) || (i==nts)      % group spikes in terms of number of spikes
         duration = (tsi-ts(id_tsb))/1E6; % in seconds
        
        %%%%%% Current normal flow vector chart %%%%%%%%
        figure(1); hold off;
        Ut = U.*(abs(U)<5000); % max THR pixel/sec
        Vt = V.*(abs(V)<5000); % max THR pixel/sec
        
        plot_flowc(X,Y,Ut,Vt,1/duration,2);
        %plot_flowc(X,Y,medfilt2(Ut),medfilt2(Vt),1/duration,2);% scale(2), sparsity
        %plot_flowc(X,Y,(Ut),-(Vt),1/duration,2);% scale, sparsity
        title (['Flow at spike ' num2str(i)]);
        
        %%%%% Current spike image in gray scale %%%%%%%%
        figure(2); hold off;
        %Img(sub2ind([ymax,xmax],y(i-spikes_per_frame+1:i)+1,x(i-spikes_per_frame+1:i)+1))=double(2*pol(i-spikes_per_frame+1:i))-1;
        Img(sub2ind([ymax,xmax],ymax+1-(y(id_tsb:i)+1),x(id_tsb:i)+1))=double(2*pol(id_tsb:i))-1;
        imagesc(Img); colormap gray;
        title('Input spikes')
        
        % load francisco datatas
%         if (~group_time && spikes_per_frame==5000)
%             load(['ball_francisco/frame_'  num2str(floor(i/spikes_per_frame),'%05d') '.mat']);
%             figure(3); hold off;
%             plot_flowc(X,Y,vx,-vy,.1,2);% scale, sparsity
%             title (['Plane fitting at spike ' num2str(i)]);          
%         end
        pause();
        
        id_tsb = i; % index of the end of the last group of spikes
        if i~= nts
            % clear the map every manipulation
            U = zeros(ymax,xmax); V = zeros(ymax,xmax); Img= zeros(ymax,xmax);
        end
    end
end
RAEE = RAEE/nts;

%avMsXSpike = toc/i*1000;
%disp(avMsXSpike);

