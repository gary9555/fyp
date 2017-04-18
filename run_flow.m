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

% Define specific operations for different datasets
translSquare = 1; % syn sample
rotBar = 2; % syn sample
translSin = 3; % real sample
rotDisk = 4;  % real sample
translBoxes = 5; % real sample
translBarX = 6; % syn sample
oscBar = 7; % syn sample

% here we assign the data set we are using 
data_name = rotDisk;
show_flow_flag = false;
%organize data in chunks
group_time = true;%true;
ms_per_frame = 80;%50; %5 %ms
id_tsb = 1; % previous group of spike index
spikes_per_frame = 400;%5000; % visualize the flow every spikes_per_frame spikes


% 3) Load data from the
if data_name == translSquare
    [x,y,pol,ts] = getDVSeventsDavis('TranslatingSquare.aedat');
elseif data_name == rotBar
    [x,y,pol,ts] = getDVSeventsDavis('RotatingBar.aedat');
elseif data_name == translSin
    [x,y,pol,ts] = getDVSeventsDavis('IMU_translSin.aedat');%,[1 1 240 180],15000);
    [t,~,~,~,~,gx,gy,gz,data] = getIMUSamplesDavis('IMU_translSin.aedat');
elseif data_name == rotDisk
    [x,y,pol,ts] = getDVSeventsDavis('IMU_rotDisk.aedat');
    [t,ax,ay,az,~,gx,gy,gz,data] = getIMUSamplesDavis('IMU_rotDisk.aedat');
elseif data_name == translBoxes
    [x,y,pol,ts] = getDVSeventsDavis('IMU_translBoxes.aedat');
    [t,ax,ay,az,temperature,gx,gy,gz,data] = getIMUSamplesDavis('IMU_translBoxes.aedat');
elseif data_name == translBarX
    [x,y,pol,ts] = getDVSeventsDavis('TranslatingBarX.aedat');
elseif data_name == oscBar
    [x,y,pol,ts] = getDVSeventsDavis('OscillatingDoubleBar.aedat');
else
    [x,y,pol,ts] = getDVSeventsDavis('RubberWhale2.aedat');
end
%[t,ax,ay,az,temperature,gx,gy,gz,data] = getIMUSamplesDavis;%('IMU_rotBars.aedat',10000);  
%ts = ts/30;

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
nt = nts;200000
% size imagedope
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

if data_name == rotBar
   load('gtRotatingBar.mat','vxGT','vyGT');
   %nts = 15000;
   vxE = zeros(180,240);
   vyE = zeros(180,240);
   update_mat = zeros(180,240);
end
if data_name == rotDisk
   tmp_cnt = zeros(180,240);
   center_x0 = 126; % initial x center of rotation
   center_y0 = 100; % initial y center of rotation
   center_vx0 = 0; % initial x speed, assuming its static at the beginning
   center_vy0 = 0;
   % apply moving average filter to smooth ax, ay
   num_avg = 40;
   nts = 170000;
   coeff = ones(1, num_avg)/num_avg;
   avg_ax = filter(coeff, 1, ax);
   avg_ay = filter(coeff, 1, ay);
   % phase shift
   avg_ax = [avg_ax(21:end) ; avg_ax(1:20)];
   avg_ay = [avg_ay(21:end) ; avg_ay(1:20)];
   avg_ax(1:19) = avg_ax(20)*ones(19,1);
   avg_ay(1:19) = avg_ay(20)*ones(19,1);
   avg_ax(end-19:end) = avg_ax(end-20)*ones(20,1);
   avg_ay(end-19:end) = avg_ay(end-20)*ones(20,1);
   omega = gz(1);
   % index of timestamp
   idx_t = 1; % for updating center
   idx_t_prev = 1; % previous idx_t considered
   idx_omega = 1; % for updating angular speed (omega)
   
end

tic
for i = 1:nts,
    % Current spike:
    xi = xs(i);
    yi = ys(i);
    tsi = ts(i);
    poli =  sign(pol(i) - .5); % +1/-1
    
    % update the belief every spike
    [u,v] = of.updateFlow(xi,yi,tsi,poli);
    
    %%%%%%%%%%%%%%%%%%% translating square %%%%%%%%%%%%%%%%%%    
    % only reasonable speed is considered
    if data_name == translSquare
        if u < speed_thres
           % u=0;
        end
        if v < speed_thres
           % v=0;
        end
        tmp = 0;
        if u<10 
            %tmp = u^2;
        else
            tmp = (20-u)^2;
        end
        if v<10
            %tmp = tmp+v^2;
        else
            tmp = tmp+ (20-v)^2; 
        end
        RAEE = RAEE + sqrt(tmp)/20;

    elseif data_name == rotBar
        if update_mat(yi,xi) == 0
            gtx = abs(vxGT(yi,xi));
            gty = abs(vyGT(yi,xi));
            if abs(gtx)> speed_thres || abs(gty)> speed_thres
                RAEE = RAEE + sqrt((abs(u)-gtx)^2+(abs(v)-gty)^2)/sqrt(gtx^2+gty^2);
            end
            update_mat(yi,xi) = 1;
            vxE(yi,xi) = u;
            vyE(yi,xi) = v;
        end

    elseif data_name == translSin


    elseif data_name == rotDisk
        tmp_cnt(yi,xi) = tmp_cnt(yi,xi)+1;
        % update center position
        if idx_t < length(t)+1 && t(idx_t) < ts(i)
            omega = gz(idx_t);
            if idx_t ~= 1 && idx_t ~= length(t)
                delta_t = double(t(idx_t) - t(idx_t-1))/1e6;
                
            elseif idx_t == 1
                delta_t = double(t(1) - ts(1)) / 1e6;
            elseif idx_t == length(t)
                if i~=nts
                    break;
                end
                delta_t = double(ts(end) - t(end)) / 1e6;
            end
            
            center_x0 = center_x0 + center_vx0*delta_t + 0.5*avg_ax(idx_t)*delta_t^2; %%%%%%%%%%%%%%
            center_y0 = center_y0 + center_vy0*delta_t + 0.5*avg_ay(idx_t)*delta_t^2; %%%%%%%%%%%%%%
            center_vx0 = center_vx0 + avg_ax(idx_t)*delta_t;
            center_vy0 = center_vy0 + avg_ay(idx_t)*delta_t;
            idx_t = idx_t + 1;
            
        end
        % calculate distance from center and then velocity
        center_x0 = 126 - 3*i/nts;
        center_y0 = 100 - 15*i/nts;
        dist = sqrt((xi-center_x0)^2 + (yi-center_y0)^2);
        
        v_gt = dist*omega; 
        theta = atan2(yi-center_y0,xi-center_x0);
        v_gtx = v_gt*cos(theta);
        v_gty = v_gt*sin(theta);       
        % update RAEE 
        RAEE = RAEE + sqrt((abs(u)-abs(v_gtx))^2+(abs(v)-abs(v_gty))^2)/sqrt(v_gtx^2+v_gty^2);
    elseif data_name == translBoxes

    end
    %%%%%%%%%%%%%%%%%%% /translating square %%%%%%%%%%%%%%%%%%    
        
    %% show flow
    if show_flow_flag == true
        U(yi,xi)= u;    
        V(yi,xi)= v;
        T(yi,xi)= tsi;

        % but only update the plot every group of spike
        if  (group_time && (tsi-ts(id_tsb) >= ms_per_frame*1E3)) || ...  % group spikes in terms of time
            (~group_time && mod(i,spikes_per_frame)==0) || (i==nts)      % group spikes in terms of number of spikes
             duration = (tsi-ts(id_tsb))/1E6; % in seconds

            %%%%%% Current normal flow vector chart %%%%%%%%
            figure(1); hold off;
            if data_name == rotDisk
                Ut = U.*(abs(U)<1000); % max THR pixel/sec
                Vt = V.*(abs(V)<1000); % max THR pixel/sec
            else
                Ut = U.*(abs(U)<5000); % max THR pixel/sec
                Vt = V.*(abs(V)<5000); % max THR pixel/sec
            end
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
                % clear the map after rendering the plot of every group of spikes
                U = zeros(ymax,xmax); V = zeros(ymax,xmax); Img= zeros(ymax,xmax);
            end
        end % update every group of spikes
    end % show_flow_flag
end % outmost loop
RAEE = RAEE/nt;
disp(['RAEE = ' num2str(RAEE)]);

%avMsXSpike = toc/i*1000;
%disp(avMsXSpike);

