function [hc, hq, hm] = plot_flowc(U,V,Du,Dv,scale,rate,Pu,Pv)
% PLOT_FLOWC
%   Du          - U-component of flow vectors.
%   Dv          - V-component of flow vectors.
%   [scale]     - Scale length of vectors.
%                 [ scale=1 ]
%   [rate]      - Subsample vector field with rate.
%                 [ rate=1 ]
%
% RETURN
%   hc      - Handle to colored plot.
%   hq      - Handle to quiver plot.
%
% DESCRIPTION
%   This function plots an flow field over a colored map coding the
%   direction of the local flow vector.
%

if nargin < 5 || isempty(scale), scale = 1; end
if nargin < 6 || isempty(rate), rate = 1; end
theta = pi;
[h w] = size(U);
C = ones(h,w,3);
C(:,:,1) = (atan2(Dv,Du) + pi)/(2*pi);
C(:,:,2) = sqrt(Du.^2+Dv.^2)/(max(max(sqrt(Du.^2+Dv.^2))) + eps);
% C(C>pi) = C(C>pi) - 2*pi;
% mapping to indizes
% Ci = round((C+pi)/(2*pi)*length(hsv));
colormap('hsv');
hc=image(U(:), V(:), hsv2rgb(C));
hold on;
axis ij; axis equal; 
axis off; axis tight;
axis manual;
U = U(rate:rate:end,1:rate:end);
V = V(rate:rate:end,1:rate:end);


Duu = Du(rate:rate:end,1:rate:end)/scale;
Dvv = Dv(rate:rate:end,1:rate:end)/scale;
hq = quiver(U,V,Duu,Dvv,2,'k');  %%%%% plot the vectors
%hq = quiver(U,V,Duu,Dvv,'k'); 
set(hq,'LineWidth',1.0);
%set(hq,'LineWidth',2.0);
set(hq,'MarkerSize',2.0);
%set(hq,'MarkerSize',4.0);

hm = plot_colormap;
xlabel('x [pixel]'); ylabel('y [pixel]');
%hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if nargin==8 && ~isempty(Pu) &&  ~isempty(Pv)    
    Pu = sqrt(Pu(rate:rate:end,1:rate:end));
    Pv = sqrt(Pv(rate:rate:end,1:rate:end));

    K= size(U,1)*size(U,2);

    thr=linspace(0,2*pi,40) ;

    allx = nan*ones(1, 40*K+(K-1)) ;
    ally = nan*ones(1, 40*K+(K-1)) ;

    allxf = nan*ones(1, 3*K) ;
    allyf = nan*ones(1, 3*K) ;

    for k=1:K
        xc = U(k) ;
        yc = V(k) ;
   
        %r= sqrt(Pu(k) + Pv(k)) ;
  
        x = Pu(k)*sin(thr) + xc ;
        y = Pv(k)*cos(thr) + yc ;
        %   th = atan2(Dv(k),Du(k)); 

        allx((k-1)*(41) + (1:40)) = x ;
        ally((k-1)*(41) + (1:40)) = y ;

        %allxf((k-1)*3 + (1:2)) = [xc xc+r*cos(th)] ;
        %allyf((k-1)*3 + (1:2)) = [yc yc+r*sin(th)] ;

    end

    h=line([allx nan allxf], [ally nan allyf], 'Color', [0.4 0.4 0.4],'LineWidth',1) ;
end
hold off;