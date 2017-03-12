classdef SpikeNormalFlow<handle
    %SPIKEOF Summary of this class goes here
    %   Compute the flow using 4 or 8 neighboors of the spike.  
    
    properties
        tdecay = 0.12*5e5;%0.004*5e5; % decay of the evidence over time
        tmax = 0.5e5; % slowest vel: 2E5 pixel/.2 sec
        tmin =  1e3; % fastest vel
        nres = 20; % Resolution of t estimates
        ndir = 4; %8
        xmax,
        ymax,
        len_buffer = 10;
        buffer,
        t_val, 
        tstd, % this parameter should be tuned
        belief,
        counter = 0, % used for debugging
        A= [1 0; 1 1; 0 1; -1 1];
          
    end
    
    methods
        %receive the size of the image
        function obj = SpikeNormalFlow(xmax, ymax)
            obj.xmax = xmax;
            obj.ymax = ymax;
            obj.buffer = cell(xmax,ymax);
            obj.belief = cell(xmax,ymax);
            obj.setPars(obj.tmin,obj.tmax,obj.nres, obj.tdecay, 0.5*obj.tmax/obj.nres);
        end
        
        function  setPars(obj,tmin,tmax,resolution, tdecay, tstd)
            obj.tmin = tmin;
            obj.tmax = tmax; 
            obj.nres = resolution;
            obj.tdecay = tdecay;
            obj.t_val = [-linspace(obj.tmax,obj.tmin,obj.nres), 0 , linspace(obj.tmin,obj.tmax,obj.nres)]; % 2*nres+1 values
            %obj.t_val = obj.tmax*(-obj.nres:obj.nres)/obj.nres;
            obj.tstd = tstd*obj.tmax/obj.nres; %obj.tstd = .5*obj.tmax/obj.nres;
            
            for x=1:obj.xmax
                for y=1:obj.ymax
                    obj.belief{x,y} = zeros(obj.ndir, size(obj.t_val,2));
                end
            end
        end
        
        function [vx,vy]=updateFlow(obj,x,y,ts,pol)
            ts  = double(ts);
            tsp = ts*pol; % Polarity encoded as positive/negative times

            % Put current event in queue for current pixel
            buf_xy = [tsp obj.buffer{x,y}];
            if (length(buf_xy) > obj.len_buffer)
                buf_xy = buf_xy(1:obj.len_buffer);
            end
            obj.buffer{x,y} = buf_xy;

            % Time since last spike:
            if length(buf_xy) > 1,
                dt_spike = abs(tsp) - abs(buf_xy(2));
            else
                dt_spike = 0;
            end

            % No need to estimate flow along outside edges of image:
            if ((x <= 1) || (x >= obj.xmax) || (y <= 1) || (y >= obj.ymax))
                vx = 0; 
                vy = 0;
                return;
            end
            % update the beliefs
            obj.belief{x,y} = obj.updateBelief(obj.belief, dt_spike,x,y,tsp,pol);
            
            %  MAP estimate
            %[tx,ty] = obj.lse(x,y);
            [tx,ty] = obj.wlseBelief(x,y);
            %[tx,ty] = obj.wlseStd(x,y);
            %% compute the flow 
            tt = max(tx^2 + ty^2, 1);
            vx = 1E6 * tx/tt; % pixel/sec
            vy = 1E6 * ty/tt;
            
            % propagate the belief
%              predict = 0.2*obj.belief{x,y};
%              if vx>0
%                  obj.belief{x+1,y} = obj.belief{x+1,y} + predict;
%                  if vy>0
%                      obj.belief{x+1,y+1} = obj.belief{x+1,y+1} + predict;
%                      obj.belief{x,y+1} = obj.belief{x,y+1} + predict;
%                  elseif vy<0
%                      obj.belief{x+1,y-1} = obj.belief{x+1,y-1} + predict;
%                      obj.belief{x,y-1} = obj.belief{x,y-1} + predict;
%                  end
%              end
%              if vy<0
%                  obj.belief{x-1,y} = obj.belief{x-1,y} + predict;
%                  if vy>0
%                      obj.belief{x-1,y+1} = obj.belief{x-1,y+1} +predict;
%                      obj.belief{x,y+1} = obj.belief{x,y+1} + predict;
%                  elseif vy<0
%                      obj.belief{x-1,y-1} = obj.belief{x-1,y-1} + predict;
%                      obj.belief{x,y-1} = obj.belief{x,y-1} + predict;
%                  end
%              end
           
        end
        
        
        
        function [tx,ty] = lse(obj,x,y)
            B =zeros(obj.ndir,1);
            bel =zeros(obj.ndir,1);
            for i=1:obj.ndir
                [val, imax] = max( obj.belief{x,y}(i,:));
                bel(i)= val;
                B(i) = obj.t_val(imax);
            end
            Txy =obj.A\B;
            %Txy = (A'*A)\(A'*B);
            tx = Txy(1);
            ty = Txy(2);
        end
        
        function [tx,ty] = wlseBelief(obj,x,y)
            B =zeros(obj.ndir,1);
            bel =zeros(obj.ndir,1);
            for i=1:obj.ndir
                [val, imax] = max( obj.belief{x,y}(i,:));
                bel(i)= val;
                B(i) = obj.t_val(imax);
            end
            
            W = diag(bel); % weights
            
            sw = sort(bel,'descend');
            
          %  disp(sw(2));
            if sw(2) > 0.7                
                Txy = (W*obj.A)\(W*B);  % gets Txy in the least square sense, such that (W*obj.A)*Txy = W*B has a least square solution Txy
                tx = Txy(1);
                ty = Txy(2);
                
            else
                %disp('zero');
                obj.counter = obj.counter+1;
                tx=0;
                ty=0;
            end
        end
        
        function [tx,ty] = wlseStd(obj,x,y)
            %f = fit(obj.t_val',obj.belief_t0{x,y}','gauss1');
            %stdx = f.c1;
            min_prob = 0.1;
            B =zeros(obj.ndir,1);
            V = max(obj.t_val).^2*ones(obj.ndir,1);
           
            for i=1:obj.ndir
                if max(obj.belief{x,y}(i,:))>=min_prob
                    n = obj.belief{x,y}(i,:)'/sum(obj.belief{x,y}(i,:));
                    B(i)  = obj.t_val*n;
                    V(i) = ( (obj.t_val-B(i) ).^2*n);
                end
            end
             
            W = diag(1./(V));
            %Txy = (A'*W*A)\(A'*W*B);
            Txy = (W*obj.A)\(W*B);
            tx = Txy(1);
            ty = Txy(2);
        end
        
        function bel = decayBelief(obj,dt_spike,bel )
            % init the belief
            if isempty(bel),
                 bel = zeros(obj.ndir, size(obj.t_val,2));
            end
            % Decay of the beliefs:
            bel = exp(-dt_spike/obj.tdecay)*bel;
        end
        
        function bel = updateBelief(obj,belief, dt_spike,x,y,tsp,pol)
            bel = belief{x,y};
           
            % size of belief{x,y} is 4*61
            % decay belief
            
            bel = decayBelief(obj,dt_spike,bel );
            
            
            %disp(bel);
            txa = tsp - pol*obj.t_val; %
            txb = tsp + pol*obj.t_val; %
            
            % for ang = 1:obj.ndir
                
            % get the two adjacent spike buffers
            % exa = zeros(obj.len_buffer,1);%*obj.ndir);
            % exb = zeros(obj.len_buffer,1);%*obj.ndir);
            nxa = zeros(obj.ndir,1); % used to store the number of buffers for each direction
            nxb = zeros(obj.ndir,1);
            
            %switch ang 
            %    case 1 % 0 deg
            xa = [x-1;x-1;x;  x+1];
            ya = [y;  y-1;y-1;y-1];
            xb = [x+1;x+1;x;  x-1];
            yb = [y;  y+1;y+1;y+1];
            %ptra = 1; ptrb =1;
            for i=1:obj.ndir
                nxa(i) = length(obj.buffer{xa(i), ya(i)});
                nxb(i) = length(obj.buffer{xb(i), yb(i)});
                if (nxa(i))>=1                   
                    %exa(ptra:ptra+ nxa(i)-1) = obj.buffer{xa(i), ya(i)}';
                    exa = obj.buffer{xa(i), ya(i)}'; % get the adjacent spike buffer
                    dt = repmat(exa,[1 length(txa)]) - repmat(txa,[nxa(i) 1]);  % create a large matrix consisting of tiling of copies of the original one
                    
                    p = exp(-.5*(dt/obj.tstd).^2); %%%%%%%%%% this parameter should be tuned 
                    bel(i,:) = bel(i,:) + sum(p,1);  % sum of each column
                end
                %ptra = ptra+nxa(i);
                if (nxb(i))>=1                    
                    %exa(ptrb:ptra+ nxb(i)-1) = obj.buffer{xb(i), yb(i)}';
                    exb = obj.buffer{xb(i), yb(i)}';
                    dt = repmat(exb,[1 length(txb)]) - repmat(txb,[nxb(i) 1]);
                    p = exp(-.5*(dt/obj.tstd).^2);        
                    bel(i,:) = bel(i,:) + sum(p,1);   
                end
                %ptrb = ptrb+nxb(i);
            end
        end
        
        
    end
    
end

