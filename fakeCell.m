classdef fakeCell
    % This is a classdef used to create a fake cell that can be used for
    % testing different analysis programs or hypothesis testing. Each fake
    % cell has:
    %       type (ex: position only, P),
    %       centre of firing field,
    %       spread of the firing field,
    %       optionally we can specify: number of bins [P V R]
    %       angle (for the Mix type)
    %     obj = fakeCell(type, dimensions, centre, spread, bins, Palpha, angle)
    % Aman Saleem
    % July 2013
    
    properties
        type;   % Options are P,V,R,PV,PR,VR,PVR,Mix,PMix
        dimensions;
        bins;
        
        Pcen; % where is the centre of the firing field
        Pspread; % Spread of the firing field
        Palpha;
        
        Vcen; % where is the centre of the firing field
        Vspread; % Spread of the firing field
        Rcen; % where is the centre of the firing field
        Rspread; % Spread of the firing field
        
        Mcen; % where is the centre of the firing field
        Mspread; % Spread of the firing field
        Mangle; % angle for the mix type of cells
        
        meanRate;
        response;  % This is the actual response of the fake cell
    end
    methods
        function obj = fakeCell(varargin)
            %   Contruct a class of fakeCell
            %   Usage: obj = fakeCell(type, centre, spread, bins, angle)
            %       type (ex: position only, P),
            %       centre of firing field (P has 5 centres),
            %       spread of the firing field,
            %       optionally we can specify: number of bins [P V R]
            %       angle (for the Mix type)
            pnames = {'type' 'dimensions' 'bins' ...
                'Pcen' 'Pspread' 'Palpha'...
                'Vcen' 'Vspread' ...
                'Rcen' 'Rspread' ...
                'Mcen' 'Mspread' 'Mangle' ...
                };
            dflts = {'P' 1 [150 30 30] ...
                75 30 1 ...
                30 15 ...
                30 15 ...
                30 15 45 ...
                };
            
            [obj.type, obj.dimensions, obj.bins, ...
                obj.Pcen, obj.Pspread, obj.Palpha, obj.Vcen, obj.Vspread, ...
                obj.Rcen, obj.Rspread, obj.Mcen, obj.Mspread, obj.Mangle] = ...
                internal.stats.parseArgs(pnames,dflts,varargin{:});
        end
        function obj = makeFakeCell(obj)
            % Calculate the P tuning
            P_bins = 1:obj.bins(1);
            if length(obj.Palpha) < length(obj.Pcen)
                obj.Palpha(end+1:length(obj.Pcen)) = obj.Palpha(1);
            elseif length(obj.Palpha) > length(obj.Pcen)
                error(message('Palpha cannot be longer than Pcen'))
            end
            if length(obj.Pspread) < length(obj.Pcen)
                obj.Pspread(end+1:length(obj.Pcen)) = obj.Pspread(1);
            elseif length(obj.Pspread) > length(obj.Pcen)
                error(message('Pspread cannot be longer than Pcen'))
            end
            
            for n = 1:length(obj.Pcen)
                P(n,:) = obj.Palpha(n).*exp((-(P_bins-obj.Pcen(n)).^2)./(2*obj.Pspread(n).^2));
            end
            P = nansum(P,1);
            
            % Calculate the V tuning
            V_bins = 1:obj.bins(2);
            if obj.Vcen & obj.Vspread
                V = exp((-(V_bins-obj.Vcen).^2)./(2*obj.Vspread^2));
            else
                V = ones(1,length(V_bins));
            end
            
            % Calculate the R tuning
            R_bins = 1:obj.bins(3);
            if obj.Rcen & obj.Rspread
                R = exp((-(R_bins-obj.Rcen).^2)./(2*obj.Rspread^2));
            else
                R = ones(1,length(R_bins));
            end
            % Calculate the M tuning
            M_bins = 1:obj.bins(3);
            if obj.Mcen & obj.Mspread
                M = exp((-(M_bins-obj.Mcen).^2)./(2*obj.Mspread^2));
            else
                M = ones(1,length(M_bins));
            end
            
            switch obj.type
                case 'P'
                    obj.response = P;
                case 'V'
                    obj.response = V;
                case 'R'
                    obj.response = R;
                case 'PV'
                    obj.response = P'*V;
                case 'PR'
                    obj.response = P'*R;
                case 'VR'
                    obj.response = V'*R;
                case 'PVR'
                    R_new = repmat(R, [bins(1)*bins(2) 1]);
                    R_new = reshape(R_new, bins(1), bins(2), bins(3));
                    
                    PV = P'*V;
                    PV_new = repmat(PV, [1 1 bins(3)]);
                    
                    obj.response = PV_new.*R_new;
                case 'M'
                    obj.response = M;
                case 'PM'
                    obj.response = P'*M;
            end
%             obj.response = obj.response./nansum(obj.response(:));
            obj.response = obj.response./25;
        end
    end
end
