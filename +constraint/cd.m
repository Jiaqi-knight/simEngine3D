classdef cd < handle
    % Filename: cd.m
    % Author:   Samuel Acuña
    % Date:     13 Oct 2016
    %
    % About:
    % This class handles instances of CD constraints. It Uses the r-p
    % formulation (euler parameters). Removes 1 degree of Freedom.
    %
    % CD constraint: 
    % the difference between the x (or y or z) coordinate of
    % point P on body i and the x (or y or z) coordinate of point Q on body
    % j assume a specified value f.
    
    properties
        system; % parent system3D object to which this constraint is a member. 
        rDOF = 1; % removes 1 degree of freedom
        cName;  % coordinate of interest {'x','y','z'}
        bodyi;  % body i
        bodyj;  % body j
        Pi;     % ID number for point P on body i
        Qj;     % ID number for point Q on body j
        f;      % prescribed constraint, often 0, but can be a function of t
        fdot;   % derivative of f
        fddot;  % derivative of fdot
    end
    properties (Dependent)
        c;      % unit vector of the coordinate of interest [3x1]
        phi;    % value of the expression of the constraint PHI^cd
        nu;     % right-hand side of the velocity equation
        gammaHat;  % right-hand side of the acceleration equation, in r-p formulation
        phi_r;  % partial derivative of constraint with respect to r
        phi_p;  % partial derivative of constraint with respect to p
    end
    
    methods
        %constructor function
        function cons = cd(system,cName,bodyi,PiID,bodyj,QjID,f,fdot,fddot,t) %constructor function
            if ~exist('f','var') || isempty(f)
                f = 0; % prescribed difference is 0
            end
            if ~exist('fdot','var') || isempty(fdot)
                fdot = 0; % derivative of ft
            end
            if ~exist('fddot','var') || isempty(fddot)
                fddot = 0; % derivative of fdt
            end
            
            cons.system = system;
            cons.cName = cName;
            cons.bodyi = bodyi;
            cons.Pi = PiID;
            cons.bodyj = bodyj;
            cons.Qj = QjID;
            cons.f = f;
            cons.fdot = fdot;
            cons.fddot = fddot;
            
            if abs(cons.phi) > 1e-4
                warning('Initial conditions for ''cd'' are not consistent. But solution will converge so constraints are satisfied.')
            end
        end
        function c = get.c(cons) % define unit vector of the coordinate of interest 
            switch cons.cName
                case 'x'
                    c = [1;0;0];
                case 'y'
                    c = [0;1;0];
                case 'z'
                    c = [0;0;1];
                otherwise
                    error('coordinate not defined in CD constraint')
            end
        end
        function phi = get.phi(cons) % value of the expression of the constraint PHI^cd
            % from ME751_f2016 slide 20 from lecture 09/26/16
            % phi : [1x1]
            
            % see if constraint is a function
            if isa(cons.f, 'function_handle')
                fVal = cons.f(cons.system.time); %evaluate at system time
            else % constant value
                fVal = cons.f;
            end
            
            phi = cons.c'*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj) - fVal;
        end
        function nu = get.nu(cons) % right-hand side of the velocity equation
            % from ME751_f2016 slide 20 from lecture 09/26/16
            % nu : [1x1]
            
            % see if constraint is a function
            if isa(cons.fdot, 'function_handle')
                fdotVal = cons.fdot(cons.system.time); %evaluate at system time
            else % constant value
                fdotVal = cons.fdot;
            end
            
            nu = fdotVal;
        end
        function gammaHat = get.gammaHat(cons) % right-hand side of the acceleration equation, in r-p formulation
            % from ME751_f2016 slide 8 from lecture 10/7/16
            % gammaHat : [1x1]
            
            % see if constraint is a function
            if isa(cons.fddot, 'function_handle')
                fddotVal = cons.fddot(cons.system.time); %evaluate at system time
            else % constant value
                fddotVal = cons.fddot;
            end
            
            gammaHat =  cons.c'*utility.Bmatrix(cons.bodyi.pdot,cons.bodyi.point{cons.Pi})*cons.bodyi.pdot + ...
                       -cons.c'*utility.Bmatrix(cons.bodyj.pdot,cons.bodyj.point{cons.Qj})*cons.bodyj.pdot + ...
                       cons.fddot;
        end
        function phi_r = get.phi_r(cons) 
            % from ME751_f2016 slide 17 from lecture 9/28/16
            % One body can be the ground. In this case, the number of columns
            % in the Jacobian is half. there are no partial derivatives with 
            % respect to r or p. Thus, we must properly dimension the size of
            % the vectors/matrices because the grounded body does not have 
            % any generalized coordinates.
            % phi_r : [1x6] normally, unless grounded, then [1x3]
            
            phi_ri = -cons.c';
            phi_rj = cons.c';
            
            if cons.bodyi.isGround
                phi_r = [phi_rj];
            elseif cons.bodyj.isGround
                phi_r = [phi_ri];
            else 
                phi_r = [phi_ri, phi_rj]; 
            end
        end
        function phi_p = get.phi_p(cons)
            % from ME751_f2016 slide 17 from lecture 9/28/16
            % One body can be the ground. In this case, the number of columns
            % in the Jacobian is half. there are no partial derivatives with 
            % respect to r or p. Thus, we must properly dimension the size of
            % the vectors/matrices because the grounded body does not have 
            % any generalized coordinates.
            % phi_p : [1x8] normally, unless grounded, then [1x4]
            
            phi_pi = -cons.c'*utility.Bmatrix(cons.bodyi.p,cons.bodyi.point{cons.Pi});
            phi_pj =  cons.c'*utility.Bmatrix(cons.bodyj.p,cons.bodyj.point{cons.Qj});
            
            if cons.bodyi.isGround
                phi_p = [phi_pj];
            elseif cons.bodyj.isGround
                phi_p = [phi_pi];
            else 
                phi_p = [phi_pi, phi_pj]; 
            end
        end
        function phiLambda_rr = phiLambda_rr(cons,lambda)
            % partial derivative of (phi_r*lambda) with respect to r (position)
            % from ME751_f2016 slide 35,39 from lecture 10/17/16
            % phiLambda_rr : [6x6] normally, unless grounded, then [3x3]
            % inputs:
            %    lambda : scalar value
            
            dim_r = length(cons.phi_r);
            phiLambda_rr = lambda*zeros(dim_r,dim_r); % [3x3] or [6x6]
        end
        function  phiLambda_rp = phiLambda_rp(cons,lambda)
            % partial derivative of (phi_r*lambda) with respect to p (orientation)
            %from ME751_f2016 slide 35,39 from lecture 10/17/16
            % phiLambda_rp : [6x8] normally, unless grounded, then [3x4]
            % inputs:
            %    lambda : scalar value
            
            dim_r = length(cons.phi_r);
            dim_p = length(cons.phi_p);
            phiLambda_rp = lambda*zeros(dim_r,dim_p); % [3x4] or [6x8]
        end
        function  phiLambda_pr = phiLambda_pr(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to r (position)
            %from ME751_f2016 slide 35,39 from lecture 10/17/16
            % phiLambda_pr : [8x6] normally, unless grounded, then [4x3]
            % inputs:
            %    lambda : scalar value
            
            dim_r = length(cons.phi_r);
            dim_p = length(cons.phi_p);
            phiLambda_pr = lambda*zeros(dim_p,dim_r); % [4x3] or [8x6]
        end
        function  phiLambda_pp = phiLambda_pp(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to p (orientation)
            %from ME751_f2016 slide 35,39 from lecture 10/17/16
            % will be [4x4] or [8x8], depends on if there is a grounded body
            % inputs:
            %    lambda : scalar value
            
            phiLambda_ppii = -utility.Kmatrix(cons.bodyi.point{cons.Pi},cons.c);
            phiLambda_ppjj =  utility.Kmatrix(cons.bodyj.point{cons.Qj},cons.c);
            
            if cons.bodyi.isGround
                phiLambda_pp = lambda*phiLambda_ppjj; % [4x4]
            elseif cons.bodyj.isGround
                phiLambda_pp = lambda*phiLambda_ppii; % [4x4]
            else % [8x8]
                phiLambda_pp = lambda*[phiLambda_ppii   zeros(4,4);
                                       zeros(4,4)      phiLambda_ppjj]; 
            end 
        end
    end
    
end