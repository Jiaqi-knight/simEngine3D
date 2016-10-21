classdef d < handle
    % Filename: d.m
    % Author:   Samuel Acu�a
    % Date:     14 Oct 2016
    %
    % About:
    % This class handles instances of D constraints. It Uses the r-p
    % formulation (euler parameters). Removes 1 degree of Freedom.
    %
    % D constraint: 
    % the distance between point P on body i and point Q on
    % body j assume a specified value f, f > 0. 
    %
    % Note: 
    % I'm using the constraint formulation from Haug eq 9.4.8, not
    % the way it was presented in class (ME751_f2016 slide 17 from lecture
    % 09/26/16)
    
    properties
        rDOF = 1; % removes 1 degree of freedom
        bodyi;  % body i
        bodyj;  % body j
        Pi;     % ID number for point P on body i
        Qj;     % ID number for point Q on body j
        f;      % prescribed distance constraint, > 0, and can be a function of t.
        fdot;   % derivative of f
        fddot;  % derivative of fdot
        t;      % system time
    end
    properties (Dependent)
        PiQj;   % vector in GLOBAL RF form point P on body i to point Q on body j
        phi;    % value of the expression of the constraint PHI^d
        nu;     % right-hand side of the velocity equation
        gammaHat;  % right-hand side of the acceleration equation, in r-p formulation
        phi_r;  % partial derivative of constraint with respect to r
        phi_p;  % partial derivative of constraint with respect to p
    end
    
    methods
        %constructor function
        function cons = d(bodyi,PiID,bodyj,QjID,f,fdot,fddot,t) %constructor function
            if exist('f','var') && (f <=0)
                error('D constraint cannot have a prescribed distance <= 0');
            end
            if ~exist('f','var') || isempty(f)
                f = 1; % prescribed constraint is 1
            end
            if ~exist('fdot','var') || isempty(fdot)
                fdot = 0; % derivative of ft
            end
            if ~exist('fddot','var') || isempty(fddot)
                fddot = 0; % derivative of fdt
            end
            if ~exist('t','var') || isempty(t)
                t = 0; % default time
            end
            
            cons.bodyi = bodyi;
            cons.Pi = PiID;
            cons.bodyj = bodyj;
            cons.Qj = QjID;
            cons.f = f;
            cons.fdot = fdot;
            cons.fddot = fddot;
            cons.t = t;
            
            if abs(cons.phi) > 1e-4
                warning('Initial conditions for ''d'' are not consistent. But solution will converge so constraints are satisfied.')
            end
        end
        function PiQj = get.PiQj(cons) % vector in GLOBAL RF of the constrained distance
            % PiQj = vector from Pi to Qj
            PiQj = utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj);
        end
        function phi = get.phi(cons) % value of the expression of the constraint PHI^d
            % from Haug 9.4.8
            % phi : [1x1]
            
            % see if constraint is a function
            if isa(cons.f, 'function_handle')
                fVal = cons.f(cons.t); %evaluate at time t
            else % constant value
                fVal = cons.f;
            end
            
            phi = utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)-(fVal)^2;
        end
        function nu = get.nu(cons) % right-hand side of the velocity equation
            % time derivative of f^2
            % nu : [1x1]
            
            if isa(cons.fdot, 'function_handle')
                fVal = cons.f(cons.t); %evaluate at time t
                fdotVal = cons.fdot(cons.t); %evaluate at time t
            else % constant value
                fVal = cons.f;
                fdotVal = cons.fdot;
            end
            
            nu = 2*fVal*fdotVal;
        end
        function gammaHat = get.gammaHat(cons) % right-hand side of the acceleration equation, in r-p formulation
            % from ME751_f2016 slide 8 from lecture 10/7/16, 
            % also, second time derivative of f^2
            % gammaHat : [1x1]
            
            % see if constraint is a function
            if isa(cons.fddot, 'function_handle')
                fVal = cons.f(cons.t); %evaluate at time t
                fdotVal = cons.fdot(cons.t); %evaluate at time t
                fddotVal = cons.fddot(cons.t); %evaluate at time t
            else % constant value
                fVal = cons.f;
                fdotVal = cons.fdot;
                fddotVal = cons.fddot;
            end
            
            
            gammaHat = -2*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.Bmatrix(cons.bodyj.pdot,cons.bodyj.point{cons.Qj})*cons.bodyj.pdot + ...
                        2*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.Bmatrix(cons.bodyi.pdot,cons.bodyi.point{cons.Pi})*cons.bodyi.pdot + ...
                       -2*utility.dijdot(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.dijdot(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj) + ...
                        2*fdotVal*fdotVal + 2*fVal*fddotVal;
        end
        function phi_r = get.phi_r(cons) 
            % from ME751_f2016 slide 16 from lecture 9/28/16
            % One body can be the ground. In this case, the number of columns
            % in the Jacobian is half. there are no partial derivatives with 
            % respect to r or p. Thus, we must properly dimension the size of
            % the vectors/matrices because the grounded body does not have 
            % any generalized coordinates.
            % phi_r : [1x6] normally, unless grounded, then [1x3]
            
            phi_ri = -2*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)';
            phi_rj =  2*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)';
            
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
            
            phi_pi = -2*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.Bmatrix(cons.bodyi.p,cons.bodyi.point{cons.Pi});
            phi_pj =  2*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.Bmatrix(cons.bodyj.p,cons.bodyj.point{cons.Qj});
            
            if cons.bodyi.isGround
                phi_p = [phi_pj];
            elseif cons.bodyj.isGround
                phi_p = [phi_pi];
            else 
                phi_p = [phi_pi, phi_pj]; 
            end
        end
    end
    
end