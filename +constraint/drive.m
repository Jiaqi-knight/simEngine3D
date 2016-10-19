classdef drive < handle
    % Filename: drive.m
    % Author:   Samuel Acuña
    % Date:     19 Oct 2016
    %
    % About:
    % This class handles instances of driving constraints. It Uses the r-p
    % formulation (euler parameters). 
    %
    % driving constraint: 
    % Specify motion of generalized coordinate. see pgs 380-382 in Haug Book.
    
    
    properties
        rDOF;   % removes 1 degree of freedom
        f;      % prescribed constraint, can be a function of t
        fdot;   % derivative of f
        fddot;  % derivative of fdot
    end
    properties (Dependent)
        phi;    % value of the expression of the constraint PHI^dp1
%         nu;     % right-hand side of the velocity equation
%         gammaHat;  % right-hand side of the acceleration equation, in r-p formulation
%         phi_r;  % partial derivative of constraint with respect to r
%         phi_p;  % partial derivative of constraint with respect to p
    end
    
    methods
        %constructor function
        function cons = drive(f,fdot,fddot) %constructor function
            if ~exist('f','var') || isempty(f)
                f = 0; % prescribed constraint is 0, indicating vectors are orthogonal
            end
            if ~exist('fdot','var') || isempty(fdot)
                fdot = 0; % derivative of ft
            end
            if ~exist('fddot','var') || isempty(fddot)
                fddot = 0; % derivative of fdt
            end
            cons.rDOF = 1; % removes 1 dof
            cons.f = f;
            cons.fdot = fdot;
            cons.fddot = fddot;
        end
        function phi = get.phi(cons) % value of the expression of the constraint PHI^dp1
            % phi : [1x1]
            phi = 0;
        end
%         function nu = get.nu(cons) % right-hand side of the velocity equation
%             % UNSURE
%             % nu : [1x1]?
%             nu = cons.fdot ?;
%         end
%         function gammaHat = get.gammaHat(cons) % right-hand side of the acceleration equation, in r-p formulation
%             % UNSURE
%             % gammaHat : [1x1] ?
%             gammaHat = ?
%         end
%         function phi_r = get.phi_r(cons) 
%             % UNSURE
%             % One body can be the ground. In this case, the number of columns
%             % in the Jacobian is half. there are no partial derivatives with 
%             % respect to r or p. Thus, we must properly dimension the size of
%             % the vectors/matrices because the grounded body does not have 
%             % any generalized coordinates.
%             % phi_r : ?? [1x6] normally, unless grounded, then [1x3]
%             
%             ?? phi_ri = zeros(1,3);
%             ?? phi_rj = zeros(1,3);
%             
%             if cons.bodyi.isGround
%                 phi_r = [phi_rj];
%             elseif cons.bodyj.isGround
%                 phi_r = [phi_ri];
%             else 
%                 phi_r = [phi_ri, phi_rj]; 
%             end
%         end
%         function phi_p = get.phi_p(cons)
%             % UNSURE
%             % One body can be the ground. In this case, the number of columns
%             % in the Jacobian is half. there are no partial derivatives with 
%             % respect to r or p. Thus, we must properly dimension the size of
%             % the vectors/matrices because the grounded body does not have 
%             % any generalized coordinates.
%             % phi_p : ??  [1x8] normally, unless grounded, then [1x4]
%             
%             ?? phi_pi = ?
%             ?? phi_pj = ?
%             
%             if cons.bodyi.isGround
%                 phi_p = [phi_pj];
%             elseif cons.bodyj.isGround
%                 phi_p = [phi_pi];
%             else 
%                 phi_p = [phi_pi, phi_pj]; 
%             end
%         end
    end
    
end