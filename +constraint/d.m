classdef d < handle
    % Filename: d.m
    % Author:   Samuel Acuña
    % Date:     14 Oct 2016
    % About:
    % This class handles instances of D constraints. It Uses the r-p
    % formulation (euler parameters). Removes 1 degree of Freedom.
    % D constraint: the distance between point P on body i and point Q on
    % body j assume a specified value f, f > 0.
    properties
        bodyi;  % body i
        bodyj;  % body j
        Pi;     % ID number for point P on body i
        Qj;     % ID number for point Q on body j
        f;      % prescribed constraint, > 0, and can be a function of t
        fdot;   % derivative of f
        fddot;  % derivative of fdot
    end
    properties (Dependent)
        phi;    % value of the expression of the constraint PHI^dp1
        nu;     % right-hand side of the velocity equation
        gammaHat;  % right-hand side of the acceleration equation, in r-p formulation
        phi_r;  % partial derivative of constraint with respect to r
        phi_p;  % partial derivative of constraint with respect to p
    end
    
    methods
        %constructor function
        function cons = d(bodyi,PiID,bodyj,QjID,f,fdot,fddot) %constructor function
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
            
            cons.bodyi = bodyi;
            cons.Pi = PiID;
            cons.bodyj = bodyj;
            cons.Qj = QjID;
            cons.f = f;
            cons.fdot = fdot;
            cons.fddot = fddot;
        end
        function phi = get.phi(cons) % value of the expression of the constraint PHI^dp1
            % from ME751_f2016 slide 17 from lecture 09/26/16
            % phi : [1x1]
            phi = utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)-cons.f;
        end
        function nu = get.nu(cons) % right-hand side of the velocity equation
            % from ME751_f2016 slide 18 from lecture 09/26/16
            % nu : [1x1]
            nu = cons.fdot;
        end
        function gammaHat = get.gammaHat(cons) % right-hand side of the acceleration equation, in r-p formulation
            % from ME751_f2016 slide 8 from lecture 10/7/16
            % gammaHat : [1x1]
            gammaHat = -2*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.Bmatrix(cons.bodyj.pdot,cons.bodyj.point{cons.Qj})*cons.bodyj.pdot + ...
                        2*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.Bmatrix(cons.bodyi.pdot,cons.bodyi.point{cons.Pi})*cons.bodyi.pdot + ...
                       -2*utility.dijdot(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.dijdot(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj) + cons.fddot;
        end
        function phi_r = get.phi_r(cons) 
            % from ME751_f2016 slide 16 from lecture 9/28/16
            % One body can be the ground. In this case, the number of columns
            % in the Jacobian is half ? there are no partial derivatives with 
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
            % in the Jacobian is half ? there are no partial derivatives with 
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