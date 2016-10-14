classdef cd < handle
    % Filename: cd.m
    % Author:   Samuel Acuña
    % Date:     13 Oct 2016
    % About:
    % This class handles instances of CD constraints. It Uses the r-p
    % formulation (euler parameters). Removes 1 degree of Freedom.
    % CD constraint: the difference between the x (or y or z) coordinate of
    % point P on body i and the x (or y or z) coordinate of point Q on body
    % j assume a specified value f.
    properties
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
        phi;    % value of the expression of the constraint PHI^dp1
        nu;     % right-hand side of the velocity equation
        gammaHat;  % right-hand side of the acceleration equation, in r-p formulation
        phi_r;  % partial derivative of constraint with respect to r
        phi_p;  % partial derivative of constraint with respect to p
    end
    
    methods
        %constructor function
        function cons = cd(cName,bodyi,PiID,bodyj,QjID,f,fdot,fddot) %constructor function
            if ~exist('f','var') || isempty(f)
                f = 0; % prescribed constraint is 0, indicating vectors are orthogonal
            end
            if ~exist('fdot','var') || isempty(fdot)
                fdot = 0; % derivative of ft
            end
            if ~exist('fddot','var') || isempty(fddot)
                fddot = 0; % derivative of fdt
            end
            
            cons.cName = cName;
            cons.bodyi = bodyi;
            cons.Pi = PiID;
            cons.bodyj = bodyj;
            cons.Qj = QjID;
            cons.f = f;
            cons.fdot = fdot;
            cons.fddot = fddot;
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
        function phi = get.phi(cons) % value of the expression of the constraint PHI^dp1
            % from ME751_f2016 slide 20 from lecture 09/26/16
            % phi : [1x1]
            phi = cons.c'*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)-cons.f;
        end
        function nu = get.nu(cons) % right-hand side of the velocity equation
            % from ME751_f2016 slide 20 from lecture 09/26/16
            % nu : [1x1]
            nu = cons.fdot;
        end
        function gammaHat = get.gammaHat(cons) % right-hand side of the acceleration equation, in r-p formulation
            % from ME751_f2016 slide 8 from lecture 10/7/16
            % gammaHat : [1x1]
            gammaHat =  cons.c'*utility.Bmatrix(cons.bodyi.pdot,cons.bodyi.point{cons.Pi})*cons.bodyi.pdot + ...
                       -cons.c'*utility.Bmatrix(cons.bodyj.pdot,cons.bodyj.point{cons.Qj})*cons.bodyj.pdot + cons.fddot;
        end
        function phi_r = get.phi_r(cons) 
            % from ME751_f2016 slide 17 from lecture 9/28/16
            % One body can be the ground. In this case, the number of columns
            % in the Jacobian is half ? there are no partial derivatives with 
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
            % in the Jacobian is half ? there are no partial derivatives with 
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
    end
    
end