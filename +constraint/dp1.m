classdef dp1 < handle
    % Filename: dp1.m
    % Author:   Samuel Acuña
    % Date:     10 Oct 2016
    % About:
    % This class handles instances of DP1 constraints. It Uses the r-p
    % formulation (euler parameters). Removes 1 degree of Freedom.
    % DP1 constraint: the dot product between vector aBari and aBarj assume
    % a specified value f. f is often 0, meaning aBari and aBarj are
    % orthogonal vectors
    properties
        bodyi;  % body i
        bodyj;  % body j
        Pi;     % ID number for point P on body i, tail of aBari vector, body i RF
        Qi;     % ID number for point Q on body i, head of aBari vector, body i RF
        Pj;     % ID number for point P on body j, tail of aBarj vector, body j RF
        Qj;     % ID number for point Q on body j, head of aBarj vector, body j RF
        f;      % prescribed constraint, often 0, but can be a function of t
        fdot;   % derivative of f
        fddot;  % derivative of fdot
    end
    properties (Dependent)
        aBari;  % vector in body i RF
        aBarj;  % vector in body j RF
        phi;    % value of the expression of the constraint PHI^dp1
        nu;     % right-hand side of the velocity equation
        gammaHat;  % right-hand side of the acceleration equation, in r-p formulation
        phi_r;  % partial derivative of constraint with respect to r
        phi_p;  % partial derivative of constraint with respect to p
    end
    
    methods
        %constructor function
        function cons = dp1(bodyi,PiID,QiID,bodyj,PjID,QjID,f,fdot,fddot) %constructor function
            if ~exist('f','var') || isempty(f)
                f = 0; % prescribed constraint is 0, indicating vectors are orthogonal
            end
            if ~exist('fdot','var') || isempty(fdot)
                fdot = 0; % derivative of ft
            end
            if ~exist('fddot','var') || isempty(fddot)
                fddot = 0; % derivative of fdt
            end
            
            cons.bodyi = bodyi;
            cons.Pi = PiID;
            cons.Qi = QiID;
            cons.bodyj = bodyj;
            cons.Pj = PjID;
            cons.Qj = QjID;
            cons.f = f;
            cons.fdot = fdot;
            cons.fddot = fddot;
        end
        function aBari = get.aBari(cons) % vector in body i RF
            % aBari = Qi - Pi
            aBari = cons.bodyi.point{cons.Qi} - cons.bodyi.point{cons.Pi};
        end
        function aBarj = get.aBarj(cons) % vector in body j RF
            % aBarj = Qj - Pj
            aBarj = cons.bodyj.point{cons.Qj} - cons.bodyj.point{cons.Pj};
        end
        function phi = get.phi(cons) % value of the expression of the constraint PHI^dp1
            % from ME751_f2016 slide 11 from lecture 09/26/16
            % phi : [1x1]
            phi = (cons.bodyi.A*cons.aBari)'*(cons.bodyj.A*cons.aBarj) - cons.f;
        end
        function nu = get.nu(cons) % right-hand side of the velocity equation
            % from ME751_f2016 slide 12 from lecture 09/26/16
            % nu : [1x1]
            nu = cons.fdot;
        end
        function gammaHat = get.gammaHat(cons) % right-hand side of the acceleration equation, in r-p formulation
            % from ME751_f2016 slide 8 from lecture 10/7/16
            % gammaHat : [1x1]
            gammaHat = -(cons.bodyi.A*cons.aBari)'*utility.Bmatrix(cons.bodyj.pdot,cons.aBarj)*cons.bodyj.pdot + ...
                       -(cons.bodyj.A*cons.aBarj)'*utility.Bmatrix(cons.bodyi.pdot,cons.aBari)*cons.bodyi.pdot + ...
                       -2*(utility.Bmatrix(cons.bodyi.p,cons.aBari)*cons.bodyi.pdot)'*(utility.Bmatrix(cons.bodyj.p,cons.aBarj)*cons.bodyj.pdot)+cons.fddot;
        end
        function phi_r = get.phi_r(cons) 
            % from ME751_f2016 slide 13 from lecture 9/28/16
            % One body can be the ground. In this case, the number of columns
            % in the Jacobian is half ? there are no partial derivatives with 
            % respect to r or p. Thus, we must properly dimension the size of
            % the vectors/matrices because the grounded body does not have 
            % any generalized coordinates.
            % phi_r : [1x6] normally, unless grounded, then [1x3]
            
            phi_ri = zeros(1,3);
            phi_rj = zeros(1,3);
            
            if cons.bodyi.isGround
                phi_r = [phi_rj];
            elseif cons.bodyj.isGround
                phi_r = [phi_ri];
            else 
                phi_r = [phi_ri, phi_rj]; 
            end
        end
        function phi_p = get.phi_p(cons)
            % from ME751_f2016 slide 13 from lecture 9/28/16
            % One body can be the ground. In this case, the number of columns
            % in the Jacobian is half ? there are no partial derivatives with 
            % respect to r or p. Thus, we must properly dimension the size of
            % the vectors/matrices because the grounded body does not have 
            % any generalized coordinates.
            % phi_p : [1x8] normally, unless grounded, then [1x4]
            
            phi_pi = (cons.bodyj.A*cons.aBarj)'*(utility.Bmatrix(cons.bodyi.p,cons.aBari));
            phi_pj = (cons.bodyi.A*cons.aBari)'*(utility.Bmatrix(cons.bodyj.p,cons.aBarj));
            
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