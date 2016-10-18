classdef dp2 < handle
    % Filename: dp2.m
    % Author:   Samuel Acuña
    % Date:     17 Oct 2016
    % About:
    % This class handles instances of DP2 constraints. It Uses the r-p
    % formulation (euler parameters). Removes 1 degree of Freedom.
    % DP2 constraint: the dot product between vector aBari and a vector
    % PiQj from body i to body j assume a specified value f. f is often 0,
    % meaning aBari and PiQj are orthogonal vectors.
    properties
        bodyi;  % body i
        bodyj;  % body j
        aBari_tail; % ID number for point aBari_tail on body i, tail of aBari vector
        aBari_head; % ID number for point aBari_head on body i, head of aBari vector
        Pi;     % ID number for point P on body i, tail of PiQj vector
        Qj;     % ID number for point Q on body j, head of PiQj vector
        f;      % prescribed constraint, often 0, but can be a function of t
        fdot;   % derivative of f
        fddot;  % derivative of fdot
    end
    properties (Dependent)
        aBari;  % vector in body i RF
        PiQj;   % vector in GLOBAL RF
        phi;    % value of the expression of the constraint PHI^dp1
        nu;     % right-hand side of the velocity equation
        gammaHat;  % right-hand side of the acceleration equation, in r-p formulation
        phi_r;  % partial derivative of constraint with respect to r
        phi_p;  % partial derivative of constraint with respect to p
    end
    
    methods
        %constructor function
        function cons = dp2(bodyi,aBari_tailID,aBari_headID,PiID,bodyj,QjID,f,fdot,fddot) %constructor function
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
            cons.aBari_tail = aBari_tailID;
            cons.aBari_head = aBari_headID;
            cons.Pi = PiID;
            cons.bodyj = bodyj;
            cons.Qj = QjID;
            cons.f = f;
            cons.fdot = fdot;
            cons.fddot = fddot;
            
            if cons.phi ~= 0
                warning('Initial conditions for ''dp2'' are not consistent. But solution will converge so constraints are satisfied.')
            end
        end
        function aBari = get.aBari(cons) % vector in body i RF
            % aBari = aBari_head - aBari_tail
            aBari = cons.bodyi.point{cons.aBari_head} - cons.bodyi.point{cons.aBari_tail};
        end
        function PiQj = get.PiQj(cons) % vector in GLOBAL RF
            % PiQj = vector from Pi to Qj
            PiQj = utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj);
        end
        function phi = get.phi(cons) % value of the expression of the constraint PHI^dp1
            % from ME751_f2016 slide 14 from lecture 09/26/16
            % phi : [1x1]
            phi = (cons.bodyi.A*cons.aBari)'*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj) - cons.f;
        end
        function nu = get.nu(cons) % right-hand side of the velocity equation
            % from ME751_f2016 slide 12 from lecture 09/26/16
            % nu : [1x1]
            nu = cons.fdot;
        end
        function gammaHat = get.gammaHat(cons) % right-hand side of the acceleration equation, in r-p formulation
            % from ME751_f2016 slide 8 from lecture 10/7/16
            % gammaHat : [1x1]
            gammaHat = -(cons.bodyi.A*cons.aBari)'*utility.Bmatrix(cons.bodyj.pdot,cons.bodyj.point{cons.Qj})*cons.bodyj.pdot + ...
                        (cons.bodyi.A*cons.aBari)'*utility.Bmatrix(cons.bodyi.pdot,cons.bodyi.point{cons.Pi})*cons.bodyi.pdot + ...
                       -(utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj))'*utility.Bmatrix(cons.bodyi.pdot,cons.aBari)*cons.bodyi.pdot + ...
                       -2*(utility.Bmatrix(cons.bodyi.p,cons.aBari)*cons.bodyi.pdot)'*utility.dijdot(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj) + cons.fddot;
        end
        function phi_r = get.phi_r(cons) 
            % from ME751_f2016 slide 15 from lecture 9/28/16
            % One body can be the ground. In this case, the number of columns
            % in the Jacobian is half. there are no partial derivatives with 
            % respect to r or p. Thus, we must properly dimension the size of
            % the vectors/matrices because the grounded body does not have 
            % any generalized coordinates.
            % phi_r : [1x6] normally, unless grounded, then [1x3]
            
            phi_ri = -(cons.bodyi.A*cons.aBari)';
            phi_rj =  (cons.bodyi.A*cons.aBari)';
            
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
            % in the Jacobian is half. there are no partial derivatives with 
            % respect to r or p. Thus, we must properly dimension the size of
            % the vectors/matrices because the grounded body does not have 
            % any generalized coordinates.
            % phi_p : [1x8] normally, unless grounded, then [1x4]
            
            phi_pi = (utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj))'*utility.Bmatrix(cons.bodyi.p,cons.aBari) + ...
                    -(cons.bodyi.A*cons.aBari)'*utility.Bmatrix(cons.bodyi.p,cons.bodyi.point{cons.Pi});
            phi_pj = (cons.bodyi.A*cons.aBari)'*utility.Bmatrix(cons.bodyj.p,cons.bodyj.point{cons.Qj});
            
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