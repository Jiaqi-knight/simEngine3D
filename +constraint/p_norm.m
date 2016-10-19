classdef p_norm < handle
    % Filename: p_norm.m
    % Author:   Samuel Acuña
    % Date:     10 Oct 2016
    %
    % About:
    % This class handles instances of the Euler Parameter (p) Normalization
    % Constraint (p_norm). removes 1 degree of freedom.
    %
    % p_norm constraint: 
    % for a given body i, the dot product between transposed vector p and p
    % is equal to 1. or, p'*p-1 = 0
    
    properties
        rDOF = 1; % removes 1 degree of freedom
        bodyi;  % body i
    end
    properties (Dependent)
        phi;    % value of the expression of the constraint PHI^p_norm
        nu;     % right-hand side of the velocity equation
        gammaHat;  % right-hand side of the acceleration equation, in r-p formulation
%         phi_r;  % partial derivative of constraint with respect to r
%         phi_p;  % partial derivative of constraint with respect to p
    end
    
    methods
        %constructor function
        function cons = p_norm(bodyi) %constructor function
            cons.bodyi = bodyi;
        end
        function phi = get.phi(cons) % value of the expression of the constraint PHI^p_norm
            % from ME751_f2016 slide 24 from lecture 9/28/16
            % phi : p'*p-1 = 0 = [1x1]
            phi = cons.bodyi.p'*cons.bodyi.p - 1;
        end
        function nu = get.nu(cons) % right-hand side of the velocity equation
            % from ME751_f2016 slide 25 from lecture 9/28/16
            % nu : [1x1]
            nu = 0;
        end
        function gammaHat = get.gammaHat(cons) % right-hand side of the acceleration equation, in r-p formulation
            % from ME751_f2016 slide 25 from lecture 9/28/16
            % gammaHat : [1x1]
            gammaHat = -2*cons.bodyi.pdot'*cons.bodyi.pdot;
        end
%         function phi_r = get.phi_r(cons) 
%             % UNSURE
%             % One body can be the ground. In this case, the number of columns
%             % in the Jacobian is half. there are no partial derivatives with 
%             % respect to r or p. Thus, we must properly dimension the size of
%             % the vectors/matrices because the grounded body does not have 
%             % any generalized coordinates.
%             % phi_r : ??  [1x6] normally, unless grounded, then [1x3]
%             
%             phi_ri = ?? zeros(1,3);
%             phi_rj = ?? zeros(1,3);
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
%             % phi_p : [1x8] normally, unless grounded, then [1x4]
%             
%             phi_pi = ?? 
%             phi_pj = ?? ;
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