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
        phi_r;  % partial derivative of constraint with respect to r
        phi_p;  % partial derivative of constraint with respect to p
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
        function phi_r = get.phi_r(cons) 
            % from Haug eq 9.6.8
            % phi_r :  [1x3]
            phi_r = zeros(1,3);
        end
        function phi_p = get.phi_p(cons)
            % from Haug eq 9.6.8
            % phi_p : [1x4]
            phi_p = 2*cons.bodyi.p';
        end
    end
    
end