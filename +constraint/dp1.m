classdef dp1 < handle
    % Filename: dp1.m
% Author:   Samuel Acuña
% Date:     10 Oct 2016
% About:    
% This class handles instances of DP1 constraints. It Uses the r-p 
% formulation (euler parameters). Removes 1 degree of Freedom.

% Outputs (any or all):
%   X the value of the expression of the constraint "PHI_dp1"
%   X the right-hand side of the velocity equation "nu"
%   - the right-hand side of the acceleration equation "gamma"
%   - the expression of the partial derivatives "PHI_r" and "PHI_p"
    properties
        bodyi;  % body i
        bodyj;  % body j
        Pi;     % ID number for point P on body i, tail of aBari vector, body i RF
        Qi;     % ID number for point Q on body i, head of aBari vector, body i RF
        Pj;     % ID number for point P on body j, tail of aBarj vector, body j RF
        Qj;     % ID number for point Q on body j, head of aBarj vector, body j RF
        f;     % prescribed constraint, often 0, but can be a function of t
        fdot;    % derivative of ft
        fddot;   % derivative of fddt
    end
    properties (Dependent)
        aBari;  % vector in body i RF
        aBarj;  % vector in body j RF
        phi;    % value of the expression of the constraint PHI^dp1
        nu;     % right-hand side of the velocity equation
        gamma;  % right-hand side of the acceleration equation, in r-p formulation
    end
    
    methods
        %constructor function
        function cons = dp1(bodyi,PiID,QiID,bodyj,PjID,QjID,f,fdot,fddot) %constructor function
            if nargin <= 6
                f = 0; % prescribed constraint is 0, indicating vectors are orthogonal
            end
            if nargin <= 7
                fdot = 0; % derivative of ft
            end
            if nargin <=8
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
            Ai = utility.p2A(cons.bodyi.p);
            Aj = utility.p2A(cons.bodyj.p);
            phi = (Ai*cons.aBari)'*(Aj*cons.aBarj) - cons.f;
        end
        function nu = get.nu(cons) % right-hand side of the velocity equation
            nu = cons.fdot;
        end
        function gamma = get.gamma(cons) % right-hand side of the acceleration equation, in r-p formulation
            
        end
    end
    
end