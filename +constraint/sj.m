classdef sj < handle
    % Filename: sj.m
    % Author:   Samuel Acuña
    % Date:     19 Oct 2016
    % About:
    % This class handles instances of spherical joint (sj) constraints. It Uses 
    % the r-p formulation (euler parameters). Removes 3 degrees of Freedom.
    %
    % spherical joint constraint: 
    % the point P on body i and point Q on body j coincide at all times.
    % This GCon is built by using 3 cd constraints.
    
    properties
        rDOF = 3; % removes 3 degree of freedom
        bodyi;  % body i
        bodyj;  % body j
        Pi;     % ID number for point P on body i
        Qj;     % ID number for point Q on body j
        subCons; % cell array of sub-constraints
    end
    properties (Dependent)
        phi;    % value of the expression of the constraint PHI^sj
        nu;     % right-hand side of the velocity equation
        gammaHat;  % right-hand side of the acceleration equation, in r-p formulation
        phi_r;  % partial derivative of constraint with respect to r
        phi_p;  % partial derivative of constraint with respect to p
    end
    
    methods
        %constructor function
        function cons = sj(bodyi,PiID,bodyj,QjID) %constructor function
            cons.bodyi = bodyi;
            cons.Pi = PiID;
            cons.bodyj = bodyj;
            cons.Qj = QjID;
            
            % create cell array of all sub constraints
            % from ME751_f2016 slide 27 from lecture 09/26/16
            cons.subCons{1} = constraint.cd('x',cons.bodyi,cons.Pi,cons.bodyj,cons.Qj);
            cons.subCons{2} = constraint.cd('y',cons.bodyi,cons.Pi,cons.bodyj,cons.Qj);
            cons.subCons{3} = constraint.cd('z',cons.bodyi,cons.Pi,cons.bodyj,cons.Qj);
                
            if abs(cons.phi) > 1e-4
                warning('Initial conditions for ''sj'' are not consistent. But solution will converge so constraints are satisfied.')
            end
        end        
        function phi = get.phi(cons) % value of the expression of the constraint PHI^sj
            % from ME751_f2016 slide 27 from lecture 09/26/16
            % phi : [3x1]
            phi = [cons.subCons{1}.phi; cons.subCons{2}.phi; cons.subCons{3}.phi];
        end
        function nu = get.nu(cons) % right-hand side of the velocity equation
            % from ME751_f2016 slide 27 from lecture 09/26/16
            % nu : [3x1]
            nu = [cons.subCons{1}.nu; cons.subCons{2}.nu; cons.subCons{3}.nu];
        end
        function gammaHat = get.gammaHat(cons) % right-hand side of the acceleration equation, in r-p formulation
            % from ME751_f2016 slide 27 from lecture 09/26/16
            % gammaHat : [3x1]
            gammaHat = [cons.subCons{1}.gammaHat; cons.subCons{2}.gammaHat; cons.subCons{3}.gammaHat];
        end
        function phi_r = get.phi_r(cons) % partial derivative of constraint with respect to r
            % from ME751_f2016 slide 27 from lecture 09/26/16
            % phi_r : [3x6] normally, unless grounded, then [3x3]
            phi_r = [cons.subCons{1}.phi_r; cons.subCons{2}.phi_r; cons.subCons{3}.phi_r];
        end
        function phi_p = get.phi_p(cons) % partial derivative of constraint with respect to p
            % from ME751_f2016 slide 27 from lecture 09/26/16
            % phi_p : [3x8] normally, unless grounded, then [3x4]
            phi_p = [cons.subCons{1}.phi_p; cons.subCons{2}.phi_p; cons.subCons{3}.phi_p];
        end
    end
    
end