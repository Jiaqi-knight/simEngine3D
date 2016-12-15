classdef rcj < handle
    % Filename: rcj.m
    % Author:   Samuel Acuña
    % Date:     17 Nov 2016
    % About:
    % This class handles instances of revolute-cylindrical composite joint (rcj) constraints. It Uses 
    % the r-p formulation (euler parameters). Removes 3 degrees of Freedom.
    %
    % (see haug book equation 9.4.32 on pgs 376-377)
    %
    % revolute cylindrical composite joint constraint: 
    % A coupler that is constrained to body i by a revolute joint about the
    % hi axis on body i and to body j through a cylindrical joint about the
    % hj axis. Vectors hi and hj are required to be orthogonal (dp1 constraint).
    
    properties
        system; % parent system3D object to which this constraint is a member. 
        rDOF = 3; % removes 3 degree of freedom
        bodyi;  % body i
        bodyj;  % body j
        hBari_tail; % ID number for point aBari_tail on body i, tail of hBari vector
        hBari_head; % ID number for point aBari_head on body i, head of hBari vector
        aBarj_tail; % ID number for point aBarj_tail on body j, tail of aBarj vector
        aBarj_head; % ID number for point aBarj_head on body j, head of aBarj vector
        bBarj_tail; % ID number for point bBarj_tail on body j, tail of bBarj vector
        bBarj_head; % ID number for point bBarj_head on body j, head of bBarj vector
        hBarj_tail; % ID number for point hBarj_tail on body j, tail of hBarj vector
        hBarj_head; % ID number for point hBarj_head on body j, head of hBarj vector
        Pi;         % ID number for point P on body i, tail of PiQj vector
        Qj;         % ID number for point Q on body j, head of PiQj vector
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
        function cons = rcj(system,bodyi,hBari_tail,hBari_head,PiID,bodyj,hBarj_tail,hBarj_head,aBarj_tail,aBarj_head,bBarj_tail,bBarj_head,QjID) %constructor function
            cons.system = system;
            cons.bodyi = bodyi;
            cons.hBari_tail = hBari_tail;
            cons.hBari_head = hBari_head;
            cons.Pi = PiID;
            cons.bodyj = bodyj;
            cons.hBarj_tail = hBarj_tail;
            cons.hBarj_head = hBarj_head;
            cons.aBarj_tail = aBarj_tail;
            cons.aBarj_head = aBarj_head;
            cons.bBarj_tail = bBarj_tail;
            cons.bBarj_head = bBarj_head;
            cons.Qj = QjID;
            
            % create cell array of all sub constraints
            
            %enforce hi and hj are orthogonal
            cons.subCons{1} = constraint.dp1(cons.system,cons.bodyi,cons.hBari_tail,cons.hBari_head,cons.bodyj,cons.hBarj_tail,cons.hBarj_head);
            
            %enforce hj passes through point Pi
            cons.subCons{2} = constraint.p2(cons.system,cons.bodyj,cons.aBarj_tail,cons.aBarj_head,cons.bBarj_tail,cons.bBarj_head,cons.Qj,cons.bodyi,cons.Pi);
            
            if abs(cons.phi) > 1e-4
                warning('Initial conditions for ''rcj'' are not consistent. But solution will converge so constraints are satisfied.')
            end
        end        
        function phi = get.phi(cons) % value of the expression of the constraint PHI^sj
            % phi : [3x1]
            phi = [cons.subCons{1}.phi; cons.subCons{2}.phi];
        end
        function nu = get.nu(cons) % right-hand side of the velocity equation
            % nu : [3x1]
            nu = [cons.subCons{1}.nu; cons.subCons{2}.nu];
        end
        function gammaHat = get.gammaHat(cons) % right-hand side of the acceleration equation, in r-p formulation
            % gammaHat : [3x1]
            gammaHat = [cons.subCons{1}.gammaHat; cons.subCons{2}.gammaHat];
        end
        function phi_r = get.phi_r(cons) % partial derivative of constraint with respect to r
            % phi_r : [3x6] normally, unless grounded, then [3x3]
            phi_r = [cons.subCons{1}.phi_r; cons.subCons{2}.phi_r];
        end
        function phi_p = get.phi_p(cons) % partial derivative of constraint with respect to p
            % phi_p : [3x8] normally, unless grounded, then [3x4]
            phi_p = [cons.subCons{1}.phi_p; cons.subCons{2}.phi_p];
        end
        function phiLambda_rr = phiLambda_rr(cons,lambda)
            % partial derivative of (phi_r*lambda) with respect to r (position)
            % following technique from get.phi_r and get.phi_p
            % phiLambda_rr : [18x6] normally, unless grounded, then [9x3]
            % inputs:
            %    lambda : [3x1] vector of lambda values, corresponding to
            %    the subconstraints

            phiLambda_rr = [cons.subCons{1}.phiLambda_rr(lambda(1)); cons.subCons{2}.phiLambda_rr(lambda(2:3))];
        end
        function phiLambda_rp = phiLambda_rp(cons,lambda)
            % partial derivative of (phi_r*lambda) with respect to p (orientation)
            % following technique from get.phi_r and get.phi_p
            % phiLambda_rp : [18x8] normally, unless grounded, then [9x4]
            % inputs:
            %    lambda : [3x1] vector of lambda values, corresponding to
            %    the subconstraints
            
            phiLambda_rp = [cons.subCons{1}.phiLambda_rp(lambda(1)); cons.subCons{2}.phiLambda_rp(lambda(2:3))];
        end
        function phiLambda_pr = phiLambda_pr(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to r (position)
            % following technique from get.phi_r and get.phi_p
            % phiLambda_pr : [24x6] normally, unless grounded, then [12x3]
            % inputs:
            %    lambda : [3x1] vector of lambda values, corresponding to
            %    the subconstraints
            
            phiLambda_pr = [cons.subCons{1}.phiLambda_pr(lambda(1)); cons.subCons{2}.phiLambda_pr(lambda(2:3))];
        end
        function phiLambda_pp = phiLambda_pp(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to p (orientation)
            % following technique from get.phi_r and get.phi_p
            % phiLambda_pp : [24x8] normally, unless grounded, then [12x4]
            % inputs:
            %    lambda : [3x1] vector of lambda values, corresponding to
            %    the subconstraints
            
            phiLambda_pp = [cons.subCons{1}.phiLambda_pp(lambda(1)); cons.subCons{2}.phiLambda_pp(lambda(2:3))];
        end
    end %methods
end %class