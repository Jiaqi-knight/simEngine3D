classdef rj < handle
    % Filename: rj.m
    % Author:   Samuel Acuña
    % Date:     19 Oct 2016
    % About:
    % This class handles instances of revolute joint (rj) constraints. It Uses 
    % the r-p formulation (euler parameters). Removes 5 degrees of Freedom.
    %
    % revolute joint constraint: 
    % if body i is fixed, body j can rotate around teh hinge axis. Body j
    % has one degree of freedome relaive to body i. This GCon is built
    % by combining a spherical joint (sj) with a perpendicular1 (p1)
    % constraint.
    
    properties
        system; % parent system3D object to which this constraint is a member. 
        rDOF = 5; % removes 5 degree of freedom
        bodyi;  % body i
        bodyj;  % body j
        Pi;     % ID number for point P on body i
        Qj;     % ID number for point Q on body j
        aBari_tail; % ID number for point aBari_tail on body i, tail of aBari vector
        aBari_head; % ID number for point aBari_head on body i, head of aBari vector
        bBari_tail; % ID number for point bBari_tail on body i, tail of bBari vector
        bBari_head; % ID number for point bBari_head on body i, head of bBari vector
        cBarj_tail; % ID number for point cBarj_tail on body j, tail of cBarj vector
        cBarj_head; % ID number for point cBarj_head on body j, head of cBarj vector
        subCons; % cell array of sub-constraints
    end
    properties (Dependent)
        phi;    % value of the expression of the constraint PHI^rj
        nu;     % right-hand side of the velocity equation
        gammaHat;  % right-hand side of the acceleration equation, in r-p formulation
        phi_r;  % partial derivative of constraint with respect to r
        phi_p;  % partial derivative of constraint with respect to p
    end
    
    methods
        %constructor function
        function cons = rj(system,bodyi,PiID,aBari_tail,aBari_head,bBari_tail,bBari_head,...
                           bodyj,QjID,cBarj_tail,cBarj_head) %constructor function
            cons.system = system;
            cons.bodyi = bodyi;
            cons.Pi = PiID;
            cons.aBari_tail = aBari_tail;
            cons.aBari_head = aBari_head;
            cons.bBari_tail = bBari_tail;
            cons.bBari_head = bBari_head;
            cons.bodyj = bodyj;
            cons.Qj = QjID;
            cons.cBarj_tail = cBarj_tail;
            cons.cBarj_head = cBarj_head;
            
            % create cell array of all sub constraints
            % from ME751_f2016 slide 33 from lecture 09/26/16
            cons.subCons{1} = constraint.sj(cons.system,cons.bodyi,cons.Pi,cons.bodyj,cons.Qj);
            cons.subCons{2} = constraint.p1(cons.system,cons.bodyi,cons.aBari_tail,cons.aBari_head,...
                cons.bBari_tail,cons.bBari_head,cons.bodyj,cons.cBarj_tail,cons.cBarj_head);
            
            if abs(cons.phi) > 1e-4
                warning('Initial conditions for ''sj'' are not consistent. But solution will converge so constraints are satisfied.')
            end
        end        
        function phi = get.phi(cons) % value of the expression of the constraint PHI^rj
            % from ME751_f2016 slide 33 from lecture 09/26/16
            % phi : [5x1]
            phi = [cons.subCons{1}.phi; cons.subCons{2}.phi];
        end
        function nu = get.nu(cons) % right-hand side of the velocity equation
            % from ME751_f2016 slide 33 from lecture 09/26/16
            % nu : [5x1]
            nu = [cons.subCons{1}.nu; cons.subCons{2}.nu];
        end
        function gammaHat = get.gammaHat(cons) % right-hand side of the acceleration equation, in r-p formulation
            % from ME751_f2016 slide 33 from lecture 09/26/16
            % gammaHat : [5x1]
            gammaHat = [cons.subCons{1}.gammaHat; cons.subCons{2}.gammaHat];
        end
        function phi_r = get.phi_r(cons) % partial derivative of constraint with respect to r
            % from ME751_f2016 slide 33 from lecture 09/26/16
            % phi_r : [5x6] normally, unless grounded, then [5x3]
            phi_r = [cons.subCons{1}.phi_r; cons.subCons{2}.phi_r];
        end
        function phi_p = get.phi_p(cons) % partial derivative of constraint with respect to p
            % from ME751_f2016 slide 33 from lecture 09/26/16
            % phi_p : [5x8] normally, unless grounded, then [5x4]
            phi_p = [cons.subCons{1}.phi_p; cons.subCons{2}.phi_p];
        end
        function phiLambda_rr = phiLambda_rr(cons,lambda)
            % partial derivative of (phi_r*lambda) with respect to r (position)
            % following technique from get.phi_r and get.phi_p
            % phiLambda_rr : [30x6] normally, unless grounded, then [15x3]
            % inputs:
            %    lambda : [5x1] vector of lambda values, corresponding to
            %    the subconstraints
            
            phiLambda_rr = [cons.subCons{1}.phiLambda_rr(lambda(1:3)); cons.subCons{2}.phiLambda_rr(lambda(4:5))];
        end
        function phiLambda_rp = phiLambda_rp(cons,lambda)
            % partial derivative of (phi_r*lambda) with respect to p (orientation)
            % following technique from get.phi_r and get.phi_p
            % phiLambda_rp : [30x8] normally, unless grounded, then [15x4]
            % inputs:
            %    lambda : [5x1] vector of lambda values, corresponding to
            %    the subconstraints
            
            phiLambda_rp = [cons.subCons{1}.phiLambda_rp(lambda(1:3)); cons.subCons{2}.phiLambda_rp(lambda(4:5))];
        end    
        function phiLambda_pr = phiLambda_pr(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to r (position)
            % following technique from get.phi_r and get.phi_p
            % phiLambda_pr : [40x6] normally, unless grounded, then [20x3]
            % inputs:
            %    lambda : [5x1] vector of lambda values, corresponding to
            %    the subconstraints
            
            phiLambda_pr = [cons.subCons{1}.phiLambda_pr(lambda(1:3)); cons.subCons{2}.phiLambda_pr(lambda(4:5))];
        end
        function phiLambda_pp = phiLambda_pp(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to p (orientation)
            % following technique from get.phi_r and get.phi_p
            % phiLambda_pp : [40x8] normally, unless grounded, then [20x4]
            % inputs:
            %    lambda : [5x1] vector of lambda values, corresponding to
            %    the subconstraints
            
            phiLambda_pp = [cons.subCons{1}.phiLambda_pp(lambda(1:3)); cons.subCons{2}.phiLambda_pp(lambda(4:5))];
        end 
    end %methods
end %class