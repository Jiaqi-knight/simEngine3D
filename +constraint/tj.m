classdef tj < handle
    % Filename: tj.m
    % Author:   Samuel Acuña
    % Date:     17 Nov 2016
    % About:
    % This class handles instances of translational joint (tj) constraints. It Uses 
    % the r-p formulation (euler parameters). Removes 5 degrees of Freedom.
    %
    % translational joint constraint: 
    % gemoetrically similar to a cylindrical joint, with the caveat that
    % that the roational degree of freedom is supressed. If body j is
    % fixed, body i can slide up and down. Therefore, the joint allows for
    % 1 degree of freedom of relative motion.
    
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
        aBarj_tail; % ID number for point aBarj_tail on body j, tail of aBarj vector
        aBarj_head; % ID number for point aBarj_head on body j, head of aBarj vector
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
        function cons = tj(system,bodyi,PiID,aBari_tail,aBari_head,bBari_tail,bBari_head,...
                           bodyj,QjID,aBarj_tail,aBarj_head,cBarj_tail,cBarj_head) %constructor function
            cons.system = system;
            cons.bodyi = bodyi;
            cons.Pi = PiID;
            cons.aBari_tail = aBari_tail;
            cons.aBari_head = aBari_head;
            cons.bBari_tail = bBari_tail;
            cons.bBari_head = bBari_head;
            cons.bodyj = bodyj;
            cons.Qj = QjID;
            cons.aBarj_tail = aBarj_tail;
            cons.aBarj_head = aBarj_head;
            cons.cBarj_tail = cBarj_tail;
            cons.cBarj_head = cBarj_head;
            
            % create cell array of all sub constraints
            % from ME751_f2016 slide 31 from lecture 09/26/16
            cons.subCons{1} = constraint.cj(cons.system,cons.bodyi,cons.Pi,cons.aBari_tail,cons.aBari_head,cons.bBari_tail,cons.bBari_head,...
                           cons.bodyj,cons.Qj,cons.cBarj_tail,cons.cBarj_head);
            cons.subCons{2} = constraint.dp1(cons.system,cons.bodyi,cons.aBari_tail,cons.aBari_head,...
                cons.bodyj,cons.aBarj_tail,cons.aBarj_head);
            
            if abs(cons.phi) > 1e-4
                warning('Initial conditions for ''tj'' are not consistent. But solution will converge so constraints are satisfied.')
            end
        end        
        function phi = get.phi(cons) % value of the expression of the constraint PHI^rj
            % from ME751_f2016 slide 31 from lecture 09/26/16
            % phi : [5x1]
            phi = [cons.subCons{1}.phi; cons.subCons{2}.phi];
        end
        function nu = get.nu(cons) % right-hand side of the velocity equation
            % from ME751_f2016 slide 31 from lecture 09/26/16
            % nu : [5x1]
            nu = [cons.subCons{1}.nu; cons.subCons{2}.nu];
        end
        function gammaHat = get.gammaHat(cons) % right-hand side of the acceleration equation, in r-p formulation
            % from ME751_f2016 slide 31 from lecture 09/26/16
            % gammaHat : [5x1]
            gammaHat = [cons.subCons{1}.gammaHat; cons.subCons{2}.gammaHat];
        end
        function phi_r = get.phi_r(cons) % partial derivative of constraint with respect to r
            % from ME751_f2016 slide 31 from lecture 09/26/16
            % phi_r : [5x6] normally, unless grounded, then [5x3]
            phi_r = [cons.subCons{1}.phi_r; cons.subCons{2}.phi_r];
        end
        function phi_p = get.phi_p(cons) % partial derivative of constraint with respect to p
            % from ME751_f2016 slide 31 from lecture 09/26/16
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
            
            phiLambda_rr = [cons.subCons{1}.phiLambda_rr(lambda(1:4)); cons.subCons{2}.phiLambda_rr(lambda(5))];
        end
        function phiLambda_rp = phiLambda_rp(cons,lambda)
            % partial derivative of (phi_r*lambda) with respect to p (orientation)
            % following technique from get.phi_r and get.phi_p
            % phiLambda_rp : [30x8] normally, unless grounded, then [15x4]
            % inputs:
            %    lambda : [5x1] vector of lambda values, corresponding to
            %    the subconstraints
            
            phiLambda_rp = [cons.subCons{1}.phiLambda_rp(lambda(1:4)); cons.subCons{2}.phiLambda_rp(lambda(5))];
        end    
        function phiLambda_pr = phiLambda_pr(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to r (position)
            % following technique from get.phi_r and get.phi_p
            % phiLambda_pr : [40x6] normally, unless grounded, then [20x3]
            % inputs:
            %    lambda : [5x1] vector of lambda values, corresponding to
            %    the subconstraints
            
            phiLambda_pr = [cons.subCons{1}.phiLambda_pr(lambda(1:4)); cons.subCons{2}.phiLambda_pr(lambda(5))];
        end
        function phiLambda_pp = phiLambda_pp(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to p (orientation)
            % following technique from get.phi_r and get.phi_p
            % phiLambda_pp : [40x8] normally, unless grounded, then [20x4]
            % inputs:
            %    lambda : [5x1] vector of lambda values, corresponding to
            %    the subconstraints
            
            phiLambda_pp = [cons.subCons{1}.phiLambda_pp(lambda(1:4)); cons.subCons{2}.phiLambda_pp(lambda(5))];
        end 
    end %methods
end %class