classdef p2 < handle
    % Filename: p2.m
    % Author:   Samuel Acu�a
    % Date:     19 Oct 2016
    % About:
    % This class handles instances of perpendicular1 (p2) constraints. It Uses 
    % the r-p formulation (euler parameters). Removes 2 degrees of Freedom.
    %
    % perpendicular2 constraint: 
    % vector PiQj from body i to body j remains perpendicular to a plane.
    % This plane is specified by two noncolinear vectors aBari and bBari 
    % that are contained in that plane. This GCon is built using dp2 twice.
    
    properties
        system;     % parent system3D object to which this constraint is a member. 
        rDOF = 2;   % removes 2 degree of freedom
        bodyi;      % body i
        bodyj;      % body j
        aBari_tail; % ID number for point aBari_tail on body i, tail of aBari vector
        aBari_head; % ID number for point aBari_head on body i, head of aBari vector
        bBari_tail; % ID number for point bBari_tail on body i, tail of bBari vector
        bBari_head; % ID number for point bBari_head on body i, head of bBari vector
        Pi;         % ID number for point P on body i, tail of PiQj vector
        Qj;         % ID number for point Q on body j, head of PiQj vector
        subCons;    % cell array of sub-constraints
    end
    properties (Dependent)
        aBari;      % vector in body i RF
        bBari;      % vector in body i RF
        PiQj;       % vector in GLOBAL RF
        phi;        % value of the expression of the constraint PHI^p2
        nu;         % right-hand side of the velocity equation
        gammaHat;   % right-hand side of the acceleration equation, in r-p formulation
        phi_r;      % partial derivative of constraint with respect to r
        phi_p;      % partial derivative of constraint with respect to p
    end
    
    methods
        %constructor function
        function cons = p2(system,bodyi,aBari_tail,aBari_head,bBari_tail,bBari_head,PiID,bodyj,QjID) %constructor function
            cons.system = system;
            cons.bodyi = bodyi;
            cons.aBari_tail = aBari_tail;
            cons.aBari_head = aBari_head;
            cons.bBari_tail = bBari_tail;
            cons.bBari_head = bBari_head;
            cons.Pi = PiID;
            cons.bodyj = bodyj;
            cons.Qj = QjID;
                            
            % create cell array of all sub constraints
            % from ME751_f2016 slide 25 from lecture 09/26/16
            cons.subCons{1} = constraint.dp2(cons.system,cons.bodyi,cons.aBari_tail,cons.aBari_head,cons.Pi,cons.bodyj,cons.Qj);
            cons.subCons{2} = constraint.dp2(cons.system,cons.bodyi,cons.bBari_tail,cons.bBari_head,cons.Pi,cons.bodyj,cons.Qj);
            
            if abs(cons.phi) > 1e-4
                warning('Initial conditions for ''p2'' are not consistent. But solution will converge so constraints are satisfied.')
            end
        end
        function aBari = get.aBari(cons) % vector in body i RF
            % aBari = aBari_head - aBari_tail
            aBari = cons.bodyi.point{cons.aBari_head} - cons.bodyi.point{cons.aBari_tail};
        end
        function bBari = get.bBari(cons) % vector in body i RF
            % bBari = bBari_head - bBari_tail
            bBari = cons.bodyi.point{cons.bBari_head} - cons.bodyi.point{cons.bBari_tail};
        end
        function PiQj = get.PiQj(cons) % vector in GLOBAL RF
            % PiQj = vector from Pi to Qj
            PiQj = utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj);
        end
        function phi = get.phi(cons) % value of the expression of the constraint PHI^p2
            % from ME751_f2016 slide 25 from lecture 09/26/16
            % phi : [2x1]
            phi = [cons.subCons{1}.phi; cons.subCons{2}.phi];
        end
        function nu = get.nu(cons) % right-hand side of the velocity equation
            % from ME751_f2016 slide 25 from lecture 09/26/16
            % nu : [2x1]
            nu = [cons.subCons{1}.nu; cons.subCons{2}.nu];
        end
        function gammaHat = get.gammaHat(cons) % right-hand side of the acceleration equation, in r-p formulation
            % from ME751_f2016 slide 25 from lecture 09/26/16
            % gammaHat : [2x1]
            gammaHat = [cons.subCons{1}.gammaHat; cons.subCons{2}.gammaHat];
        end
        function phi_r = get.phi_r(cons) % partial derivative of constraint with respect to r
            % from ME751_f2016 slide 25 from lecture 09/26/16
            % phi_r : [2x6] normally, unless grounded, then [2x3]
            phi_r = [cons.subCons{1}.phi_r; cons.subCons{2}.phi_r];
        end
        function phi_p = get.phi_p(cons) % partial derivative of constraint with respect to p
            % from ME751_f2016 slide 25 from lecture 09/26/16
            % phi_p : [2x8] normally, unless grounded, then [2x4]
            phi_p = [cons.subCons{1}.phi_p; cons.subCons{2}.phi_p];
        end
        function phiLambda_rr = phiLambda_rr(cons,lambda)
            % partial derivative of (phi_r*lambda) with respect to r (position)
            % following technique from get.phi_r and get.phi_p
            % phiLambda_rr : [12x6] normally, unless grounded, then [6x3]
            % inputs:
            %    lambda : [2x1] vector of lambda values, corresponding to
            %    the subconstraints
            
            phiLambda_rr = [cons.subCons{1}.phiLambda_rr(lambda(1)); cons.subCons{2}.phiLambda_rr(lambda(2))];
        end
        function phiLambda_rp = phiLambda_rp(cons,lambda)
            % partial derivative of (phi_r*lambda) with respect to p (orientation)
            % following technique from get.phi_r and get.phi_p
            % phiLambda_rp : [12x8] normally, unless grounded, then [6x4]
            % inputs:
            %    lambda : [2x1] vector of lambda values, corresponding to
            %    the subconstraints
            
            phiLambda_rp = [cons.subCons{1}.phiLambda_rp(lambda(1)); cons.subCons{2}.phiLambda_rp(lambda(2))];
        end    
        function phiLambda_pr = phiLambda_pr(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to r (position)
            % following technique from get.phi_r and get.phi_p
            % phiLambda_pr : [16x6] normally, unless grounded, then [8x3]
            % inputs:
            %    lambda : [2x1] vector of lambda values, corresponding to
            %    the subconstraints
            
            phiLambda_pr = [cons.subCons{1}.phiLambda_pr(lambda(1)); cons.subCons{2}.phiLambda_pr(lambda(2))];
        end
        function phiLambda_pp = phiLambda_pp(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to p (orientation)
            % following technique from get.phi_r and get.phi_p
            % phiLambda_pp : [16x8] normally, unless grounded, then [8x4]
            % inputs:
            %    lambda : [2x1] vector of lambda values, corresponding to
            %    the subconstraints
            
            phiLambda_pp = [cons.subCons{1}.phiLambda_pp(lambda(1)); cons.subCons{2}.phiLambda_pp(lambda(2))];
        end 
    end %methods
end %class