classdef d < handle
    % Filename: d.m
    % Author:   Samuel Acuña
    % Date:     14 Oct 2016
    %
    % About:
    % This class handles instances of D constraints. It Uses the r-p
    % formulation (euler parameters). Removes 1 degree of Freedom.
    %
    % D constraint: 
    % the distance between point P on body i and point Q on
    % body j assume a specified value f, f > 0. 
    %
    % Note: 
    % I'm using the constraint formulation from Haug eq 9.4.8, not
    % the way it was presented in class (ME751_f2016 slide 17 from lecture
    % 09/26/16)
    
    properties
        system; % parent system3D object to which this constraint is a member. 
        rDOF = 1; % removes 1 degree of freedom
        bodyi;  % body i
        bodyj;  % body j
        Pi;     % ID number for point P on body i
        Qj;     % ID number for point Q on body j
        f;      % prescribed distance constraint, > 0, and can be a function of t.
        fdot;   % derivative of f
        fddot;  % derivative of fdot
    end
    properties (Dependent)
        PiQj;   % vector in GLOBAL RF form point P on body i to point Q on body j
        phi;    % value of the expression of the constraint PHI^d
        nu;     % right-hand side of the velocity equation
        gammaHat;  % right-hand side of the acceleration equation, in r-p formulation
        phi_r;  % partial derivative of constraint with respect to r
        phi_p;  % partial derivative of constraint with respect to p
    end
    
    methods
        %constructor function
        function cons = d(system,bodyi,PiID,bodyj,QjID,f,fdot,fddot,t) %constructor function
            if exist('f','var') && (f <=0)
                error('D constraint cannot have a prescribed distance <= 0');
            end
            if ~exist('f','var') || isempty(f)
                f = 1; % prescribed constraint is 1
            end
            if ~exist('fdot','var') || isempty(fdot)
                fdot = 0; % derivative of ft
            end
            if ~exist('fddot','var') || isempty(fddot)
                fddot = 0; % derivative of fdt
            end

            cons.system = system;
            cons.bodyi = bodyi;
            cons.Pi = PiID;
            cons.bodyj = bodyj;
            cons.Qj = QjID;
            cons.f = f;
            cons.fdot = fdot;
            cons.fddot = fddot;
            
            if abs(cons.phi) > 1e-4
                warning('Initial conditions for ''d'' are not consistent. But solution will converge so constraints are satisfied.')
            end
        end
        function PiQj = get.PiQj(cons) % vector in GLOBAL RF of the constrained distance
            % PiQj = vector from Pi to Qj
            PiQj = utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj);
        end
        function phi = get.phi(cons) % value of the expression of the constraint PHI^d
            % from Haug 9.4.8
            % phi : [1x1]
            
            % see if constraint is a function
            if isa(cons.f, 'function_handle')
                fVal = cons.f(cons.system.time); %evaluate at system time
            else % constant value
                fVal = cons.f;
            end
            
            phi = utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)-(fVal)^2;
        end
        function nu = get.nu(cons) % right-hand side of the velocity equation
            % time derivative of f^2
            % nu : [1x1]
            
            if isa(cons.fdot, 'function_handle')
                fVal = cons.f(cons.t); %evaluate at time t
                fdotVal = cons.fdot(cons.system.time); %evaluate at system time
            else % constant value
                fVal = cons.f;
                fdotVal = cons.fdot;
            end
            
            nu = 2*fVal*fdotVal;
        end
        function gammaHat = get.gammaHat(cons) % right-hand side of the acceleration equation, in r-p formulation
            % from ME751_f2016 slide 8 from lecture 10/7/16, 
            % also, second time derivative of f^2
            % gammaHat : [1x1]
            
            % see if constraint is a function
            if isa(cons.fddot, 'function_handle')
                fVal = cons.f(cons.system.time); %evaluate at system time
                fdotVal = cons.fdot(cons.system.time); %evaluate at system time
                fddotVal = cons.fddot(cons.system.time); %evaluate at system time
            else % constant value
                fVal = cons.f;
                fdotVal = cons.fdot;
                fddotVal = cons.fddot;
            end
            
            
            gammaHat = -2*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.Bmatrix(cons.bodyj.pdot,cons.bodyj.point{cons.Qj})*cons.bodyj.pdot + ...
                        2*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.Bmatrix(cons.bodyi.pdot,cons.bodyi.point{cons.Pi})*cons.bodyi.pdot + ...
                       -2*utility.dijdot(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.dijdot(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj) + ...
                        2*fdotVal*fdotVal + 2*fVal*fddotVal;
        end
        function phi_r = get.phi_r(cons) 
            % from ME751_f2016 slide 16 from lecture 9/28/16
            % One body can be the ground. In this case, the number of columns
            % in the Jacobian is half. there are no partial derivatives with 
            % respect to r or p. Thus, we must properly dimension the size of
            % the vectors/matrices because the grounded body does not have 
            % any generalized coordinates.
            % phi_r : [1x6] normally, unless grounded, then [1x3]
            
            phi_ri = -2*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)';
            phi_rj =  2*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)';
            
            if cons.bodyi.isGround
                phi_r = [phi_rj];
            elseif cons.bodyj.isGround
                phi_r = [phi_ri];
            else 
                phi_r = [phi_ri, phi_rj]; 
            end
        end
        function phi_p = get.phi_p(cons)
            % from ME751_f2016 slide 17 from lecture 9/28/16
            % One body can be the ground. In this case, the number of columns
            % in the Jacobian is half. there are no partial derivatives with 
            % respect to r or p. Thus, we must properly dimension the size of
            % the vectors/matrices because the grounded body does not have 
            % any generalized coordinates.
            % phi_p : [1x8] normally, unless grounded, then [1x4]
            
            phi_pi = -2*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.Bmatrix(cons.bodyi.p,cons.bodyi.point{cons.Pi});
            phi_pj =  2*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj)'*utility.Bmatrix(cons.bodyj.p,cons.bodyj.point{cons.Qj});
            
            if cons.bodyi.isGround
                phi_p = [phi_pj];
            elseif cons.bodyj.isGround
                phi_p = [phi_pi];
            else 
                phi_p = [phi_pi, phi_pj]; 
            end
        end
        function  phiLambda_rr = phiLambda_rr(cons,lambda)
            % partial derivative of (phi_r*lambda) with respect to r (position)
            % from ME751_f2016 slide 35,38 from lecture 10/17/16
            % phiLambda_rr : [6x6] normally, unless grounded, then [3x3]
            % inputs:
            %    lambda : scalar value
            
            if cons.bodyi.isGround
                phiLambda_rr = 2*lambda*eye(3); % [3x3]
            elseif cons.bodyj.isGround
                phiLambda_rr = 2*lambda*eye(3); % [3x3]
            else % [6x6]
                phiLambda_rr = 2*lambda*[eye(3) -eye(3);
                                        -eye(3)  eye(3)];
            end
        end
        function  phiLambda_rp = phiLambda_rp(cons,lambda)
            % partial derivative of (phi_r*lambda) with respect to p (orientation)
            %from ME751_f2016 slide 35,38 from lecture 10/17/16
            % phiLambda_rp : [6x8] normally, unless grounded, then [3x4]
            % inputs:
            %    lambda : scalar value
            
            phiLambda_rpii =  utility.Bmatrix(cons.bodyi.p,cons.bodyi.point{cons.Pi});
            phiLambda_rpjj =  utility.Bmatrix(cons.bodyj.p,cons.bodyj.point{cons.Qj});
            phiLambda_rpij = -utility.Bmatrix(cons.bodyj.p,cons.bodyj.point{cons.Qj});
            phiLambda_rpji = -utility.Bmatrix(cons.bodyi.p,cons.bodyi.point{cons.Pi});
            
            if cons.bodyi.isGround
                phiLambda_rp = 2*lambda*phiLambda_rpjj; % [3x4]
            elseif cons.bodyj.isGround
                phiLambda_rp = 2*lambda*phiLambda_rpii; % [3x4]
            else % [6x8]
                phiLambda_rp = 2*lambda*[phiLambda_rpii  phiLambda_rpij;
                                         phiLambda_rpji  phiLambda_rpjj];
            end
        end
        function  phiLambda_pr = phiLambda_pr(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to r (position)
            %from ME751_f2016 slide 35,38 from lecture 10/17/16
            % phiLambda_pr : [8x6] normally, unless grounded, then [4x3]
            % inputs:
            %    lambda : scalar value
            
            phiLambda_prii =  utility.Bmatrix(cons.bodyi.p,cons.bodyi.point{cons.Pi})';
            phiLambda_prjj =  utility.Bmatrix(cons.bodyj.p,cons.bodyj.point{cons.Qj});
            phiLambda_prij = -utility.Bmatrix(cons.bodyi.p,cons.bodyi.point{cons.Pi})';
            phiLambda_prji = -utility.Bmatrix(cons.bodyj.p,cons.bodyj.point{cons.Qj});
            
            if cons.bodyi.isGround
                phiLambda_pr = 2*lambda*phiLambda_prjj; % [3x4]
            elseif cons.bodyj.isGround
                phiLambda_pr = 2*lambda*phiLambda_prii; % [3x4]
            else % [6x8]
                phiLambda_pr = 2*lambda*[phiLambda_prii  phiLambda_prij;
                                         phiLambda_prji  phiLambda_prjj];
            end
        end
        function  phiLambda_pp = phiLambda_pp(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to p (orientation)
            %from ME751_f2016 slide 35,38 from lecture 10/17/16
            % will be [4x4] or [8x8], depends on if there is a grounded body
            % inputs:
            %    lambda : scalar value
            
            Dij = utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj);
            phiLambda_ppii = -utility.Kmatrix(cons.bodyi.point{cons.Pi},Dij) + ...
                              utility.Bmatrix(cons.bodyi.p,cons.bodyi.point{cons.Pi})'*utility.Bmatrix(cons.bodyi.p,cons.bodyi.point{cons.Pi});
            phiLambda_ppjj =  utility.Kmatrix(cons.bodyj.point{cons.Qj},Dij) + ...
                              utility.Bmatrix(cons.bodyj.p,cons.bodyj.point{cons.Qj})'*utility.Bmatrix(cons.bodyj.p,cons.bodyj.point{cons.Qj});
            phiLambda_ppij = -utility.Bmatrix(cons.bodyi.p,cons.bodyi.point{cons.Pi})'*utility.Bmatrix(cons.bodyj.p,cons.bodyj.point{cons.Qj});
            phiLambda_ppji = -utility.Bmatrix(cons.bodyj.p,cons.bodyj.point{cons.Qj})'*utility.Bmatrix(cons.bodyi.p,cons.bodyi.point{cons.Pi});
            
            if cons.bodyi.isGround
                phiLambda_pp = 2*lambda*phiLambda_prjj; % [3x4]
            elseif cons.bodyj.isGround
                phiLambda_pp = 2*lambda*phiLambda_prii; % [3x4]
            else % [6x8]
                phiLambda_pp = 2*lambda*[phiLambda_ppii  phiLambda_ppij;
                                         phiLambda_ppji  phiLambda_ppjj];
            end
        end
    end % methods
end % class