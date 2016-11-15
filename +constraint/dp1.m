classdef dp1 < handle
    % Filename: dp1.m
    % Author:   Samuel Acuña
    % Date:     10 Oct 2016
    %
    % About:
    % This class handles instances of DP1 constraints. It Uses the r-p
    % formulation (euler parameters). Removes 1 degree of Freedom.
    %
    % DP1 constraint: 
    % the dot product between vector aBari and aBarj assumes
    % a specified value f. f is often 0, meaning aBari and aBarj are
    % orthogonal vectors
    
    properties
        system; % parent system3D object to which this constraint is a member. 
        rDOF = 1; % removes 1 degree of freedom
        bodyi;  % body i
        bodyj;  % body j
        Pi;     % ID number for point P on body i, tail of aBari vector, body i RF
        Qi;     % ID number for point Q on body i, head of aBari vector, body i RF
        Pj;     % ID number for point P on body j, tail of aBarj vector, body j RF
        Qj;     % ID number for point Q on body j, head of aBarj vector, body j RF
        f;      % prescribed constraint, often 0, but can be a function of t
        fdot;   % derivative of f
        fddot;  % derivative of fdot
    end
    properties (Dependent)
        aBari;  % vector in body i RF
        aBarj;  % vector in body j RF
        phi;    % value of the expression of the constraint PHI^dp1
        nu;     % right-hand side of the velocity equation
        gammaHat;  % right-hand side of the acceleration equation, in r-p formulation
        phi_r;  % partial derivative of constraint with respect to r
        phi_p;  % partial derivative of constraint with respect to p
    end
    
    methods
        %constructor function
        function cons = dp1(system,bodyi,PiID,QiID,bodyj,PjID,QjID,f,fdot,fddot) %constructor function
            if ~exist('f','var') || isempty(f)
                f = 0; % prescribed constraint is 0, indicating vectors are orthogonal
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
            cons.Qi = QiID;
            cons.bodyj = bodyj;
            cons.Pj = PjID;
            cons.Qj = QjID;
            cons.f = f;
            cons.fdot = fdot;
            cons.fddot = fddot;
            
               
            if abs(cons.phi) > 1e-4
                warning('Initial conditions for ''dp1'' are not consistent. But solution will converge so constraints are satisfied.')
            end
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
            % from ME751_f2016 slide 11 from lecture 09/26/16
            % phi : [1x1]
            
            % see if constraint is a function
            if isa(cons.f, 'function_handle')
                fVal = cons.f(cons.system.time); %evaluate at system time
            else % constant value
                fVal = cons.f;
            end
            
            phi = (cons.bodyi.A*cons.aBari)'*(cons.bodyj.A*cons.aBarj) - fVal;
        end
        function nu = get.nu(cons) % right-hand side of the velocity equation
            % from ME751_f2016 slide 12 from lecture 09/26/16
            % nu : [1x1]
            
            % see if constraint is a function
            if isa(cons.fdot, 'function_handle')
                fdotVal = cons.fdot(cons.system.time); %evaluate at system time
            else % constant value
                fdotVal = cons.fdot;
            end
            
            nu = fdotVal;
        end
        function gammaHat = get.gammaHat(cons) % right-hand side of the acceleration equation, in r-p formulation
            % from ME751_f2016 slide 8 from lecture 10/7/16
            % gammaHat : [1x1]
            
            % see if constraint is a function
            if isa(cons.fddot, 'function_handle')
                fddotVal = cons.fddot(cons.system.time); %evaluate at system time
            else % constant value
                fddotVal = cons.fddot;
            end
            
            gammaHat = -(cons.bodyi.A*cons.aBari)'*utility.Bmatrix(cons.bodyj.pdot,cons.aBarj)*cons.bodyj.pdot + ...
                       -(cons.bodyj.A*cons.aBarj)'*utility.Bmatrix(cons.bodyi.pdot,cons.aBari)*cons.bodyi.pdot + ...
                       -2*(utility.Bmatrix(cons.bodyi.p,cons.aBari)*cons.bodyi.pdot)'*(utility.Bmatrix(cons.bodyj.p,cons.aBarj)*cons.bodyj.pdot) + ...
                       fddotVal;
        end
        function phi_r = get.phi_r(cons) % get jacobian with respect to position (r)
            % from ME751_f2016 slide 13 from lecture 9/28/16
            % One body can be the ground. In this case, the number of columns
            % in the Jacobian is half. there are no partial derivatives with 
            % respect to r or p. Thus, we must properly dimension the size of
            % the vectors/matrices because the grounded body does not have 
            % any generalized coordinates.
            % phi_r : [1x6] normally, unless grounded, then [1x3]
            
            phi_ri = zeros(1,3);
            phi_rj = zeros(1,3);
            
            if cons.bodyi.isGround
                phi_r = [phi_rj];
            elseif cons.bodyj.isGround
                phi_r = [phi_ri];
            else 
                phi_r = [phi_ri, phi_rj]; 
            end
        end 
        function phi_p = get.phi_p(cons) % get jacobian with respect to orientation (p)
            % from ME751_f2016 slide 13 from lecture 9/28/16
            % One body can be the ground. In this case, the number of columns
            % in the Jacobian is half. there are no partial derivatives with 
            % respect to r or p. Thus, we must properly dimension the size of
            % the vectors/matrices because the grounded body does not have 
            % any generalized coordinates.
            % phi_p : [1x8] normally, unless grounded, then [1x4]
            
            phi_pi = (cons.bodyj.A*cons.aBarj)'*(utility.Bmatrix(cons.bodyi.p,cons.aBari));
            phi_pj = (cons.bodyi.A*cons.aBari)'*(utility.Bmatrix(cons.bodyj.p,cons.aBarj));
            
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
            % from ME751_f2016 slide 35,36 from lecture 10/17/16
            % phiLambda_rr : [6x6] normally, unless grounded, then [3x3]
            % inputs:
            %    lambda : scalar value
            
            dim_r = length(cons.phi_r);
            phiLambda_rr = lambda*zeros(dim_r,dim_r); % [3x3] or [6x6]
        end
        function  phiLambda_rp = phiLambda_rp(cons,lambda)
            % partial derivative of (phi_r*lambda) with respect to p (orientation)
            %from ME751_f2016 slide 35,36 from lecture 10/17/16
            % phiLambda_rp : [6x8] normally, unless grounded, then [3x4]
            % inputs:
            %    lambda : scalar value
            
            dim_r = length(cons.phi_r);
            dim_p = length(cons.phi_p);
            phiLambda_rp = lambda*zeros(dim_r,dim_p); % [3x4] or [6x8]
        end
        function  phiLambda_pr = phiLambda_pr(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to r (position)
            %from ME751_f2016 slide 35,36 from lecture 10/17/16
            % phiLambda_pr : [8x6] normally, unless grounded, then [4x3]
            % inputs:
            %    lambda : scalar value
            
            dim_r = length(cons.phi_r);
            dim_p = length(cons.phi_p);
            phiLambda_pr = lambda*zeros(dim_p,dim_r); % [4x3] or [8x6]
        end
        function  phiLambda_pp = phiLambda_pp(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to p (orientation)
            %from ME751_f2016 slide 35,36 from lecture 10/17/16
            % will be [4x4] or [8x8], depends on if there is a grounded body
            % inputs:
            %    lambda : scalar value
            
            phiLambda_ppii = utility.Kmatrix(cons.aBari,cons.bodyj.A*cons.aBarj);
            phiLambda_ppjj = utility.Kmatrix(cons.aBarj,cons.bodyi.A*cons.aBari);
            phiLambda_ppij = utility.Bmatrix(cons.bodyi.p,cons.aBari)'*utility.Bmatrix(cons.bodyj.p,cons.aBarj);
            phiLambda_ppji = utility.Bmatrix(cons.bodyj.p,cons.aBarj)'*utility.Bmatrix(cons.bodyi.p,cons.aBari);
            
            if cons.bodyi.isGround
                phiLambda_pp = lambda*phiLambda_ppjj; % [4x4]
            elseif cons.bodyj.isGround
                phiLambda_pp = lambda*phiLambda_ppii; % [4x4]
            else % [8x8]
                phiLambda_pp = lambda*[phiLambda_ppii  phiLambda_ppij;
                                       phiLambda_ppji phiLambda_ppjj]; 
            end 
        end
    end %methods
end %class