classdef dp2 < handle
    % Filename: dp2.m
    % Author:   Samuel Acuña
    % Date:     17 Oct 2016
    %
    % About:
    % This class handles instances of DP2 constraints. It Uses the r-p
    % formulation (euler parameters). Removes 1 degree of Freedom.
    %
    % DP2 constraint: 
    % the dot product between vector aBari and a vector PiQj 
    % from body i to body j assume a specified value f. f is often 0,
    % meaning aBari and PiQj are orthogonal vectors.
    
    properties
        system;     % parent system3D object to which this constraint is a member. 
        rDOF = 1;   % removes 1 degree of freedom
        bodyi;      % body i
        bodyj;      % body j
        aBari_tail; % ID number for point aBari_tail on body i, tail of aBari vector
        aBari_head; % ID number for point aBari_head on body i, head of aBari vector
        Pi;         % ID number for point P on body i, tail of PiQj vector
        Qj;         % ID number for point Q on body j, head of PiQj vector
        f;          % prescribed constraint, often 0, but can be a function of t
        fdot;       % derivative of f
        fddot;      % derivative of fdot
    end
    properties (Dependent)
        aBari;      % vector in body i RF
        PiQj;       % vector in GLOBAL RF
        phi;        % value of the expression of the constraint PHI^dp2
        nu;         % right-hand side of the velocity equation
        gammaHat;   % right-hand side of the acceleration equation, in r-p formulation
        phi_r;      % partial derivative of constraint with respect to r
        phi_p;      % partial derivative of constraint with respect to p
    end
    
    methods
        %constructor function
        function cons = dp2(system, bodyi,aBari_tailID,aBari_headID,PiID,bodyj,QjID,f,fdot,fddot,t) %constructor function
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
            cons.aBari_tail = aBari_tailID;
            cons.aBari_head = aBari_headID;
            cons.Pi = PiID;
            cons.bodyj = bodyj;
            cons.Qj = QjID;
            cons.f = f;
            cons.fdot = fdot;
            cons.fddot = fddot;
            
            if abs(cons.phi) > 1e-4
                warning('Initial conditions for ''dp2'' are not consistent. But solution will converge so constraints are satisfied.')
            end
        end
        function aBari = get.aBari(cons) % vector in body i RF
            % aBari = aBari_head - aBari_tail
            aBari = cons.bodyi.point{cons.aBari_head} - cons.bodyi.point{cons.aBari_tail};
        end
        function PiQj = get.PiQj(cons) % vector in GLOBAL RF
            % PiQj = vector from Pi to Qj
            PiQj = utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj);
        end
        function phi = get.phi(cons) % value of the expression of the constraint PHI^dp2
            % from ME751_f2016 slide 14 from lecture 09/26/16
            % phi : [1x1]
            
            % see if constraint is a function
            if isa(cons.f, 'function_handle')
                fVal = cons.f(cons.system.time); %evaluate at system time
            else % constant value
                fVal = cons.f;
            end
            
            phi = (cons.bodyi.A*cons.aBari)'*utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj) - fVal;
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
            
            gammaHat = -(cons.bodyi.A*cons.aBari)'*utility.Bmatrix(cons.bodyj.pdot,cons.bodyj.point{cons.Qj})*cons.bodyj.pdot + ...
                        (cons.bodyi.A*cons.aBari)'*utility.Bmatrix(cons.bodyi.pdot,cons.bodyi.point{cons.Pi})*cons.bodyi.pdot + ...
                       -(utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj))'*utility.Bmatrix(cons.bodyi.pdot,cons.aBari)*cons.bodyi.pdot + ...
                       -2*(utility.Bmatrix(cons.bodyi.p,cons.aBari)*cons.bodyi.pdot)'*utility.dijdot(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj) + fddotVal;
        end
        function phi_r = get.phi_r(cons) 
            % from ME751_f2016 slide 15 from lecture 9/28/16
            % One body can be the ground. In this case, the number of columns
            % in the Jacobian is half. there are no partial derivatives with 
            % respect to r or p. Thus, we must properly dimension the size of
            % the vectors/matrices because the grounded body does not have 
            % any generalized coordinates.
            % phi_r : [1x6] normally, unless grounded, then [1x3]
            
            phi_ri = -(cons.bodyi.A*cons.aBari)';
            phi_rj =  (cons.bodyi.A*cons.aBari)';
            
            if cons.bodyi.isGround
                phi_r = [phi_rj];
            elseif cons.bodyj.isGround
                phi_r = [phi_ri];
            else 
                phi_r = [phi_ri, phi_rj]; 
            end
        end
        function phi_p = get.phi_p(cons)
            % from ME751_f2016 slide 13 from lecture 9/28/16
            % One body can be the ground. In this case, the number of columns
            % in the Jacobian is half. there are no partial derivatives with 
            % respect to r or p. Thus, we must properly dimension the size of
            % the vectors/matrices because the grounded body does not have 
            % any generalized coordinates.
            % phi_p : [1x8] normally, unless grounded, then [1x4]
            
            phi_pi = (utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj))'*utility.Bmatrix(cons.bodyi.p,cons.aBari) + ...
                    -(cons.bodyi.A*cons.aBari)'*utility.Bmatrix(cons.bodyi.p,cons.bodyi.point{cons.Pi});
            phi_pj = (cons.bodyi.A*cons.aBari)'*utility.Bmatrix(cons.bodyj.p,cons.bodyj.point{cons.Qj});
            
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
            % from ME751_f2016 slide 35,37 from lecture 10/17/16
            % phiLambda_rr : [6x6] normally, unless grounded, then [3x3]
            % inputs:
            %    lambda : scalar value
            
            dim_r = length(cons.phi_r);
            phiLambda_rr = lambda*zeros(dim_r,dim_r); % [3x3] or [6x6]
        end
        function  phiLambda_rp = phiLambda_rp(cons,lambda)
            % partial derivative of (phi_r*lambda) with respect to p (orientation)
            %from ME751_f2016 slide 35,37 from lecture 10/17/16
            % phiLambda_rp : [6x8] normally, unless grounded, then [3x4]
            % inputs:
            %    lambda : scalar value
            
            phiLambda_rpii = -utility.Bmatrix(cons.bodyi.p,cons.aBari);
            phiLambda_rpjj = zeros(3,4);
            phiLambda_rpij = zeros(3,4);
            phiLambda_rpji = utility.Bmatrix(cons.bodyi.p,cons.aBari);
            
            if cons.bodyi.isGround
                phiLambda_rp = lambda*phiLambda_rpjj; % [3x4]
            elseif cons.bodyj.isGround
                phiLambda_rp = lambda*phiLambda_rpii; % [3x4]
            else % [6x8]
                phiLambda_rp = lambda*[phiLambda_rpii  phiLambda_rpij;
                                       phiLambda_rpji  phiLambda_rpjj]; 
            end 
        end
        function  phiLambda_pr = phiLambda_pr(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to r (position)
            %from ME751_f2016 slide 35,37 from lecture 10/17/16
            % phiLambda_pr : [8x6] normally, unless grounded, then [4x3]
            % inputs:
            %    lambda : scalar value
            
            phiLambda_prii = -utility.Bmatrix(cons.bodyi.p,cons.aBari)';
            phiLambda_prjj = zeros(4,3);
            phiLambda_prij = utility.Bmatrix(cons.bodyi.p,cons.aBari)';
            phiLambda_prji = zeros(4,3);
            
            if cons.bodyi.isGround
                phiLambda_pr = lambda*phiLambda_prjj; % [3x4]
            elseif cons.bodyj.isGround
                phiLambda_pr = lambda*phiLambda_prii; % [3x4]
            else % [6x8]
                phiLambda_pr = lambda*[phiLambda_prii  phiLambda_prij;
                                       phiLambda_prji  phiLambda_prjj];
            end
        end
        function  phiLambda_pp = phiLambda_pp(cons,lambda)
            % partial derivative of (phi_p*lambda) with respect to p (orientation)
            %from ME751_f2016 slide 35,37 from lecture 10/17/16
            % will be [4x4] or [8x8], depends on if there is a grounded body
            % inputs:
            %    lambda : scalar value
            Dij = utility.dij(cons.bodyi,cons.Pi,cons.bodyj,cons.Qj);
            phiLambda_ppii = utility.Kmatrix(cons.aBari,Dij) + ...
                            -utility.Kmatrix(cons.bodyi.point{cons.Pi},cons.bodyi.A*cons.aBari) + ...
                            -utility.Bmatrix(cons.bodyi.p,cons.aBari)'*utility.Bmatrix(cons.bodyi.p,cons.bodyi.point{cons.Pi}) + ...
                            -utility.Bmatrix(cons.bodyi.p,cons.bodyi.point{cons.Pi})'*utility.Bmatrix(cons.bodyi.p,cons.aBari);
            phiLambda_ppjj = utility.Kmatrix(cons.bodyj.point{cons.Qj},cons.bodyi.A*cons.aBari);
            phiLambda_ppij = utility.Bmatrix(cons.bodyi.p,cons.aBari)'*utility.Bmatrix(cons.bodyj.p,cons.bodyj.point{cons.Qj});
            phiLambda_ppji = utility.Bmatrix(cons.bodyj.p,cons.bodyj.point{cons.Qj})'*utility.Bmatrix(cons.bodyi.p,cons.aBari);
            
            if cons.bodyi.isGround
                phiLambda_pp = lambda*phiLambda_ppjj; % [4x4]
            elseif cons.bodyj.isGround
                phiLambda_pp = lambda*phiLambda_ppii; % [4x4]
            else % [8x8]
                phiLambda_pp = lambda*[phiLambda_ppii  phiLambda_ppij;
                                       phiLambda_ppji  phiLambda_ppjj]; 
            end 
        end
    end
    
end