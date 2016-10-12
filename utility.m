classdef utility
    % Filename: utility.m
    % Author:   Samuel Acuña
    % Date:     11 Oct 2016
    % About:
    % a collection of helpful static functions, these are not
    % associated with any specific class.
    % 
    % invoke functions as, for example: utility.p2A(p)
    
    methods (Static)
        function A = p2A(p) % orientation matrix A from euler parameters p
            % from ME751_f2016 slide 14 from lecture 9/21/16
            % input:
            %   p = [e0;e1;e2;e3] euler parameters
            % output:
            %   A = [3x3] orientation matrix
            e0 = p(1);
            e1 = p(2);
            e2 = p(3);
            e3 = p(4);
            A = 2*[e0^2+e1^2-0.5, e1*e2-e0*e3, e1*e3+e0*e2;
                e1*e2+e0*e3, e0^2+e2^2-0.5, e2*e3-e0*e1;
                e1*e3-e0*e2, e2*e3+e0*e1, e0^2+e3^2-0.5];
        end
        function p = A2p(A) % euler parameters p from orientation matrix A 
            % from ME751_f2016 slides 20-25 from lecture 9/21/16
            % input:
            %   A = [3x3] orientation matrix
            % output:
            %   p = [e0;e1;e2;e3] euler parameters

            e0 = sqrt((trace(A)+1)/4); % the sign of e0 is arbitrary
            
            if e0 ~= 0
                e1 = (A(3,2)-A(2,3))/(4*e0);
                e2 = (A(1,3)-A(3,1))/(4*e0);
                e3 = (A(2,1)-A(1,2))/(4*e0);

            elseif e0 == 0
                disp('Not implemented yet, see slide 25 (9/21/16) or page 341')
                return
            end
            
            p =[e0;e1;e2;e3];
        end
        function R1 = R1(theta) % rotation matrix: about X axis by theta degrees
            theta = theta*pi/180; %convert to radians
            R1 = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
        end
        function R2 = R2(theta) % rotation matrix: about Y axis by theta degrees
            theta = theta*pi/180; %convert to radians
            R2 = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
        end
        function R3 = R3(theta) % rotation matrix: about Z axis by theta degrees
            theta = theta*pi/180; %convert to radians
            R3 = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
        end
        function B = Bmatrix(p,abar) % B matrix, r-p formulation
            % from ME751_f2016 slide 12 from lecture 9/28/16
            % input:
            %   p = [e0;e1;e2;e3] euler parameters
            %   abar = [3x1], position of point expressed in BODY RF
            % output:
            %   B = [3x4] matrix
            e0 = p(1);
            e = p(2:4);
        end
    end
end
