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
        
    end
end
