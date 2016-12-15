classdef tsda3D < handle
    % Filename: tsda3D.m
    % Author:   Samuel Acuña
    % Date:     14 Dec 2016
    %
    % About:
    % This class handles instances of translational spring-damper-actuators
    % (TSDA).  It Uses the r-p formulation (euler parameters). 
    %
    % This is described in lecture 10/07/16 ME751_f2016
    
    properties
        system; % parent system3D object to which this TSDA is a member. 
        bodyi;  % body i
        bodyj;  % body j
        Pi;     % ID number for point P on body i
        Qj;     % ID number for point Q on body j
        k;      % spring stiffness
        l0;     % zero tension length of the spring
        c;      % damping coefficient
        h;      % actuation force, anonymous function, h(l,ldot,time)
        Fi;     % forces3D object to the body i end of the TSDA
        Fj;     % forces3D object to the body j end of the TSDA
    end
    properties (Dependent)
        force;  % expression of the overall TSDA force
        l;      % distance between points Pi and Qj
        ldot;   % rate of change of l
        dij;    % vector from Pi to Qj
        eij;    % dij/l, direction of force
    end
    methods
        function tsda = tsda3D(system,bodyi,PiID,bodyj,QjID,k,l0,c,h) % constructor function
            tsda.system = system;
            tsda.bodyi = bodyi;
            tsda.Pi = PiID;
            tsda.bodyj = bodyj;
            tsda.Qj = QjID;
            tsda.k = k;
            tsda.l0 = l0;
            tsda.c = c;
            tsda.h = h;
            
            tsda.Fi = forces3D(tsda.bodyi,'tsda',tsda.Pi,tsda.force, tsda.eij,tsda);
            tsda.Fj = forces3D(tsda.bodyj,'tsda',tsda.Qj,tsda.force,-tsda.eij,tsda);
        end
        function updateForces(tsda) % update forces3d objects
            tsda.Fi.updateForces(tsda.force, tsda.eij)
            tsda.Fj.updateForces(tsda.force,-tsda.eij)
        end
        function l = get.l(tsda) % distance Pi to Qj
            l = norm(tsda.dij);
        end
        function ldot = get.ldot(tsda) % rate of change of l
            % from haug, eq 11.4.3
            ldot = (tsda.eij)'*utility.dijdot(tsda.bodyi,tsda.Pi,tsda.bodyj,tsda.Qj);
        end
        function force = get.force(tsda) % main expression of force for TSDA element
            force = tsda.k*(tsda.l-tsda.l0) + tsda.c*tsda.ldot + tsda.h(tsda.l,tsda.ldot,tsda.system.time);
        end
        function dij = get.dij(tsda) % vector from Pi to Qj
            dij = utility.dij(tsda.bodyi,tsda.Pi,tsda.bodyj,tsda.Qj);
        end
        function eij = get.eij(tsda) % dij/l
            eij = tsda.dij/tsda.l;
        end
    end %methods
end %class

