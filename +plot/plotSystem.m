function plotSystem(sys,frames)
% plotSystem  plots the bodies in the simEngine3D system.
%
%	DRAWFRAME(sys)
%	DRAWFRAME(sys, frames)
%
% sys        : created with the system3D class.
% frames = 1 : Show body reference frames.
%        = 0 : Don't show body reference frames.
%

if nargin == 1
    frames = 0;
end



figure()

%% PLOT GLOBAL REFERNCE FRAME

plot.drawframe(sys.r,sys.p,1,1) % plot GLOBAL RF
xlabel('X');
ylabel('Y');
zlabel('Z');

%% PLOT BODIES
if sys.nBodies == 0; disp('No bodies in the system.'); return; end;


hold on;
for i = 1:sys.nBodies % plot bodies in system
    
    % plot marker for every body
    r = [sys.body{i}.r(1);sys.body{i}.r(2);sys.body{i}.r(3)];
    scatter3(r(1),r(2),r(3),500);
    
    %optionally plot body reference frames
    if frames 
        if sys.body{i}.isGround
            plot.drawframe(sys.body{i}.r,sys.body{i}.p,[],2)
        else
            plot.drawframe(sys.body{i}.r,sys.body{i}.p)
        end
    end
    
    % add body labels
    r_text = r + 0.1; %adjust so label is not right over body origin
    if sys.body{i}.isGround
        text(r_text(1),r_text(2),r_text(3),['Body ' num2str(sys.body{i}.ID) ', ground'],'FontWeight','bold');
    else
        text(r_text(1),r_text(2),r_text(3),['Body ' num2str(sys.body{i}.ID)],'FontWeight','bold');
    end
    
    % IF POINTS ON BODY, PLOT THEM
    if sys.body{i}.nPoints ~= 0
        for j = 1:sys.body{i}.nPoints % plot bodies in system
            sbar = sys.body{i}.point{j}; % local position of point
            Asbar = utility.p2A(sys.body{i}.p)*sbar; % rotated to global rf
            rAsbar = r + Asbar; % global position of point
            
            plot3([r(1); rAsbar(1)],[r(2); rAsbar(2)],[r(3); rAsbar(3)],'ks-') % line from BODY RF to point
            
            % add point labels
            rAsbar_text = rAsbar + 0.05; %adjust so label is not right over point
            text(rAsbar_text(1),rAsbar_text(2),rAsbar_text(3),num2str(j));
        end
    end
end
hold off;
axis equal
end