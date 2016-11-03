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
fig = gcf;
fig.Color = [1 1 1]; % set background color to white

%% PLOT GLOBAL REFERENCE FRAME

plot.drawframe(sys.r_global,sys.p_global,0.5,1) % plot GLOBAL RF
xlabel('X');
ylabel('Y');
zlabel('Z');

%% PLOT BODIES
if sys.nBodies == 0; disp('No bodies in the system.'); return; end;

texts = {}; % holds text objects
hold on;
for i = 1:sys.nBodies % plot bodies in system
    
    % plot marker for every body
    r = [sys.body{i}.r(1);sys.body{i}.r(2);sys.body{i}.r(3)];
    if ~sys.body{i}.isGround
        scatter3(r(1),r(2),r(3),500,sys.body{i}.color,'filled');
    end
    
    %optionally plot body reference frames
    if frames 
        if sys.body{i}.isGround
            plot.drawframe(sys.body{i}.r,sys.body{i}.p,[],2)
        else
            plot.drawframe(sys.body{i}.r,sys.body{i}.p)
        end
    end
    
    % add body labels
    r_text = r + 0.3; %adjust so label is not right over body origin
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
            
            plot3([r(1); rAsbar(1)],[r(2); rAsbar(2)],[r(3); rAsbar(3)],'ks-','MarkerSize',12) % line from BODY RF to point
            
            % add point labels
            rAsbar_text = rAsbar + 0.1; %adjust so label is not right over point
            textLabel = ['B' num2str(i) 'P' num2str(j)];
            if length(texts) > 0
                for k = 1:length(texts)
                    if isequal(round(texts{k}.Position,3),round(rAsbar_text',3)) % if text exists at this point, append to label
                        textLabel = [textLabel ',' texts{k}.String];
                        texts{k}.String = '';
                        break;
                    end
                end
            end
            texts{length(texts)+1} = text(rAsbar_text(1),rAsbar_text(2),rAsbar_text(3),textLabel);
        end
    end
end
hold off;
axis equal
end