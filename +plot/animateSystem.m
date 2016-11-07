function animateSystem(sys,state)
% simple animation of system3D.
% INPUTS:
%   sys   : system3D object
%   state : cell array of system states throughout time, this is generated
%           by kinematicsAnalysis, inverseDynamicsAnalysis, and dynamicsAnalysis functions

% determine time step
timeStep = state{2}.time-state{1}.time;

% determine number of free bodies in simulation
nFreeBodies = length(state{1}.r)/3;

% determine size of axes
xMin = -1; xMax = 1; yMin = -1; yMax = 1; zMin = -1; zMax = 1;
for i = 1:length(state)
    for j = 1:nFreeBodies % plot free bodies
            x = state{i}.r(3*j-2); y = state{i}.r(3*j-1); z = state{i}.r(3*j);
    end
    if x < xMin; xMin = x; end
    if x > xMax; xMax = x; end
    if y < yMin; yMin = y; end
    if y > yMax; yMax = y; end
    if z < zMin; zMin = z; end
    if z > zMax; zMax = z; end
end
border = 0.5;
axisWindow = [xMin-border xMax+border yMin-border yMax+border zMin-border zMax+border];

% allocate for body handles
bodyHandles = cell(nFreeBodies,1);


% iterate over the time grid 
% (each state{i} is a snapshot of the system in time) 
for i = 1:length(state)
    phIndex = 1;
    if state{i}.time == state{1}.time % if first time step, set up figure.
        figure();
        fig = gcf;
        fig.Color = [1 1 1]; % set background color to white
        hold on
        for j = 1:sys.nBodies % plot grounded bodies as frames
            if sys.body{j}.isGround
                 plot.drawframe(sys.body{i}.r,sys.body{i}.p,[],2)
            end
        end
        for j = 1:nFreeBodies % plot free bodies
            bodyID = sys.bodyIDs(j); % pull current free-body ID
            
            % plot Center of Mass
            r = [state{i}.r(3*j-2);state{i}.r(3*j-1);state{i}.r(3*j)];
            p = [state{i}.p(4*j-3);state{i}.p(4*j-2);state{i}.p(4*j-1);state{i}.p(4*j)];
            A = utility.p2A(p);
            bodyHandles{j} = scatter3(r(1),r(2),r(3),500,sys.body{bodyID}.color,'filled');
            %bodyHandles{j} = plot3([0;state{i}.r(3*j-2)],[0;state{i}.r(3*j-1)],[0;state{i}.r(3*j)],'b-','LineWidth',2);
            
            
            % Plot points on the body
            if sys.body{bodyID}.nPoints ~= 0
                for k = 1:sys.body{bodyID}.nPoints % plot points of bodies in system
                    sbar = sys.body{bodyID}.point{k}; % local position of point
                    Asbar = A*sbar; % rotated to global rf
                    rAsbar = r + Asbar; % global position of point
                    pointHandles{phIndex} = plot3([r(1); rAsbar(1)],[r(2); rAsbar(2)],[r(3); rAsbar(3)],'ks-','MarkerSize',12); % line from BODY RF to point
                    phIndex = phIndex+1;
                end
            end
        end
        hold off
        view(98,12); % default viewing angle
        axis equal
        axis(axisWindow); % set axes size
        pause;
    else % update points to this moment of time
        for j = 1:nFreeBodies % update free bodies
            r = [state{i}.r(3*j-2);state{i}.r(3*j-1);state{i}.r(3*j)];
            p = [state{i}.p(4*j-3);state{i}.p(4*j-2);state{i}.p(4*j-1);state{i}.p(4*j)];
            A = utility.p2A(p);
            set(bodyHandles{j},'XData',r(1),'YData',r(2),'ZData',r(3));
            %set(bodyHandles{j},'XData',[0;state{i}.r(3*j-2)],'YData',[0;state{i}.r(3*j-1)],'ZData',[0;state{i}.r(3*j)]);
            
            bodyID = sys.bodyIDs(j); % pull current free-body ID
            % update points on the body
            if sys.body{bodyID}.nPoints ~= 0
                for k = 1:sys.body{bodyID}.nPoints % plot points of bodies in system
                    sbar = sys.body{bodyID}.point{k}; % local position of point
                    Asbar = A*sbar; % rotated to global rf
                    rAsbar = r + Asbar; % global position of point
                    set(pointHandles{phIndex},'XData',[r(1); rAsbar(1)],'YData',[r(2); rAsbar(2)],'ZData',[r(3); rAsbar(3)]);
                    phIndex = phIndex+1;
                end
            end
        end
        
        drawnow; % redraw those updated points on the plot
    end
    pause(timeStep); % wait an actual time step
end

end