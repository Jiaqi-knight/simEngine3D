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
border = 0.05;
axisWindow = [xMin-border xMax+border yMin-border yMax+border zMin-border zMax+border];

% allocate for body handles
bodyHandles = cell(nFreeBodies,1);

% iterate over the time grid 
% (each state{i} is a snapshot of the system in time) 
for i = 1:length(state)
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
            bodyHandles{j} = scatter3(state{i}.r(3*j-2),state{i}.r(3*j-1),state{i}.r(3*j),500);
            %bodyHandles{j} = plot3([0;state{i}.r(3*j-2)],[0;state{i}.r(3*j-1)],[0;state{i}.r(3*j)],'b-','LineWidth',2);
        end
        hold off
        view(98,12);
        axis equal
        axis(axisWindow); % set axes size
    else % update points to this moment of time
        for j = 1:nFreeBodies % update free bodies
            set(bodyHandles{j},'XData',state{i}.r(3*j-2),'YData',state{i}.r(3*j-1),'ZData',state{i}.r(3*j));
            %set(bodyHandles{j},'XData',[0;state{i}.r(3*j-2)],'YData',[0;state{i}.r(3*j-1)],'ZData',[0;state{i}.r(3*j)]);
        end
        drawnow; % redraw those updated points on the plot
    end
    
    pause(timeStep); % wait an actual time step
end


end