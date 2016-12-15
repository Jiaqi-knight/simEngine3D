function animateSystem(sys,state,viewAngle,frameStep)
% simple animation of system3D.
% INPUTS:
%   sys   : system3D object
%   state : cell array of system states throughout time, this is generated
%           by kinematicsAnalysis, inverseDynamicsAnalysis, and dynamicsAnalysis functions
%viewAngle: [Az, El] to set viewing angle of animation figure, optional

if ~exist('viewAngle','var') || isempty(viewAngle)
    viewAngle = [98,12];
end
if ~exist('frameStep','var') || isempty(frameStep)
    frameStep = 1;
end
% determine time step
timeStep = state{2}.time-state{1}.time;

% determine number of free bodies in simulation
nFreeBodies = length(state{1}.r)/3;

% determine number of points on bodies
nFreePoints = 0;
for i = 1:nFreeBodies
    bodyID = sys.bodyIDs(i); % pull current free-body ID
    nFreePoints = nFreePoints + sys.body{bodyID}.nPoints;
end


% allocate for speed
bodyHandles = gobjects(nFreeBodies,1); % empty graphics handles array
pointHandles = gobjects(nFreePoints,1);
bodyColors = zeros(nFreeBodies,3);
rBodies = cell(length(state),1); % body positions
rPoints = cell(length(state),1); % point positions (as lines)

%% PULL DATA
% iterate over the time grid 
% (each state{i} is a snapshot of the system in time) 
for i = 1:length(state)
    pointNum = 1;
    rBody = cell(nFreeBodies,3);
    rPoint = cell(nFreePoints,3);
    for j = 1:nFreeBodies % plot free bodies
        
        % save color of body
        bodyID = sys.bodyIDs(j); % pull current free-body ID
        bodyColors(j,:) = sys.body{bodyID}.color;
        
        % save position of body
        rBody(j,:) = {state{i}.r(3*j-2), state{i}.r(3*j-1), state{i}.r(3*j)};
        r = cell2mat(rBody(j,:));
        
        % save points as lines
        p = [state{i}.p(4*j-3);state{i}.p(4*j-2);state{i}.p(4*j-1);state{i}.p(4*j)];
        A = utility.p2A(p);
        if sys.body{bodyID}.nPoints ~= 0
            for k = 1:sys.body{bodyID}.nPoints
                sbar = sys.body{bodyID}.point{k}; % local position of point
                Asbar = A*sbar; % rotated to global rf
                rAsbar = r + Asbar'; % global position of point
                rPoint(pointNum,:) = {[r(1);rAsbar(1)],[r(2);rAsbar(2)],[r(3);rAsbar(3)]};
                pointNum = pointNum + 1;
            end
        end
    end
    rBodies{i} = rBody;  % save body  positions for this timestep
    rPoints{i} = rPoint; % save point positions for this timestep
end

%% ANIMATE DATA

% set frame
figure();
fig = gcf;
fig.Color = [1 1 1]; % set background color to white

% plot ground bodies 
hold on
for j = 1:sys.nBodies 
    if sys.body{j}.isGround 
        % plot grounded bodies as frames
        plot.drawframe(sys.body{j}.r,sys.body{j}.p,[],2) 
        
        % plot points on ground bodies
        A = utility.p2A(sys.body{j}.p);
        r = sys.body{j}.r;
        if sys.body{j}.nPoints ~= 0
            for k = 1:sys.body{j}.nPoints 
                sbar = sys.body{j}.point{k}; % local position of point
                Asbar = A*sbar; % rotated to global rf
                rAsbar = r + Asbar; % global position of point
                plot3([r(1); rAsbar(1)],[r(2); rAsbar(2)],[r(3); rAsbar(3)],'ks-','MarkerSize',12); % line from BODY RF to point
            end
        end        
    end
end

% plot initial bodies
for j = 1:nFreeBodies 
    bodyHandles(j) = scatter3(rBodies{1}{j,1},rBodies{1}{j,2},rBodies{1}{j,3},500,bodyColors(j,:),'filled');
end

% plot initial points
for j = 1:nFreePoints
    pointHandles(j) = plot3(rPoints{1}{j,1},rPoints{1}{j,2},rPoints{1}{j,3},'ks-','MarkerSize',12); % line from BODY RF to point
end
hold off

% determine size of axes
maxs = [0,0,0]; mins = [0,0,0];
for i = 1:length(state) % find max and min positions over time
    maxs = max([maxs; cell2mat(rBodies{i}); cell2mat(rPoints{i})]);
    mins = min([mins; cell2mat(rBodies{i}); cell2mat(rPoints{i})]);
end
border = 0.5;
axisWindow = [mins(1)-border maxs(1)+border mins(2)-border maxs(2)+border mins(3)-border maxs(3)+border];

% set frame properties
axis equal
axis(axisWindow); % set axes size
view(viewAngle(1),viewAngle(2)); 
pause;

% animate through time
propertyCell = {'XData','YData','ZData'};
for i = 1:frameStep:length(state)
    % update bodies
    set(bodyHandles,propertyCell,rBodies{i});
    
    % update points on bodies
    set(pointHandles,propertyCell,rPoints{i});
    
    drawnow; % redraw updated points on the plot
    %pause(timeStep); % wait an actual time step
end


end