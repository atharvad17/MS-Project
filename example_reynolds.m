function example_reynolds(start_t, dur, att_id, vic_id, dev_y, sed, pos_csv, dist_csv, col_csv, info_csv)
%% Clear console and workspace and add project root to path

close all;
%clearvars -except app;

project_root = strcat(extractBefore(mfilename('fullpath'),mfilename),'../..');
addpath(genpath(project_root));


%% Simulation options

DRONE_TYPE = "point_mass"; % swarming mode supports only quadcopter and point_mass
ACTIVE_ENVIRONMENT = true;
DEBUG = true;
VIDEO = true;
CENTER_VIEW_ON_SWARM = false;
SWARM_ALGORITHM = "reynolds"; % either vasarhelyi or olfati_saber

if DEBUG || VIDEO
    results_dirname = strcat('results/results_swarm');
    date_string = datestr(now,'yyyy_mm_dd_HH_MM_SS');
    subfolder = strcat(erase(mfilename,"example_"), '_', date_string);
    results_dirname = strcat(results_dirname, '/', subfolder);
    if ~exist(results_dirname, 'dir')
        mkdir(results_dirname)
    end
end

fontsize = 12;

%% Get changes from GUI

if exist('app', 'var')
    % Simulation parameters
    p_sim.end_time = app.sim_time;
    
    % Drone parameters
    DRONE_TYPE = app.drone_type;
    
    % Swarming parameters
    SWARMING_ALGO = "vasarhelyi";
    p_swarm.nb_agents = app.nb_agents;
    p_swarm.d_ref = app.d_ref;
    p_swarm.v_ref = app.v_ref;
    
    % Map parameters
    ACTIVE_ENVIRONMENT = app.active_environment;
    
    % Debug plot
    DEBUG = app.debug_plot;
end

if DRONE_TYPE == "point_mass"
   SWARM_VIEWER_TYPE = "agent";
elseif DRONE_TYPE == "quadcopter"
   SWARM_VIEWER_TYPE = "drone";
end


%% Call parameters files
% set the seed
p_swarm.seed = sed;

run('param_sim');
run('param_battery');
run('param_physics');
if DRONE_TYPE == "fixed_wing" || DRONE_TYPE == "quadcopter"
    run('param_drone'); 
elseif DRONE_TYPE == "point_mass"
    run('param_drone'); 
end
run('param_map'); % creates map: struct for map params
run('param_swarm');


%% Init Swarm object, Wind, Viewer and other variables


% Init swarm and set positions
swarm = Swarm();
swarm.algorithm = SWARM_ALGORITHM;

for i = 1 : p_swarm.nb_agents
    swarm.add_drone(DRONE_TYPE, p_drone, p_battery, p_sim, p_physics,...
         map);
end
swarm.set_pos(p_swarm.Pos0);

% Init wind
wind = zeros(6,1); % steady wind (1:3), wind gusts (3:6)

% Init variables for history
x0 = [p_swarm.Pos0; zeros(3,p_swarm.nb_agents)];
x_history(1,:) = x0(:);

VIDEO = 1;
DEBUG = 1;
if VIDEO    
    video_filename = strcat(erase(mfilename, "example_"), '_', date_string);
    video_filepath = strcat(results_dirname, '/', video_filename);
    video = VideoWriterWithRate(video_filepath, p_sim.dt_video);
end

% Init viewer
swarm_viewer = SwarmViewer(p_sim.dt_plot, CENTER_VIEW_ON_SWARM);
swarm_viewer.viewer_type = SWARM_VIEWER_TYPE;
states_handle = [];


%% Main simulation loop

disp('Type CTRL-C to exit');
for time = p_sim.start_time:p_sim.dt:p_sim.end_time
    
    % Check if program terminated from GUI
    if exist('app', 'var')
        switch app.StartsimulationSwitch.Value
            case 'Off'
                close all;
                return;
        end
    end
    
    % Get change from GUI
    if exist('app', 'var')
        % Wind parameters
        wind_active = app.wind;
        wind_gust_active = app.wind_gust;
        wind_level = app.wind_level;
        wind_gust_level = app.wind_gust_level;
        
        % Debug plot
        debug_plot = app.debug_plot;
        
        % Orientation of swarm migration
        orientation = app.orientation;
        p_swarm.u_ref = [-cosd(orientation), -sind(orientation), 0]';
    end
    
% ----------------------
%     attack(time);
% ----------------------
    % Compute velocity commands from swarming algorithm
    [vel_c,collisions] = swarm.update_command(p_swarm, p_swarm.r_coll, p_sim.dt, time, start_t, dur, att_id, vic_id, dev_y, pos_csv, dist_csv, col_csv, info_csv);
    
    % Introduce cohesion between drones
    cohesion_command = compute_cohesion(swarm, p_swarm);
    vel_c = vel_c + cohesion_command;
    
    % Update swarm states and plot the drones
    swarm.update_state(wind, time);
   
    % Plot state variables for debugging
    if DEBUG
        swarm.plot_state(time, p_sim.end_time, ...
        1, p_sim.dt_plot, collisions, p_swarm.r_coll/2);
    end
    
    % Update video
    if VIDEO
        swarm_viewer.update(time, swarm, map);
        video.update(time, swarm_viewer.figure_handle);  
    end
    
end

if VIDEO
    video.close(); 
end

% Close all plots
close all;

disp('Simulation completed successfully');

end




function cohesion_command = compute_cohesion(swarm, p_swarm)
    % Compute cohesion between drones
    nb_agents = swarm.nb_agents;
    pos = swarm.get_pos_ned();
    cohesion_command = zeros(3, nb_agents);
    
    % Adjust parameters as needed
    cohesion_radius = 5; % Adjust cohesion radius as needed
    
    for i = 1:nb_agents
        % Compute cohesion for each agent
        for j = 1:nb_agents
            if i ~= j
                % Compute relative position
                rel_pos = pos(:, j) - pos(:, i);
                
                % Compute distance between agents
                dist_ij = norm(rel_pos);
                
                % Check if the agent is within the cohesion radius
                if dist_ij < cohesion_radius
                    % Compute cohesion force
                    k_cohesion = 0.1; % Adjust cohesion force coefficient as needed
                    cohesion_force = k_cohesion * rel_pos;
                    
                    % Accumulate cohesion force
                    cohesion_command(:, i) = cohesion_command(:, i) + cohesion_force;
                end
            end
        end
    end
end
