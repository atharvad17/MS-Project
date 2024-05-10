function [vel_command, collisions] = compute_vel_reynolds(self, p_swarm, r_agent, dt, time, start_t, dur, att_id, vic_id, dev_y, pos_csv, dist_csv, col_csv, info_csv)
    
    % Modif:    a cohesion term has been added to make the agents get
    %           closer when they are farther than r0_rep.
    %
    % Inputs:
    %   self: handle to the drone swarm
    %   p_swarm: swarm parameters
    %   r_agent: safety radius of agents
    %   dt: time step
    %   time: current simulation time
    %   start_t: start time for GPS spoofing
    %   dur: duration of GPS spoofing
    %   att_id: ID of the drone to be attacked
    %   vic_id: ID of the victim drone
    %   dev_y: deviation in y-direction for GPS spoofing
    %   pos_csv: file path for position data CSV
    %   dist_csv: file path for distance data CSV
    %   col_csv: file path for collision data CSV
    %   info_csv: file path for additional info CSV
    %
    % Outputs:
    %   vel_command: commanded velocities for every agent
    %   collisions: [nb_agent_collisions nb_obs_collisions min_dist_obs]
    %
    
    
    %% Initialize variables
    
    pos = self.get_pos_ned();
    vel = self.get_vel_ned(); 

    % Initialize variables
    nb_agents = self.nb_agents;  
    M = zeros(nb_agents, nb_agents);    % Neighborhood matrix 
    dist_mat = zeros(nb_agents, nb_agents);    % Distance matrix 
    vel_rep = zeros(3, nb_agents);      % Repulsion velocity 
    vel_fric = zeros(3, nb_agents);     % Velocity matching velocity  
    vel_wall = zeros(3, nb_agents);     % Arena repulsion velocity  
    vel_obs = zeros(3, nb_agents);      % Obstacle repulsion velocity  
    vel_command = zeros(3, nb_agents);  % Total commanded velocity 
    
    % Initialization 
    nb_agent_collisions = 0; % Nb of collisions among agents
    nb_obs_collisions = 0; % Nb of collisions against obstacles
    min_dist_obs = 20; % Minimal distance between any agent and the obstacle. Set it to a relatively high value in the beginning

    
    %% Compute velocity commands

    % GPS spoofing
    if (time >= start_t) && (time < (dur + start_t))
        pos(1:2, att_id) = pos(1:2, att_id) + [0; dev_y];
    end
    
    % Compute cohesion
    cohesion_command = compute_cohesion(self, p_swarm);
    
    % Compute separation
    separation_command = compute_separation(self, p_swarm, r_agent);
    
    % Compute alignment
    alignment_command = compute_alignment(self, p_swarm, r_agent);
    
    for agent = 1:nb_agents
            
        % Compute agent-agent distance matrix
        p_rel = pos - pos(:, agent); 
        dist = sqrt(sum((p_rel.^2), 1)); 
        dist_mat(agent, :) = dist; 

        % Define neighbours list
        neig_list = (1:nb_agents)'; 
        neig_list = neig_list(dist ~= 0);
        
        % Count collisions      
        nb_agent_collisions = nb_agent_collisions + sum(dist < 2 * r_agent) - 1; 

        % Initialize number of neighbours
        nb_neig = nb_agents - 1;

        % Constraint on neighborhood given by the euclidean distance
        if isfield(p_swarm, 'r') 
            neig_list = neig_list(dist(neig_list) < p_swarm.r); 
            nb_neig = length(neig_list); 
        end

        % Constraint on neighborhood given by the topological distance
        if isfield(p_swarm, 'max_neig')
            if nb_neig > p_swarm.max_neig
                [~, idx] = sort(dist(neig_list));
                neig_list = neig_list(idx(1:p_swarm.max_neig));
                nb_neig = p_swarm.max_neig;
            end
        end

        % Adjacency matrix (asymmetric in case of limited fov)
        M(agent, neig_list) = 1; 

        
        %% Compute different contributions

        if nb_neig ~= 0
            v_rel = vel - vel(:, agent); 
            v_rel_norm = sqrt(sum((v_rel.^2), 1)); 

            % Compute vel and pos unit vector between two agents
            p_rel_u = -p_rel ./ dist; 
            v_rel_u = -v_rel ./ v_rel_norm;

            for agent2 = neig_list' 
                
                % Repulsion and attraction
                if dist(agent2) < p_swarm.r0_rep   % repulsion
                    vel_rep(:, agent) = vel_rep(:, agent) + ...
                        p_swarm.p_rep * (p_swarm.r0_rep - dist(agent2)) * p_rel_u(:, agent2);
                else  % attraction
                    vel_rep(:, agent) = vel_rep(:, agent) + ...
                        p_swarm.p_rep * (dist(agent2) - p_swarm.r0_rep) * -p_rel_u(:, agent2);
                end

                % Velocity alignement
                v_fric_max = get_v_max(p_swarm.v_fric, dist(agent2) - p_swarm.r0_fric, p_swarm.a_fric, p_swarm.p_fric);
                
                if v_rel_norm(agent2) > v_fric_max
                    vel_fric(:, agent) = vel_fric(:, agent) + ...
                        p_swarm.C_fric * (v_rel_norm(agent2) - v_fric_max) * v_rel_u(:, agent2);
                end

            end
        end
        
        
        %% Wall and obstacle avoidance
        
        % Add arena repulsion effect
        if p_swarm.is_active_arena
            unit = eye(3); 
            %On each axis we have the two repulsions
            for axis = 1:3
                %On each axis there is two forces (each side of the arena)
                for dir = 1:2
                    dist_ab = abs(pos(axis, agent) - p_swarm.x_arena(axis, dir));

                    %Compute velocity of wall shill agent toward center of the arena
                    v_wall_virtual = unit(:, axis) .* p_swarm.v_shill;
                    
                    if dir == 2
                        v_wall_virtual = -v_wall_virtual; % replusion from the opposite direction
                    end

                    %Compute relative velocity (Wall - Agent)
                    vel_ab = sqrt(sum((vel(:, agent) - v_wall_virtual).^2)); 

                    v_wall_max = get_v_max(0, dist_ab - p_swarm.r0_shill, p_swarm.a_shill, p_swarm.p_shill);

                    if vel_ab > v_wall_max
                        vel_wall(:, agent) = vel_wall(:, agent) + ...
                            (vel_ab - v_wall_max) * (v_wall_virtual - vel(:, agent)) ./ vel_ab;
                    end
                end
            end
        end

        % Compute spheric effect
        if p_swarm.is_active_spheres

            for obs = 1:p_swarm.n_spheres
                % Get obstacle center and radius
                c_obs = p_swarm.spheres(1:3, obs); 
                r_obs = p_swarm.spheres(4, obs);

                % Compute distance agent(a)-obstacle(b)
                dist_ab = sqrt(sum((pos(:, agent) - c_obs).^2)) - r_obs;
                nb_obs_collisions = nb_obs_collisions + sum(dist_ab < r_agent);

                % Set the virtual speed of the obstacle direction out of
                % the obstacle
                v_obs_virtual = (pos(:, agent) - c_obs) / (dist_ab + r_obs) * p_swarm.v_shill;

                % Compute relative velocity agent-obstacle
                vel_ab = sqrt(sum((vel(:, agent) - v_obs_virtual).^2));
                
                if dist_ab < min_dist_obs
                    min_dist_obs = dist_ab;
                end
      
                v_obs_max = get_v_max(0, dist_ab - p_swarm.r0_shill, p_swarm.a_shill, p_swarm.p_shill);

                if vel_ab > v_obs_max
                    vel_obs(:, agent) = vel_obs(:, agent) + (vel_ab - v_obs_max) * (v_obs_virtual - vel(:, agent)) ./ vel_ab;
                end
            end
        end

        % Compute cylindric effect 
        if p_swarm.is_active_cyl

            for obs = 1:p_swarm.n_cyl 
                % Get obstacle center and radius
                c_obs = p_swarm.cylinders(1:2, obs); 
                r_obs = p_swarm.cylinders(3, obs); 

                % Compute distance agent(a)-obstacle(b)
                dist_ab = sqrt(sum((pos(1:2, agent) - c_obs).^2)) - r_obs; 
                
                % Record distance between agent and obstacle
                dist_all = [time, agent, dist_ab];
                writematrix(dist_all, dist_csv, 'Delimiter', ',', 'WriteMode', 'append');
                
                % Record collisions between victim drone and obstacle
                if vic_id == 0
                    nb_obs_collisions = nb_obs_collisions + sum(dist_ab < r_agent); 
                else
                    if agent == vic_id
                        nb_obs_collisions = nb_obs_collisions + sum(dist_ab < r_agent); 
                    end
                end
                
                % Record distance between victim drone and obstacle
                if agent == vic_id
                    res_dist_mat = [time, dist_ab];
                    writematrix(res_dist_mat, info_csv, 'Delimiter', ',', 'WriteMode', 'append');
                end
                
                % Set the virtual speed of the obstacle direction out of
                % the obstacle
                v_obs_virtual = (pos(1:2, agent) - c_obs) / (dist_ab + r_obs) * p_swarm.v_shill;

                % Compute relative velocity agent-obstacle
                vel_ab = sqrt(sum((vel(1:2, agent) - v_obs_virtual).^2));
                
                if dist_ab < min_dist_obs
                    min_dist_obs = dist_ab;
                end
                
                v_obs_max = get_v_max(0, dist_ab - p_swarm.r0_shill, p_swarm.a_shill, p_swarm.p_shill);

                if vel_ab > v_obs_max
                    vel_obs(1:2, agent) = vel_obs(1:2, agent) + (vel_ab - v_obs_max) * (v_obs_virtual - vel(1:2, agent)) ./ vel_ab;
                end

            end

        end
        
        %% Sum agent-agent and obstacle contributions
        vel_command(:, agent) = vel_rep(:, agent) + vel_fric(:, agent) + vel_obs(:, agent) + vel_wall(:, agent) + cohesion_command(:, agent) + separation_command(:, agent) + alignment_command(:, agent);

        % Add self propulsion OR migration term
        v_norm = sqrt(sum((vel(:, agent).^2), 1));

        if p_swarm.is_active_migration % migration 
            vel_command(:, agent) = vel_command(:, agent) + p_swarm.v_ref * p_swarm.u_ref; 
        elseif p_swarm.is_active_goal
            x_goal_rel = p_swarm.x_goal(:, agent) - pos(:, agent);
            u_goal = x_goal_rel / norm(x_goal_rel);
            vel_command(:, agent) = vel_command(:, agent) + p_swarm.v_ref * u_goal;
        else
            % self-propulsion
            if v_norm > 0
                vel_command(:, agent) = vel_command(:, agent) + p_swarm.v_ref * vel(:, agent) / v_norm;
            end
        end
    end
    
    % Record velocity
    vel_record = reshape(vel_command(1:2, :), 1, []);
    writematrix(vel_record, "./velAttack.csv", 'Delimiter', ',', 'WriteMode', 'append');

    %% Compute collisions and bound velocities and accelerations

    % Total number of collisions per time step
    nb_agent_collisions = nb_agent_collisions / 2; % reciprocal 
    collisions = [nb_agent_collisions nb_obs_collisions min_dist_obs];

    % Add random effect on velocities 
    if isfield(p_swarm, 'c_r')
        vel_command = vel_command + p_swarm.c_r * randn(3, nb_agents);
    end

    % Bound velocities and acceleration
    if ~isempty(p_swarm.max_v) 
        vel_cmd_norm = sqrt(sum((vel_command.^2), 1)); 
        v_norm = sqrt(sum((vel.^2), 1));
        
        idx_to_bound = (vel_cmd_norm > p_swarm.max_v); 
        if sum(idx_to_bound) > 0 
            vel_command(:, idx_to_bound) = p_swarm.max_v * ...
                vel_command(:, idx_to_bound) ./ repmat(vel_cmd_norm(idx_to_bound), 3, 1);

        end
    end
    if ~isempty(p_swarm.max_a) 
        accel_cmd = (vel_command - vel) / dt;
        accel_cmd_norm = sqrt(sum(accel_cmd.^2, 1));
        idx_to_bound = (accel_cmd_norm > p_swarm.max_a | accel_cmd_norm < - p_swarm.max_a);
        if sum(idx_to_bound) > 0
            vel_command(:, idx_to_bound) = vel(:, idx_to_bound) + ...
                dt * p_swarm.max_a * accel_cmd(:, idx_to_bound) ./ ...
                repmat(accel_cmd_norm(idx_to_bound), 3, 1);
        end
    end
     
    % Record position
    if (start_t == 0 && dur == 0)
        pos_mat = zeros(1, 1 + 2 * p_swarm.nb_agents);
        pos_mat(1) = time;
        for i = 1:p_swarm.nb_agents
            pos_mat(2 * i) = pos(1, i);
            pos_mat(2 * i + 1) = pos(2, i);
        end
        writematrix(pos_mat, pos_csv, 'Delimiter', ',', 'WriteMode', 'append');
    end
   
    % Record collisions
    res_col_mat = [time, nb_obs_collisions, start_t, dur];
    writematrix(res_col_mat, col_csv, 'Delimiter', ',', 'WriteMode', 'append');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate V fric max
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [v_fricmax] = get_v_max(v_fric, r, a, p)

    if r < 0
        v_fricmax = 0;
    elseif r * p > 0 && r * p < a / p
        v_fricmax = r * p;
    else
        v_fricmax = sqrt(2 * a * r - a^2 / p^2);
    end
    if v_fricmax < v_fric
        v_fricmax = v_fric;
    end
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

function separation_command = compute_separation(swarm, p_swarm, r_agent)
    % Compute separation between drones
    nb_agents = swarm.nb_agents;
    pos = swarm.get_pos_ned();
    separation_command = zeros(3, nb_agents);
    
    % Adjust parameters as needed
    separation_radius = 2 * r_agent; % Adjust separation radius as needed
    
    for i = 1:nb_agents
        % Compute separation for each agent
        for j = 1:nb_agents
            if i ~= j
                % Compute relative position
                rel_pos = pos(:, j) - pos(:, i);
                
                % Compute distance between agents
                dist_ij = norm(rel_pos);
                
                % Check if the agent is within the separation radius
                if dist_ij < separation_radius
                    % Compute separation force
                    k_separation = 0.5; % Adjust separation force coefficient as needed
                    separation_force = k_separation * (1 / dist_ij - 1 / separation_radius) * rel_pos;
                    
                    % Accumulate separation force
                    separation_command(:, i) = separation_command(:, i) + separation_force;
                end
            end
        end
    end
end

function alignment_command = compute_alignment(swarm, p_swarm, r_agent)
    % Compute alignment between drones
    nb_agents = swarm.nb_agents;
    vel = swarm.get_vel_ned();
    alignment_command = zeros(3, nb_agents);
    
    % Adjust parameters as needed
    alignment_radius = 5 * r_agent; % Adjust alignment radius as needed
    
    for i = 1:nb_agents
        % Compute alignment for each agent
        for j = 1:nb_agents
            if i ~= j
                % Compute relative velocity
                rel_vel = vel(:, j) - vel(:, i);
                
                % Compute distance between agents
                dist_ij = norm(rel_vel);
                
                % Check if the agent is within the alignment radius
                if dist_ij < alignment_radius
                    % Compute alignment force
                    k_alignment = 0.1; % Adjust alignment force coefficient as needed
                    alignment_force = k_alignment * rel_vel;
                    
                    % Accumulate alignment force
                    alignment_command(:, i) = alignment_command(:, i) + alignment_force;
                end
            end
        end
    end
end
