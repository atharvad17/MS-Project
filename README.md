Clone the two repositories.

git clone https://github.com/lis-epfl/swarmlab.git
git clone https://github.com/DependableSystemsLab/SwarmFuzz.git


We have to copy SwarmFuzz into the Swarmlab folder and also change some settings in Swarmlab to record the swarm's state.

Run the following commands

cp -r SwarmFuzz/fuzz/ swarmlab
cp SwarmFuzz/changed_swarmlab/compute_vel_vasarhelyi.m swarmlab/@Swarm
cp SwarmFuzz/changed_swarmlab/param_swarm.m swarmlab/parameters/
cp SwarmFuzz/changed_swarmlab/example_vasarhelyi.m swarmlab/examples/examples_swarm/
cp SwarmFuzz/changed_swarmlab/Swarm.m swarmlab/@Swarm
cp SwarmFuzz/changed_swarmlab/create_shifted_buildings.m swarmlab/graphics/graphics_map/


cd swarmlab
matlab

Now, replace the files in the codebase with the files from this repository.

Run swarmfuzz(200,300,10,10).
