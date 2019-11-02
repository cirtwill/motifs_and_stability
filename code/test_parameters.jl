#######################################################
### Build a food web, simulate the dynamics,
### remove a species, rinse and repeat
#######################################################


# Run from command line (after setting working directory to that below)
# include(“..\\..\\code\\simulations.jl”) 
# Or, from Linux:
# include("..//..//code//simulations.jl")

## Commands for getting set up
using Pkg
#Pkg.add("BioEnergeticFoodWebs")  # install the package
#Pkg.add("JSON")
# Pkg.add("JLD2")
#Pkg.update()  # check that it's the most up to date one

using BioEnergeticFoodWebs  # tell Julia we want to use the package
using JSON
using JLD
using LightGraphs
#Do we want to allow rewiring after extinctions? If so which type?


# # set working directory
# cd("C:\\Users\\kewo0002\\gitstuff\\motif_role_stability\\data\\Julia\\")

# initialize params
maxreps = 10000; # Giving up after this many
reps = 100 # How many of each combination do we want?
S = 10; # how many species do we want?
Smin=50
Smax=100
Cmin=0.02
Cmax=0.2
		 ## NOTE: C gets converted into links for the function "nichemodel()". If in the conversion the value is < 0, it will throw an error
		 ## Sometimes it can work a bunch of times before we get an error.
start = 0; # where should we start our simulations?
stop = 10; # and how long should they run?
ext_threshold = 0.00001; # what's our threshold for when a species should be considered extinct?

function save_json(net,S,C,i)
	filename = "../data/networks/pre_disturbance/$S/$C/initial_net_$i.json";
	f=open(filename,"w");
	JSON.print(f,net)
	close(f)
end

# check that we don't have any species that are neither preyed upon nor prey apon
# returns false if there are NO isolated species and true if there are TRUlY isolated species
# If true, reject the network
function check_isolated(network)
	colsums = sum(network, dims = 1);
	rowsums = sum(network, dims = 2);
	sumsum  = sum(hcat(colsums',rowsums),dims = 2);
	if 0 in sumsum
		return true;
	else
		return false;
	end
end

# We also don't want any isolated components of 2+ species. 
# Like check_isolated, returns false if there are NO isolated components<N
# If true, reject the network
function check_components(network)
	g = SimpleDiGraph(network)
	# strong=LightGraphs.strongly_connected_components(g)
	# I think we want weakly connected components. Strong means all components are connected along the link direction, so two basals (e.g.) will never be connected.
	# Weakly connected means we ignore the direction.
	weak=LightGraphs.weakly_connected_components(g)
	# Hmm. 
	if length(weak) > 1
		return true
	else
		return false
	end
end

# Want to check whether the random biomasses are allowing all species to persist
# returns true if all species survive
function check_extinction(sim)
	species_below = length(findall(sim[:B].<ext_threshold))
	if species_below > 0
		return false
	else 
		return true
	end
end


# Seperating the building of the web and isolation checking from the sim.
# If sim fails to get coexistence, want to send it back to web building

function web_build(S,C,)
	foodweb = nichemodel(S,C); #has predators in rows, and preys in columns. It can only have 0 and 1.

	# Check we don't have any isolated species
	let iter = 1;
		while check_isolated(foodweb) & (iter < 1000)
			foodweb = nichemodel(S,C)
			iter += 1
			# println(iter)
		end
		if iter == 1000
      		error("ERROR: 100 tries and still getting isolated species")
       	end
	end
	return foodweb
end


# Run through the species range
for S in Smin:10:Smax
	try 
		mkdir("fw_dstrbnc/$S")
	catch
	end
	# And the C range
	for C in Cmin:0.02:Cmax
		println(S," ",C)
		try 
			mkdir("fw_dstrbnc/$S/$C")
		catch
		try 
			mkdir("../data/networks/pre_disturbance/$S/$C")
		catch
		end
		end
		j = 0 # Number of good reps
		let i = 0;
			# Run through up to maxreps iteractions, trying to get 100 good ones
			while (i<maxreps) & (j<reps)
				# Build a niche network with S species and C connectance
				foodweb=web_build(S,C)
				let com_iter = 1;
					while check_components(foodweb) & (com_iter<100)
						foodweb=web_build(S,C)
						check_components(foodweb)
						com_iter += 1
					end
				end
				# Use the default parameters 
				params = model_parameters(foodweb);
				# Set initial biomasses randomly 
				bm = rand(size(foodweb,1));

				# run the simulation
				sim = simulate(params, bm, start = start, stop = stop,use=:nonstiff);
				# save the output - this will be fixed in the next release of BioEnergeticFoodWebs, for now a work-around

				# Check that all of the species can survive
				let iter = 1;
					while !check_extinction(sim) & (iter<100)
						bm = rand(size(foodweb,1))
						sim = simulate(params, bm, start=start, stop=stop, use=:nonstiff);
						iter += 1
					end
					# if iter == 100
						# warning("ERROR: This configuration doesn't result in all species surviving")
					# end
				end


				# Julia is STILL unable to reload data in the same form that it was saved. Alyssa hates this language. Help is not googleable.
				# Attempting to save final biomasses in sim[:p] because that seems to work while sim[:B] is converted to an impossible data format.
				sim[:p][:B]=sim[:B][end,1:end]
				# Reformatting sim[:B] may also work. Going to do both because fuuuuuuuck Julia.
				sim[:B]=collect(sim[:B])
				#BioEnergeticFoodWebs.save(filename = "init_fw_"*string(i), as = json)
				JLDname = "../data/networks/pre_disturbance/$S/$C/initial_$j.jld";
				# Note! If the file already exists, this will fail.
				if check_extinction(sim)
					save(JLDname,"sim",sim)
					net=params[:A]
					save_json(net,S,C,j)
					j += 1
				end
			end
			i += 1
		end
	end
end

# To access stuff from the simulation, use the following commands

#sim[:B]  # biomasses for every species (columns) for every timestep (rows)
#sim[:p]  # all parameters
#sim[:t]  # list of timesteps

