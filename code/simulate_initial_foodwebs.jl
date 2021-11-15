#######################################################
### 1. Build a food web
### 2. Simulate the dynamics
### 3. Check for extinctions
### 4. If any extinctions, try again 
### 5. Repeat for all combinations of S and C
#######################################################


# Run from command line
# include(“simulations.jl”) 

## Install packages if necessary
#Pkg.add("BioEnergeticFoodWebs")  # install the package
#Pkg.add("JSON")
#Pkg.add("JLD2")
#Pkg.update()  # check that it's the most up to date one

## Load packages
using Pkg
using BioEnergeticFoodWebs  # tell Julia we want to use the package
using LightGraphs
using CSV
using DataFrames

# Save network as json
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
	# We want weakly connected components. Strong means all components are connected along the link direction, so two basals (e.g.) will never be connected.
	# Weakly connected means we ignore the direction.
	weak=LightGraphs.weakly_connected_components(g)

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


# Separating the building of the web and isolation checking from the sim.
# If sim fails to get coexistence, want to send it back to web building
function web_build(S,C,)
	foodweb = nichemodel(S,C); #has predators in rows, and preys in columns. It can only have 0 and 1.

	# Check we don't have any isolated species
	let iter = 1;
		while check_isolated(foodweb) & check_components(foodweb) & (iter < 1000)    # Try 1000 times to get a fully connected network
			foodweb = nichemodel(S,C)
			iter += 1
			# println(iter)
		end
		if iter == 1000
      		error("ERROR: 1000 tries and still getting isolated species")
       	end
	end
	return foodweb
end

# initialize params
maxreps = 10000; # Giving up after this many
reps = 5 # How many of each combination do we want?
S = 10; # how many species do we want?
Smin=50;
Smax=100;
Cmin=0.02;
Cmax=0.2;
		 ## NOTE: C gets converted into links for the function "nichemodel()". If in the conversion the value is < 0, it will throw an error
		 ## Sometimes it can work a bunch of times before we get an error.
start = 0; # where should we start our simulations?
stop = 10; # and how long should they run?
ext_threshold = 0.00001; # what's our threshold for when a species should be considered extinct?


try 
	mkdir("data/networks/pre_disturbance")
catch
end

# Run through the species range
for S in Smin:10:Smax
	try 
		mkdir("data/networks/pre_disturbance/$S/")
	catch
	end
	# And the C range
	for C in Cmin:0.02:Cmax
		println(S," ",C)

		try 
			mkdir("../data/networks/pre_disturbance/$S/$C")
		catch
		end

		j = 0 # Number of good reps
		let i = 0;
			# Run through up to maxreps iteractions
			while (i<maxreps) & (j<reps)
				# Build a niche network with S species and C connectance
				foodweb=web_build(S,C);

				# Use the default parameters
				# EXCEPT! Z=3.065, all non-basal species are vertebrates
				# Z is vertebrate average from Table 1, Brose et al 2006.
				Fwsum = sum(foodweb, dims =2)[1:end,1];
				verts = Bool[ x > 0 for x in Fwsum];
				params = model_parameters(foodweb, vertebrates = verts, Z=Float64(3.065));
				# Set initial biomasses randomly 
				bm = rand(size(foodweb,1));

				# run the simulation
				sim = simulate(params, bm, start = start, stop = stop, use=:nonstiff);
				
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


				# It looks like the best thing for *this* iteration of Julia is to save things as .csv's and read them in individually. All params except web structure are static and don't need to be saved.

				sim[:p][:B]=sim[:B][end,1:end]   # save the final values of biomass
				sim[:B]=collect(sim[:B])
				## Note! If the file already exists, this will fail.
				if check_extinction(sim)
					CSV.write("../data/networks/pre_disturbance/$S/$C/initial_B_$j.csv",DataFrame(sim[:B],:auto),writeheader=true)
					CSV.write("../data/networks/pre_disturbance/$S/$C/initial_t_$j.csv",DataFrame(reshape(sim[:t],length(sim[:t]),1),:auto),writeheader=true)
					CSV.write("../data/networks/pre_disturbance/$S/$C/initial_net_$j.csv",DataFrame(sim[:p][:A],:auto),writeheader=true)

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

