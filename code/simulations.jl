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
#Pkg.update()  # check that it's the most up to date one

using BioEnergeticFoodWebs  # tell Julia we want to use the package
using JSON
#Do we want to allow rewiring after extinctions? If so which type?


# # set working directory
# cd("C:\\Users\\kewo0002\\gitstuff\\motif_role_stability\\data\\Julia\\")

# initialize params
reps = 20; # how many replicates do we want?
S = 50; # how many species do we want?
C = 0.1; # what connectance do we want?
		 ## NOTE: C gets converted into links for the function "nichemodel()". If in the conversion the value is < 0, it will throw an error
		 ## Sometimes it can work a bunch of times before we get an error.
start = 0; # where should we start our simulations?
stop = 10; # and how long should they run?
ext_threshold = 0.00001; # what's our threshold for when a species should be considered extinct?
percent = 10; # What percentage of species can go extinct in a web before we throw it out?
magns = collect(0:0.1:1); # Set which magnitudes of disturbances we want to run (for complete removals set to 1)
trmnt = "bm"; # do we disturb biomass ("bm") or growth rate ("gr")?


function save_json(file,i,j,S,trmnt)
	filename = "fw_dstrbnc/S_$S/fw_$i/$(trmnt)_$j.json";
	f=open(filename,"w");
	JSON.print(f,file);
	close(f)
end

# check that we don't have any species that are neither preyed upon nor prey apon
# returns false if there are NO isolated species and true if there are TRUlY isolated species
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

# Make the food web and check that we don't have isolated species
# Repeat up to 100 times to get a food web that has no isolated species
function build_foodweb(S,C)
	# Build a niche network with S species and C connectance
	foodweb = nichemodel(S,C); #has predators in rows, and preys in columns. It can only have 0 and 1.
	#
	#Check we don't have any isolated species
	let iter2 = 1;
		while check_isolated(foodweb) & (iter2 < 100)
			foodweb = nichemodel(S,C)
			iter2 += 1;
			println(iter2)
		end
		if iter2 == 100
      		error("ERROR: 100 tries and still getting isolated species")
       	end
	end

	return foodweb;
end

# Check if any species went extinct
function check_extinctions(sim, ext_threshold, percent)
	last = size(sim[:B])[1]
	#
	println(sim[:B][last,:])
	#
	#if maximum(sim[:B][last,:].<ext_threshold)
	#what happens if we let some go extinct?
	if sum(sim[:B][last,:].<ext_threshold) > (percent/10*S)
		rerun = true;
	else 
		rerun = false;
	end

	return rerun;
end

# whichsp is column vector of zeros (not disturbed) and ones (disturbed)
# magn is a scalar of the magnitude of the disturbance
# trmnt is a string indicating whether it's biomass ("bm") we're decreasing or growth rate ("gr")
function perturb(whichsp, magn, trmnt, sim, params, i, j, S, start, stop, print = true)
	
	# reset biomass
	if trmnt == "bm"
		bm_dist  = copy(sim[:B][size(sim[:B])[1],:]);
		# decrease biomass of species specified in whichsp by magn
		bm_dist  = bm_dist.*(1 .- magn * whichsp);
		sim_dist = simulate(params, bm_dist, start = start, stop = stop);
	end

	if print
		# need to figure out what to do about j (can't handle a vector)
		save_json(sim_dist, i, j, S, trmnt);
	end


	# how do I reset growth rate?
	if trmnt == "gr"
		bm_dist = copy(sim[:B][size(sim[:B])[1],:])
		println("WARNING: Can't change growth rate yet")
	end

	return sim_dist
end 


try 
	mkdir("fw_dstrbnc/S_$S")
catch
end


# Run through as many iterations as we want
for i in 1:reps
	try 
		mkdir("fw_dstrbnc/S_$S/fw_$i")
	catch
	end

	rerun = true;
	let iter1 = 1;
		while (rerun == true) & (iter1 < 5000)
			iter1 += 1;
			# set j (for use in saving network)
			j = 0;

			# Get a foodweb and check for isolated species
			foodweb = build_foodweb(S,C);

			# Use the default parameters 
			params = model_parameters(foodweb);
			# Set initial biomasses randomly 
			bm = rand(size(foodweb,1));
			#bm = fill(1,S)
			# run the simulation
			sim = simulate(params, bm, start = start, stop = stop);
			#
			rerun = check_extinctions(sim, ext_threshold, percent);
		end
		if iter1 == 5000
			error("ERROR: 5000 tries and still getting species extinctions. Rep = $i")
		end
	end

	# save the output - this will be fixed in the next release of BioEnergeticFoodWebs, for now a work-around
	save_json(sim,i,j,S,"initial")

	#BioEnergeticFoodWebs.save(filename = "init_fw_"*string(i), as = json)
	#JLD.save(filename,"sim",sim)

	# start with removing no species
	nosp = fill(0,S);

	# run through all magnitude levels:
	for magn in magns
		# to perturb one by one:
		for j in 1:S
			whichsp  = copy(nosp);
			# set which species to remove
			whichsp[j] = 1;
			# run the disturbance
			sim_dist = perturb(whichsp, magn, trmnt, sim, params, i, j, S, start, stop);
		end

		# or all at once:
		whichsp  = fill(1,S);
		sim_dist = perturb(whichsp, magn, trmnt, sim, params, i, j = "all", S, start, stop);

	end




end


# To access stuff from the simulation, use the following commands

#sim[:B]  # biomasses for every species (columns) for every timestep (rows)
#sim[:p]  # all parameters
#sim[:t]  # list of timesteps

