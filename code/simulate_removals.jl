#######################################################
### Build a food web, simulate the dynamics,
### remove a species, rinse and repeat
#######################################################

# This code saves the biomasses at each time step after simulating the extinction of a species.
# Easier to track extinction order in R or python than Julia.
# Does the extinction threshold work?

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
using JLD # Need to use this one to pull up the sim
using CSV
using DataFrames
#Do we want to allow rewiring after extinctions? If so which type?


# # set working directory
# cd("C:\\Users\\kewo0002\\gitstuff\\motif_role_stability\\data\\Julia\\")

# initialize params
reps = 1
0; # how many replicates do we want?
Smin=50
Smax=100
Cmin=0.02
Cmax=0.2

C = 0.1; # what connectance do we want?
		 ## NOTE: C gets converted into links for the function "nichemodel()". If in the conversion the value is < 0, it will throw an error
		 ## Sometimes it can work a bunch of times before we get an error.
start = 0; # where should we start our simulations?
stop = 50; # and how long should they run?
ext_threshold = 0.00001; # what's our threshold for when a species should be considered extinct?
percent = 10; # What percentage of species can go extinct in a web before we throw it out?
# For complete removal, magnitude of biomass reduction=1
magns = collect(0:0.1:1); # Set which magnitudes of disturbances we want to run (for complete removals set to 1)
trmnt = "bm"; # do we disturb biomass ("bm") or growth rate ("gr")?

# Note: names with colons don't play super well with the terminal.
function save_biomasses(grid,S,C,i,j)
	filename="../data/networks/removals/$S/$C/network_$i:removal_$j.json";
	f=open(filename,"w")
	JSON.print(f,grid)
	close(f)
end

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

# whichsp is column vector of zeros (not disturbed) and ones (disturbed)
# magn is a scalar of the magnitude of the disturbance
# trmnt is a string indicating whether it's biomass ("bm") we're decreasing or growth rate ("gr")
function perturb(whichsp, magn, trmnt, params, B, i, j, S, start, stop, print = true)
	# reset biomass
	if trmnt == "bm"
		removal_bms = zeros(51,S)
		removal_bms[1,1:end] = B[end,1:end]
		bm_dist  = copy(B[size(B)[1],:]);
		# decrease biomass of species specified in whichsp by magn
		bm_dist  = bm_dist.*(1 .- magn * whichsp);
		# Breaking the simulation into smaller chunks so that species which go extinct stay extinct
		# One "step" in simulate = 10 timesteps. Taking only last row to start the next one.
		let timestep = 1;
			while (timestep < 50) # Will be 500 simulated steps total.
				sim_dist = simulate(params, bm_dist, start = start, stop = 1);
				final_bms = sim_dist[:B][end,1:end]
				for q in findall(final_bms .< ext_threshold);
					final_bms[q]=float(0)
				end
				bm_dist = final_bms
				removal_bms[timestep+1,1:end]=final_bms
				timestep += 1
			end
		end
	end

	# # how do I reset growth rate?
	# if trmnt == "gr"
	# 	bm_dist = copy(sim[:B][size(sim[:B])[1],:])
	# 	println("WARNING: Can't change growth rate yet")
	# end

	return removal_bms
end 


# Run through the species range
for S in Smin:10:Smax
	# And the C range
	try 
		mkdir("../data/networks/removals/$S")
	catch
	end
	# And the C range
	for C in Cmin:0.02:Cmax
		println(S," ",C)
		try 
			mkdir("../data/networks/removals/$S/$C")
		catch
		end
		# Run through all the webs
		for i in 0:99
			# Load in the initial web
			B=convert(Matrix,CSV.read("../data/networks/pre_disturbance/$S/$C/initial_B_$i.csv"))
			t=CSV.read("../data/networks/pre_disturbance/$S/$C/initial_t_$i.csv")
			A=CSV.read("../data/networks/pre_disturbance/$S/$C/initial_net_$i.csv")

			# Extract the initial food web, parameters
			foodweb = convert(Matrix,A)
			Fwsum = sum(foodweb, dims =2)[1:end,1]
			verts = Bool[ x > 0 for x in Fwsum]
			params = model_parameters(foodweb, vertebrates = verts, Z=Float64(3.065))

			# Biomasses are the last post-equilibrium biomasses
			bm = B[end,1:end]

			# Don't need to re-simulate.

			# remove each species in turn
			# start with removing no species
			nosp = fill(0,S);
			println("Network ",i)
			for j in 1:S
				#println("Species ",j)
				magn=1
				whichsp  = copy(nosp);
				# set which species to remove
				whichsp[j] = 1;
				# run the disturbance, save the biomasses
				removal_bms = perturb(whichsp, magn, trmnt, params, B, i, j, S, start, stop);
				save_biomasses(removal_bms,S,C,i,j)
			end

		end
	end
end

# To access stuff from the simulation, use the following commands

#sim[:B]  # biomasses for every species (columns) for every timestep (rows)
#sim[:p]  # all parameters
#sim[:t]  # list of timesteps

