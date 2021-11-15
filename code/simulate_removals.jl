#######################################################
### 1. Load stable food webs
### 2. Remove a species, simulate dynamics
### 3. Repeat for all species
### 4. Repeat for all webs
#######################################################

# Run from command line 
# include(“simulations.jl”) 

## Install packages if necessary
#Pkg.add("BioEnergeticFoodWebs")  # install the package
#Pkg.add("JSON")
#Pkg.update()  # check that it's the most up to date one

## Load packages
using Pkg
using BioEnergeticFoodWebs  # tell Julia we want to use the package
using JSON
using JLD # Need to use this one to pull up the sim
using CSV
using DataFrames


# initialize params
Smin=50;    # minimum number of species
Smax=100;   # maximum number of species
Cmin=0.02;	# minimum connectance
Cmax=0.2;	# maximum connectance
magn=1;     # how much do we want to decrease biomass by in the disturbance? (1 = full removal)

## NOTE: C gets converted into links for the function "nichemodel()". If in the conversion the value is < 0, it will throw an error
## Sometimes it can work a bunch of times before we get an error.

start = 0; # where should we start our simulations?
stop = 50; # and how long should they run (in steps of 10, so total is stop*10)?
ext_threshold = 0.00001; # what's our threshold for when a species should be considered extinct?


# Note: names with colons don't play super well with the terminal.
function save_biomasses(grid,S,C,i,j)
	filename="../data/networks/removals/$S/$C/network_$i:removal_$j.json";
	f=open(filename,"w")
	JSON.print(f,grid)
	close(f)
end


# whichsp is column vector of zeros (not disturbed) and ones (disturbed)
function perturb(whichsp, magn, params, B, i, j, S, start, stop, print = true)
	# reset biomass

	removal_bms = zeros(stop+1,S)   # matrix for output. Rows are timesteps, columns are species
	removal_bms[1,1:end] = B[end,1:end]  # start biomasses with biomasses from end of initial simulation
	bm_dist  = copy(B[size(B)[1],:]);   
	# decrease biomass of species specified in whichsp by magn
	bm_dist  = bm_dist.*(1 .- magn * whichsp);
	# Breaking the simulation into smaller chunks so that species which go extinct stay extinct
	# One "step" in simulate = 10 timesteps. Taking only last row to start the next one.
	let timestep = 1;
		while (timestep < stop) # Will be 500 simulated steps total.
			sim_dist = simulate(params, bm_dist, start = start, stop = 10);   # run the simulation
			final_bms = sim_dist[:B][end,1:end]                               # what are the biomasses at the end of the simulation?
			for q in findall(final_bms .< ext_threshold);                     # Find those species that went extinct (below the extinction threshold)
				final_bms[q]=float(0)                                         # And set them to zero so we don't get zombies
			end
			bm_dist = final_bms                                               # Final biomasses from the disturbance become the biomasses we start the next round from
			removal_bms[timestep+1,1:end]=final_bms
			timestep += 1
		end
	end

	return removal_bms
end 


# Make directory to save data (may also need to make data/ and networks/)
try 
	mkdir("../data/networks/removals")
catch
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
		for i in 0:4
			# Load in the initial web
			B=Matrix(CSV.read("../data/networks/pre_disturbance/$S/$C/initial_B_$i.csv",DataFrame))      # biomasses
			t=CSV.read("../data/networks/pre_disturbance/$S/$C/initial_t_$i.csv",DataFrame)              # 
			A=CSV.read("../data/networks/pre_disturbance/$S/$C/initial_net_$i.csv",DataFrame)            # network

			# Extract the initial food web, parameters
			foodweb = Matrix(A)
			Fwsum = sum(foodweb, dims =2)[1:end,1]                                                       # How many prey does each species have?
			verts = Bool[ x > 0 for x in Fwsum]                                                          # Set all non-basal species to vertebrates
			params = model_parameters(foodweb, vertebrates = verts, Z=Float64(3.065))                    # calculate the parameters

			# Biomasses are the last post-equilibrium biomasses
			bm = B[end,1:end]

			# remove each species in turn
			println("Network ",i)
			for j in 1:S
				whichsp  = fill(0,S);
				# set which species to remove
				whichsp[j] = 1;
				# run the disturbance, save the biomasses
				removal_bms = perturb(whichsp, magn, params, B, i, j, S, start, stop);
				save_biomasses(removal_bms,S,C,i,j)
			end

		end
	end
end

# To access stuff from the simulation, use the following commands

#sim[:B]  # biomasses for every species (columns) for every timestep (rows)
#sim[:p]  # all parameters
#sim[:t]  # list of timesteps

