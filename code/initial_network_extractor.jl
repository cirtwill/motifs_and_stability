#######################################################
### Extract the initial network from the .jld
### Saves space on github by only having one file per initial web
#######################################################

## Commands for getting set up
using Pkg
#Pkg.add("JSON")
#Pkg.add("JLD")
#Pkg.update()  # check that it's the most up to date one

using JSON
using JLD # Need to use this one to pull up the sim
#Do we want to allow rewiring after extinctions? If so which type?


Smin=50
Smax=100
Cmin=0.02
Cmax=0.2

function save_network(network,S,C,i)
	filename="../data/networks/pre_disturbance/$S/$C/initial_net_$i.json";
	f=open(filename,"w")
	JSON.print(f,network)
	close(f)
end


# Run through the species range
for S in Smin:10:Smax
	# And the C range
	for C in Cmin:0.02:Cmax
		# Run through all the webs
		for i in 0:99
			# Load in the initial web
			dat=load("../data/networks/pre_disturbance/$S/$C/initial_$i.jld")
			sim=dat[:"sim"]

			# Extract the initial food web, parameters
			foodweb = sim[:p][:A]
			save_network(foodweb,S,C,i)
		end
	end
end
