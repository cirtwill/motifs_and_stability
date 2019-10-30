import pymfinder as py
import networkx as nx
import sys
import os

# List of positions identified as (motif, npreds, npreys)
all_positions=py.roles.STOUFFER_ROLE_IDS.keys()

# Read in an interaction list network
def web_builder(webfile):
  predprey=[]
  f=open(webfile)
  for line in f:
    (pred,prey)=line.split()
    predprey.append((pred,prey))
  f.close()

  spp=[]
  for (pred,prey) in predprey:
    spp.append(pred)
    spp.append(prey)
  spp_sorted=sorted(spp)   #list of species
  shortlist=[]             #new list with each species entered only once
  for sp in spp_sorted:
    if not sp in shortlist:
      shortlist.append(sp) 
  G=nx.DiGraph() #Directed graph. To undirect nx.Graph()
  for sp in shortlist:
    G.add_node(sp) 
  for pair in predprey:
    G.add_edge(*pair)  
  return G

# Extract all triplets in web
# Do we allow cannibalism? Can we prevent it? That would be simpler.
# It appears that cannibalism is allowed. Will need to double-check that this doesn't cause havoc with pymfinder.
# Seems to be safely ignored in one of the small test-drive webs :)
def extract_triplets(G):
	triplets = []
	for s1 in set(G.nodes()):
		partners=set(G.predecessors(s1)).union(set(G.successors(s1))) # Anyone who interacts with s1
		for s2 in partners:
			if s2!=s1:
				partners2=set(G.predecessors(s2)).union(set(G.successors(s2)))
				for s3 in partners2:
					if s3 not in [s1,s2]:
						triplist=sorted((s1,s2,s3))
						if triplist not in triplets:
							triplets.append(triplist)

	return triplets

# Extract a triplet and write to tiny.web
def create_tiny_web(G,triplet):
	links=[]
	for (s1,s2) in set(G.edges()):
		if s1 in triplet:
			if s2 in triplet:
				links.append((s1,s2))

	f=open('tiny.web','w')
	for (s1,s2) in links:
		f.write(s1+'\t'+s2+'\n')
	f.close()

# Get the motif role of each species in tiny.web
# Returns a dict of motif: position
def extract_pymfinder_output(res):
	roledict={}
	for m in res.nodes:
		name=str(res.nodes[m].id)
		role=sorted(res.nodes[m].roles.items())
		for (position, val) in role:
			if val==1:
				roledict[name]=position
	return roledict

# Fold back into a dict for the whole web.
def get_network_roles(net):
	G=web_builder(net)
	network_roles={}
	for sp in set(G.nodes()):
		network_roles[sp]={}
	for trip in extract_triplets(G):
		create_tiny_web(G,trip)
		res=py.pymfinder('tiny.web',motifsize=3)
		roledict=extract_pymfinder_output(res)
		for sp in roledict:
			try:
				network_roles[sp][roledict[sp]]+=1
			except KeyError:
				network_roles[sp][roledict[sp]]=1
	return network_roles

# Aaaaand ... now what?

def main():
	netdir='../../data/networks/'
	for net in os.listdir(netdir):
		print net
		network_roles=get_network_roles(netdir+net)
		print 'roles done'

		
if __name__ == '__main__':
  main()
