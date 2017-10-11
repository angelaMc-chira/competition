# Functions for converting ape-format phyogenies to
# 'peth': a format better suited for time-forward simulation.
# magnusclarke@gmail.com
# modified 2016

library(ape)
library(TESS)

#--------------------------------------------------------------------------------------#
#---- Random ultrametric tree ---------------------------------------------------------#
#--------------------------------------------------------------------------------------#
rand_umt    = function(nt, lambda=1, mu=0, max=1e5)
{
    # TESS versions 1 and 2 have different function names and parameters.
    if(unlist (packageVersion('TESS'))[1] == 1)
    {
        tree = sim.globalBiDe.taxa(n=1, nTaxa=nt, lambda=lambda, mu=mu)
    } else {
        tree = tess.sim.taxa(n=1, nTaxa=nt, max=max, lambda=lambda, mu=mu)
    }
    return( tree[[1]] )
}

#--------------------------------------------------------------------------------------#
#---- Constructor for peth tree class. ------------------------------------------------#
#--------------------------------------------------------------------------------------#
pethtree = function(map, dord, nts, tts) 
{
    out = list(map, dord, nts, tts)
    class(out) = 'pethtree'
    names(out) = c('map', 'data_order', 'splitting_nodes', 'times')
    invisible(out)
}

#--------------------------------------------------------------------------------------#
# Convert a tree in ape format to one in new format (return a new 'peth' tree object). #
# Only works for ntip > 4. ------------------------------------------------------------#
#--------------------------------------------------------------------------------------#
ape2peth = function(tree)
{
    edge = tree$edge
    edge.length = tree$edge.length

    nedge = length(edge[,1])
    ntip = (nedge+2)/2

    # Set up ape-to-peth map and vector to keep track of running branches.
    running = rep(FALSE, 2*ntip-2)
    map = matrix(ncol=2, nrow=nedge+1)

    # Assign root
    root = edge[1,1]
    map[1,1] = root
    map[1,2] = 1

    # splitting gives the NODE, as labelled in edge, which is speciating
    splitting = root

    edge_counter = 1

    # Keep track of how many nodes we've already assigned.
    node_counter = 1

    # Produce vector of nodes to split for future simulations (peth numbering)
    nts = c()

    # Produce vector of times between speciation events
    tts = c()

    # Loop over speciation events
    for(i in seq(2,2*(ntip-1), by=2))
    {
        # 'splitting' is the splitting node in ape format
        # 'peth_split' is the splitting node in peth format
        peth_split = map[ which(map[,1] == splitting), 2 ][1]
        nts = c(nts, peth_split)

        # The branches coming from this node are now running
        running[ which(edge[,1]==splitting) ] = TRUE
        running[ which(edge[,2]==splitting) ] = FALSE

        # We're at a speciation event. Find out which nodes come from it.
        new_edges = which(edge[,1] == splitting)
        new_nodes = edge[new_edges, 2]

        # Get which of the new nodes comes from the shortest branch
        # We want this one to inherit the label-number of its parent.
        new_node1 = new_nodes[order(edge.length[new_edges])[1]]
        new_node2 = new_nodes[order(edge.length[new_edges])[2]]

        # First new node keeps ancestor number (splitting node), 
        # other new node takes node_counter+1
        map[i,1] = new_node1
        map[i,2] = peth_split
        map[i+1,1] = new_node2
        map[i+1,2] = node_counter + 1

        node_counter = node_counter + 1

        # time_to_split = time to next split = shortest branch running.
        time_to_split = min(edge.length[running])
        tts = c(tts, time_to_split)

        # Subtract time to next split from running branches.
        edge.length[running] = edge.length[running] - time_to_split

        # Which node splits next? 
        # Of the running branches, take the shortest one and then the splitting node 
        # is the second elemtn of edge for that branch (i.e. the child)
        el_temp = edge.length
        el_temp[!running] = 9e99
        splitting = edge[ which.min(el_temp),2 ]
    }
    
    # Need to know what order simulated data will be in, relative to ape's expectatons
    ordered_map = order(map[,1])
    want = ordered_map[1:ntip]
    dord=rev(map[want,2])

    # Need to know what order simulated data will be in, relative to ape's expectatons
    # This is given by working from bottom up through map
    dord = 1:ntip
    for(i in 1:ntip)
    {
        want = which(map[,1]==i)
        dord[i] = map[want , 2]
    }

    ptree = pethtree(map, dord, nts, tts)
    return(ptree)
}

#--------------------------------------------------------------------------------------#
#----- Convert ape to peth, or pass on peth tree, or return error if incorrect. -------#
#--------------------------------------------------------------------------------------#
checktree = function(t)
{
    if(class(t)=="phylo")
    {
        t = ape2peth(t)
    } else if(class(t)!="pethtree") {
        stop("Tree incorrectly formatted.")
    }
    t
}
