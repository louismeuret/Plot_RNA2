import os
import sys
import numpy as np
import mdtraj as md


def topology_contacts( topology, r0, delta ):
    traj = md.load( topology )
    dist, res_pairs = md.compute_contacts(traj, contacts = "all", scheme = "closest-heavy", ignore_nonprotein = False, periodic = True)

    indexes = [ ]
    for i in range( len( dist[0] ) ):
        if ( dist[0,i] <= r0 and res_pairs[i,1] - res_pairs[i,0] > delta ):
            indexes.append(i)

    return indexes

#-------------------------------------------------------------------------------------

def Q( traj_file, out_file, topology, indexes, r0, delta ):
    traj = md.load( traj_file, top = topology )


    dist, res_pairs = md.compute_contacts( traj, contacts = "all", scheme = "closest-heavy", ignore_nonprotein = TFalse, periodic = False )

    dim = len(indexes)
    out = open(out_file, 'w')
            
    t = 0
    beta= 5.0
    la= 1.5
    
    for i in range( len(dist) ):
        frac = 0
        for j in indexes: #scorro gli indici dei contatti nativi
#            if( dist[i,j] <= r0 ) :
#                frac += 1
            xx = beta*(dist[i,j]-la*r0)
            yy = 1/(1+np.exp(xx))
            frac += yy

        frac /= dim
        string = str(i) + '\t' + str(frac) + '\n'
        out.write(string)

    out.close()
#-------------------------------------------------------------------------------------


if( __name__ == '__main__' ):

    ###############
    # Free params #
    ###############

#    r0 = 0.75
    r0=4.0
    delta = 3

    ###############

    path = input("Path to the trajectories folder (must contain a topology named em.pdb): ").rstrip("/ ")+ "/"

    os.system( "rm -r " + path + "Qs" )
    ll = os.listdir( path )
    os.system( "mkdir " + path + "Qs" )

    ref = path + "em.pdb"

    indexes = topology_contacts( ref, r0, delta )

    for i in ll:

        if ( i.startswith(".") == False and i != "em.pdb" ):
            print( '# Analyzing traj ' + i )
            flnm = path + "Qs/" + i.rstrip(".xtc").rstrip(".trr").rstrip(".pdb") + ".txt"
            Q( path + i, flnm, ref, indexes, r0, delta )
