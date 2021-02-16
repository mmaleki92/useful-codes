def iterateWolff(mat,T,calcE=False):
    # Wolff algorithm for a square lattice
    """
    Parameters:
        mat: numpy array containing the current spin configuration
        T: Temperature for the simulation
        calcE: Flag to calculate and return the average spin energy for the final configuration.
    """
    L,_ = mat.shape
    tracker = np.zeros((L,L)) # Keep track of which spins have already be added to the cluster.
    
    i,j = np.random.randint(0,L,size=2)
    spin = mat[i,j]
    stack = [(i,j)]
    tracker[i,j]=1
    
    cluster = [(i,j)]
    while len(stack)>0:
        i,j = stack.pop()
        neighbors = [(i,(j+1)%L),(i,(j-1)%L),((i+1)%L,j),((i-1)%L,j)]
        for pair in neighbors:
            l,m = pair
            if (mat[l,m]==spin and tracker[l,m]==0 and np.random.random()< (1.0-np.exp(-2.0/T))):
                cluster.append((l,m))
                stack.append((l,m))
                tracker[l,m]=1
            
    flipCluster(mat,cluster)
            
    if calcE:
        avgE=0.0
        for i in range(L):
            for j in range(L):
                spin_final = mat[i,j]
                neighbor_sum = mat[(i+1)%L,j]+mat[(i-1)%L,j]+mat[i,(j+1)%L]+mat[i,(j-1)%L]
                E_final = -spin_final*neighbor_sum
                avgE+=float(E_final)
    
        avgE/=float(L**2)                   
        return avgE
