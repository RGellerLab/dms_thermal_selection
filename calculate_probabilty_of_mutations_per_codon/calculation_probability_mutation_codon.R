## calculation of the probability of X mutations ocuring in  the same codon

cds= 6555 # coding region of CVB3
mu= 1e-4 # mutation rate

# probability of each event is:
one=dbinom(x = 1,prob = mu,size = cds) # can happen in any position in the CDS
second=dbinom(x = 1,prob = mu,size = 2)# can happen in the remaining 2 nt in the codon
third=dbinom(x = 1,prob = mu,size = 1) # can happen only in the remaining 1 position

# probabilty of two mutations per codon:
one*second  #  6.807331e-05
# probability of three  mutations per codon:
one*second*third # 6.8e9


# via poisson... similar results
# prob of genome with one mutation
dpois(1, lambda = cds*mu) #0.34

# chance of 2 mutations
dpois(1, lambda = cds*mu)*dpois(1, lambda = 2*mu) #6.8e-5


#chance of 3 mutations
dpois(1, lambda = cds*mu)*dpois(1, lambda = 1*mu)*dpois(1, lambda = 1*mu)# 3.403096e-09 
