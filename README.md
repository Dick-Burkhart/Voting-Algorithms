My repository (https://github.com/Dick-Burkhart/Voting-Algorithms) contains many Fortran routines (*.f90) that are used by 3 projects for the paper "A Clustering Algorithm for Proportional Representation"

My code is in Fortran 2003, using the 2022 Intel classic Fortran compiler.

The first project does preprocessing (Pre_process.f90) of the ballots from 40 electoral districts. Run input parameters in Pre_opt.txt, with data parameters Vote_files.txt. The actual ballots for each district are in the file with the district name.

The second project computes the proportional representation for these 40 elections by clustering (PR_clustering.f90). Run input parameters in PR_opt.txt, with data parameters in Data.txt, RndSeed.txt, and Vote_files.txt. The ballots for each district are in the file with the district name followed by 1 for strong ranking, 2 for weak ranking, and 3 for rating data. These were produced by Pre_process.f90.

The third project tests the sensitivity of this clustering to changes in the domain of candidates (Domain_sens.f90). Run input parameters in Domain_opt.txt, with data parameters in Data.txt, Top_dat1.txt, and RndSeed.txt. The ballots for each district are those produced by Pre_process.f90.

Note: The run input parameter files always include the name of the print out file, the print out level, and the range of districts. The Vote_files.txt contains the electoral data for each district except that the actual ballots are in district named files


First project overview:

The preprocessing on district ballot files transforms them into a standard format 
for ranking or rating.  It also may perform ballot consolidation, both to reduce the computational load by eliminating marginal candidates and to restrict the number of candidates ranked or rated by each ballot if that is specified by the user, or to reduce the computation.  

For testing weak ranking in PR_clustering.f90, the strongly ranked ballots are perturbed so that sometimes more than 1 ballot may be assigned the same ranking.  Likewise, these strongly ranked ballots may be perturbed to produce ballots with ratings, given specified rating levels with point values from 10.0 to -10.0, which are adjusted to yield ballots whose rating point values average out close to 0.

    
Second project overview:

The clusters are the voting blocks needed by PR.  The number of candidates that represent a voting block should equal the size of the voting block, where the total vote size is scaled to the number of candidates to be elected.  So if this number is an integer, the proportionality can only be approximate.  In addition, if the voting block itself is not predefined, then it can only be known after the clustering how well a particular candidate represents a particular voting block, with much ambiguity possible.

To make this problem solvable, at least in an optimization sense, we compute how well the voters in a voting block rank or rate each candidate.  This is done by computing, for each cluster, a mean ranking or rating vector over the candidates, by specifying a function which converts ranking or rating levels to point values (real numbers normalized to a maximum value of 10.0 with a value of 0.0 for unranked candidates or neutrally rated candidates).

These mean point vectors are also used to distinguish the clusters, since strongly correlated mean vectors indicate strongly overlapped clusters, which should be merged.  In addition, clusters that represent too few voters are deleted by merging them into a group of "independents", along with unclustered ballots.  This leaves the "regular clusters" (after mering and deleting) to represent the voting blocks. 

Each ballot may be partially represented by more than one regular cluster, depending on the correlation of that ballot's point vector with the mean point vector of the cluster.  The resulting cluster memberships of each ballot are limited so that their sum over the regular clusters does not exceed 1.0, and if it is < 1.0 then the shortfall in membership is assigned to the independents.  Thus the size of a voting block is the sum of memberships in its cluster over all ballots, with any shortfall in the sum of these voting block sizes going to the independents. 

This way the group of independents has its own mean point vector so that it gives at least a minimal representation to each candidate, with this representation possibly significant if only 1 or 2 voting blocks have been identified.  When this happens the independents contribute to the determination of proportionality.

The key fact here is that the determination of the regular clusters is non-linear process since the mean point vector of a cluster is computed from the memberships of the ballots in that cluster, which are in turn revised by the correlations of the ballot mean vectors with the new cluster mean vector.  This iterative process normally converges (called k-means in the clustering literature), but it may require a damped Newton method in certain cases, or may fail due to stagnation or divergence in a few cases. 

In addition this whole process requires an initial set of clusters to start it off.  In practice this could be some simple method of approximate clustering, or some randomized guesses.  For each electoral district, we start over 20 times 15 of them random initial cluster sets. However all this is still too much computation, even for ballots that have been consolidated in pre-processing.  Thus we do a further consolidation to 100 or so "slate ballots" (still a partition), centered around slates of 1, 2, or 3 top ranked candidates with sufficient summed ballot weight to form "slate clusters".

The 20 converged cluster sets may not be all the same, so an objective function is needed to evaluate them.  The key determinant is the number of clusters in a converged cluster set, but penalty functions allow for a nuanced evaluation.  Thus a small cluster picks up a size penalty as its size approaches the deletion level. Likewise overlapping clusters pick up an overlap penalty as their strongest correlation approaches the merge level. Otherwise the objective function tends to maximize the sum of the sizes of the clusters in the set.

Likewise an objective function is needed to evaluate a possible set of elected candidates for a given cluster set.  To compute the proportionality, we already have the sizes of the clusters in the cluster set, but still need to figure out how well a set of candidates matches a cluster size, based on the mean values of a cluster for the elected candidates.  This is given by our concept of a "portion", which is simply the product of the candidate's mean cluster value times the size of the cluster, normalized so that the sum over all the clusters, plus the independents, is 1.0. 

This number represents the portion, or fraction, of that candidate's voting power which derives from that cluster (the higher the rating and bigger the size of the cluster, the stronger its representation for that candidate).  Thus good proportionality means that the sum of the portions over the elected candidates for each cluster should closely match the size of that cluster.  Our penalties begin only when this mismatch exceeds 1/5 of a candidate. In the absence of proportionality penalties, our objective is to compute the cluster-size averaged mean value (including independents) for each candidate, then sum 
over the elected candidates, much like the Borda Count.  All possible sets of  elected are examined for each of the top tier of cluster sets if there is more than one.

Third project overview:

This program implements uses PR_clustering to test the sensitivity of clustering and non-clustering methods to perturbations or changes in the domain as a subset of the full set of candidates 1...nc. This is related  to Nic Tideman's concept for a pairwise tournament between competing sets of candidates, with winning set to be elected. Except that my clustering full objective functions are used to compare the pairs of sets instead of STV concepts. 

The domains may be either subsets of size nc – 2 or subsets of np+1 through 2*np, where np candidates are to be elected.  The second method is like Tideman’s Tournament algorithm, where 2 competing sets of size np are compared on the domain consisting of their union. 

The point of a sensitivity test is that the best elected set may differ when it is computed for a domain of candidates that is a proper subset of the full set because the proportionality may be affected.  To demonstrate this I compute the full objective function for the elected set on the domain ‘D’ in question. 

To compare 2 subsets A’ and ‘B’ of ‘D’ I use the ratio Rd(A,B) = Fd(A) / Fd(B) of their full objective functions Fd for the top cluster set of ‘D’. Then the overall objective function for ‘A’ becomes a weighted average Sd(A) of a function of Rd(A,B) over all distinct sets ‘B’ over all domains ‘D’ containing both ‘A’ and B’, weighted by the cluster set objective for the top cluster set of ‘D’

However, it is surprisingly effective to just use N(A) = number of domains for which ‘A’ is the top elected set. In fact, to determine the “winner”, we  restrict our search to sets ‘A’ which are the top elected set (= the set with the best fitness) for at least one domain. Then Sd(A) could be used as a tie breaker.
