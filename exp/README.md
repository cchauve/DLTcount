The directory contains the results dicussed in the paper, together with the code to obtain these results.

To repeat the experiments:
> cd scripts
> run_exp1.sh # Generate random species trees, results in the directory trees/
> run_exp2.sh # Count DL and DLT histories for the generated trees, considered as unranked, results in unranked/
> run_exp3.sh # Count DL and DLT histories for the generated trees eqch with 10 different random rankings, results in ranked/
> run_exp4.sh # Compute the asymototics growth factor for the generated trees, in the DL model
> run_exp5.sh # Sample 10000 histories of size 30, for 50 random species trees, in the DL and DLT models, both unranked and ranked (random ranking), results in sampling/