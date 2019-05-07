# Counting and sampling gene family evolutionary histories in the duplication-loss and duplication-loss-transfer models

### Cedric Chauve, Yann Ponty, Michael Wallner, May 2, 2019.

The directory contains the results dicussed in the paper, together with the scripts to obtain these results.

The directory *scripts/* contains the python and shell scripts to repeat the experiments. To repeat the experiments:  
> \> cd scripts  
> \> run_exp1.sh # Generate random species trees, results in the directory trees/  
> \> run_exp2.sh # Count DL and DLT histories for the generated trees, considered as unranked, results in unranked/  
> \> run_exp3.sh # Count DL and DLT histories for the generated trees each with 10 different random rankings, results in ranked/  
> \> run_exp4.sh # Compute the asymototics growth factor for the generated trees, in the DL model  
> \> run_exp5.sh # Sample 10000 histories of size 30, for 50 random species trees, in the DL and DLT models, both unranked and ranked (random ranking), results in sampling/  

The directory *trees/* contains the 100 unranked species trees generated for experiments in unranked models; the trees of size $k$ are in the file *trees_unranked_k*. It also contains the 10 random rankings generated for each such trees; the trees and rankings for species trees of size $k$ are in the file *trees_ranked_k*.

The directory *unranked/* contains the results of counting the number of histories of size up to $50$ for each unranked species tree; the results for the species trees of size $k$ are in the file *results_unranked_k.gz*.

The directory *ranked/* contains the results of counting the number of histories of size up to $50$ for each ranked species tree; the results for the species trees of size $k$ are in the file *results_ranked_k.gz*.

The directory *asymptotics/* contains the exponential growth factor for each unranked tree in the unranked DL-model; the results for species trees of size *k* are in the file *asymptotics_k_DL*.

Last the directory *sampling/* contains the results of sampling histories of size $n=30$ in $50$ random species trees of size $k=16$ in the four models unranked DL (UDL), unranked DLT (UDLT), ranked DL (RDL) and ranked DLT (RDLT). For each sampled history, we record the number of genes per extant species, the number of losses, duplications and HGTs.
