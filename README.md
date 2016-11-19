mrpmatrix
=========

This java program builds, quickly, 
a matrix representation (MR) on a set of (input) source trees. 

Each of the (input) source trees is expected to be represented 
in Newick tree format.

It uses some tricks to use very little time and memory.


```
Usage: <trees_file> <output> <output_format> [-dna] [-randomize seed]
                <trees_file>: A file containing Newick trees, one tree per line
                <output>: The name of the output matrix representation (MR) file
                <output_format>: use NEXUS for nexus, PHYLIP for phylip, or FASTA for fasta fromatted output
                -dna: output As and Ts instead of 0 and 1
                -randomize: randomize 0-1 codings, the seed number is optional
```

The randomize option is useful when doing MRL analyses, but is not needed for MRP. 
