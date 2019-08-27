# PRG from MSA

## TL; DR

Rationale: build a sequence graph, starting with one node per aligned character.

The **key parameter** is `num_max_incidents`. Num_incidents is the number of bubbles (fork points)
that end at a certain sequence node. When this number is large, serialisation to a nested PRG string
 is extremely redundant.
 
 Instead, when `num_max_incidents` is reached, the outermost bubble is rewritten, with sequence node
 length doubled. 
 
 Recommended value: 2.

### One Depth form

### Nested form

## Unit tests
Can use strings directly, or files.

* successiveDeletions: is here to illustrate:
    * a case where an empty allele can occur; but this makes genotyping a deletion difficult (using coverage)
    * adjacent markers cannot be avoided: either you do empty alleles, or have site exit and entry next to each other.

## TODOs

- [x] Debugging mode for writing out useful information using boost
- [x] Google test framework
- [ ] Control nesting level
- [ ] Fasta & vcf parsing using ~~htslib~~ [SeqAn](https://seqan.readthedocs.io/en/master/index.html)
- [ ] Dot graph production possibility- see vg.
- [x] Use ~~unique~~ shared pointers for `auto_Node`s
- [ ] Prg string serialisation: probably assign site numbers based on topological ordering + 
add the site number to the bubble start node (map of site nb -> auto_Node).
- [ ] Indel-tied haplotype expansion + 
make a test case of multiple deletions inside large incidence bubble.
- [ ] Make a .gfa of built graph
- [ ] Commit leading nucleotide to Is it possible to always avoid empty site (indel). 
- [ ] Consider enumerating (DFS-ing) all paths within small bubbles; 
gramtools will appreciate this (concurrent allele querying). Cf test `SNP_rewriting`