# PRG from MSA

### One Depth form

### Nested form

## Test data
Can use strings directly, or files.

### adjacent_markers.fasta

Is here to illustrate a case where adjacent markers cannot be avoided. 
A choice needs to be made between adjacent variant site markers vs adjacent site and allele markers (ie an empty allele record).

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