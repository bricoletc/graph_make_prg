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
- [] Control nesting level
- [] Fasta & vcf parsing using ~~htslib~~ [SeqAn](https://seqan.readthedocs.io/en/master/index.html)
- [] Dot graph production possibility- see vg.
- [x] Use ~~unique~~ shared pointers for `auto_Node`s
- [] Prg string serialisation: probably assign site numbers based on topological ordering + 
add the site number to the bubble start node (map of site nb -> auto_Node).
- [] Can have a large inc. fixed point inside another: delete the former, using uncleared 
fixed_point_numbers map.
- [] Indel-tied haplotype expansion + 
make a test case of multiple deletions inside large incidence bubble.
