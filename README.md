# PRG from MSA

## TL; DR

Rationale: build a sequence graph, starting with one node per aligned character. This initially allows for as
 much recombination as possible.

Then areas of the graph where this approach is not viable are *haplotype expanded*: they are rebuilt, with nodes
having sequence size that gets doubled.

Define 'not viable' as :

    * No bubble can be traversed & used more than once **within another bubble**. This preserved variant site uniqueness.
    * There is a limit on how many bubbles can end at the same sequence node: parameter `num_max_incidents`

When `num_max_incidents` is reached, the outermost bubble is rewritten, with sequence node
 length doubled. Current value: 2.

### One Depth form

### Nested form

## Unit tests
Can use strings directly, or files.

* successiveDeletions: is here to illustrate:
    * a case where an empty allele can occur; but this makes genotyping a deletion difficult (using coverage)
    * adjacent markers cannot be avoided: either you do empty alleles, or have site exit and entry next to each other.

## TODOs

- [ ] Control nesting level
- [ ] Fasta & vcf parsing using ~~htslib~~ [SeqAn](https://seqan.readthedocs.io/en/master/index.html)
- [ ] Dot graph production possibility- see vg.
- [ ] Make a .gfa of built graph
- [ ] Consider enumerating (DFS-ing) all paths within small bubbles; 
gramtools will appreciate this (concurrent allele querying). Cf test `SNP_rewriting`
- [ ] Prg string serialisation: probably assign site numbers based on topological ordering + 
- [ ] Refactor graph code into smaller functions
- [ ] ^^Unit test those functions^^
- [x] Debugging mode for writing out useful information using boost
- [x] Google test framework
add the site number to the bubble start node (map of site nb -> auto_Node).
- [x] Use ~~unique~~ shared pointers for `auto_Node`s
- [x] Indel-tied haplotype expansion + 
make a test case of multiple deletions inside large incidence bubble.
- [x] Commit leading nucleotide to Is it possible to always avoid empty site (indel). 
