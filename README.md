**The new multi-threaded C++ implementation will be available shortly. Thanks for your patience!**

# FindFRs
FindFRs is a Java implementation of the Frequented Regions algorithm presented at the ACM BCB 2017 conference: "Exploring Frequented Regions in Pan-Genoimc Graphs".
A _Frequented Region_ (FR) is a region in a pan-genome de Bruijn graph that is frequently traversed by a subset of the genome paths in the graph.
A path that contributes to an FR being frequent is called a _supporting path_.
The algorithm works by iteratively constructing FRs via hierarchical aglomerative clsutering and then traversing the hierarchy and selecting nodes that qualify as clusters according to the given parameters.

## Parameters
FindFRs has two required parameters: `alpha` and `kappa`.
`alpha` is the minimum fraction of the nodes in an FR that each subpath must contain.
This is referred to as the _penetrance_.
`kappa` is maximum insertion length (measured in base-pairs) that any supporting path may have.
This is referred to as the _maximum insertion_.

Additionally, there are two optional parameters: `minsup` and `minsize`.
`minsup` is the minimum number of genome paths that must meet the other parameters in order for a region to be considered frequent.
This is referred to as the _minimum support_.
`minsize` is the minimum size (measured in de Bruijn nodes) that an FR that meets the other parameters must be in order to be considered frequent.
This is referred to as the _minimum size_.

## de Bruijn Graphs
FrinFRs consumes de Bruijn graphs in the `dot` file format.
A `dot` file representation of a pan-genome de Bruijn graph can be constructed from a `fasta` using the one of the programs presented in the following works:
* "SplitMEM: a graphical algorithm for pan-genome analysis with suffix skips"
* "Efficient Construction of a Compressed de Bruijn Graph for Pan-Genome Analysis"

## Running
FindFRs can be run as follows:
```
    $ java FindFRs <dotFile> <faFile> <K> <alpha> <kappa> [<minsup> <minsize>]
```
where `<dotFile>` is the pan-genome de Bruijn graph constructed for the `<faFile>` with _k_-value `<K>`.
