ACPC v1.0
=========

Don't hesitate to contact Francois Berenger in case you have problems running
the software or discover a bug.

Installation
------------

To install ACPC, you first need to install and configure
the OCaml package manager (OPAM). Cf. http://opam.ocaml.org/

Once this is done, you can automatically install ACPC
with the following command:

$ opam install acpc

Recommended usage (protocol)
----------------------------

ACPC was designed to be rotation and translation invariant.
ACPC is not invariant to the conformer of a molecule,
neither to the charge model that was used to assign partial charges to it.

The following protocol has been validated:
the query molecule(s) _AND_ the database to screen must be prepared in
the same way. The same software with same parameters must be used
to assign partial charges and generate conformers for _ALL_ molecules.

The recommended charge models are:
MOE's MMFF94x or as a fallback Open Babel's Gasteiger.

Reading the related research article is highly recommended:
"A rotation-translation invariant molecular descriptor of
partial charges and its use in ligand-based virtual screening".
Citing the article is kindly asked from users of the software.

Examples
--------

1) one query on a small database

$ acpc -q query.mol2 -db database.mol2

2) same but storing the top 10 molecules

$ acpc -q query.mol2 -db database.mol2 -top 10 -o ten_best.mol2

3) one query on a large database

$ acpc_big -q query.mol2 -db database.mol2

4) separate each molecule from a mol2 file into separate files

$ acpc_mol2tool some_molecules.mol2

Get some help
-------------

$ acpc -h
  -cmp {CC|Tani|Tref|Tdb} LBAC+/- comparison method (default: CC)
  -htq                    list molecules scoring Higher Than the Query with itself
  -q query.mol2           query (incompatible with -qf)
  -qf f                   file containing a list of mol2 files (incompatible with -q)
  -db db.mol2             database
  -dx float               X axis discretization (default: 0.005000)
  -v                      output intermediate results
  -nopp                   don't rm duplicate molecules
  -np nprocs              max CPUs to use (default: 1)
  -ng                     no gnuplot
  -nr                     no ROC curve (also sets -ng)
  -o output.mol2          output file (also requires -top, incompatible with -qf)
  -top N                  nb. best scoring molecules to output (also requires -o)
  -help                   Display this list of options
  --help                  Display this list of options

# acpc_big -h
  -q query.mol2 query
  -db db.mol2   database
  -dx float     X axis discretization (default: 0.005000)
  -help         Display this list of options
  --help        Display this list of options

WARNING
-------

Don't do two queries at the same time on the same computer or on top of NFS
with a same query molecule file (-q SOME_QUERY.mol2), this may overwrite
result files in a strange way (SOME_QUERY.ranks, SOME_QUERY.scores and
SOME_QUERY.scored-label).
