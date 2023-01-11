(*
  operations on the pl File Format from Francois Berenger

  pl stands for "position" "logP contribution (per heavy atom)"

  those files can be created with a patched version of Open Babel
  https://github.com/UnixJunkie/openbabel/tree/logP_contrib_per_heavy_atom
  plus the script ../bin/to_pl_file.sh

  Example from ace_ligdecs.uniq.1conf.pl:
  ---
BEGIN
MOL active_InChIKeyABBSOQIXYPZCKO-NHCYSSNCSA-N
ATOM 8.3894 8.41 -2.139 0.2811
ATOM 8.3573 7.6993 -3.3387 0.2811
ATOM 7.4475 8.1439 -1.1452 0.2811
ATOM 7.3832 6.7224 -3.5448 0.2811
ATOM 6.4733 7.167 -1.3512 0.2811
ATOM 6.4412 6.4563 -2.551 0.136
ATOM 3.8681 3.6178 -1.7283 -0.2783
ATOM 3.9727 0.2529 -3.279 -0.2783
ATOM 4.1412 6.0412 -3.4199 0.5131
ATOM 2.3728 2.0638 -4.1551 0.5131
ATOM 6.2364 3.9058 -0.8576 0.0425
ATOM 5.3767 5.3896 -2.7755 0.2423
ATOM 5.0092 4.6138 -1.4715 0.123
ATOM 3.2404 1.5769 -3.0039 -0.0821
ATOM 4.1734 2.6162 -2.6335 -0.1045
ATOM 5.0817 0.0917 -2.6983 -0.1526
ATOM 2.7886 3.7344 -1.1478 -0.1526
ATOM 3.3627 -0.5286 -4.0661 0.0087
ATOM 5.8469 3.0778 0.719 0.3805
END
[...]
  --- *)

open Batteries
open Printf

module A    = Array
module At   = Atom
module Dacc = Lacc.Dacc
module F    = Filename
module L    = List
module MU   = My_utils
module Mol  = Molecule
module Pqr  = Pqr_parser
module S    = String

let mol_start_str   = "BEGIN"
let mol_start_tag   = Str.regexp_string mol_start_str
let atoms_start_tag = Str.regexp_string "MOL "
let atoms_end_tag   = Str.regexp_string "END"

(* will throw End_of_file once there is no more to read *)
let read_one_molecule counter input =
  let res = Dacc.create () in
  (* skip lines until we start to read a new molecule *)
  while not (Str.string_match mol_start_tag (Legacy.input_line input) 0)
  do () done;
  (* the line after the mol_start_tag one contains the molecule name *)
  let molecule_name = Pqr.extract_name (Legacy.input_line input) in
  (* we are on the atoms now *)
  let rec read_atoms acc =
    let line = Legacy.input_line input in
    if not (Str.string_match atoms_end_tag line 0) then
      let atom = At.of_pl_line line in
      read_atoms (Dacc.accum acc atom)
    else
      ()
  in
  read_atoms res;
  let i = !counter in
  incr counter;
  Mol.create molecule_name i (Dacc.return res)

(* read all molecules from fn
   returns: [(name, [atoms]); ...] *)
let read_molecules fn =
  if not (F.check_suffix fn ".pl")
  then Log.warn "%s not a .pl file" fn;
  let nb_molecules = ref 0 in
  let molecules, _eof =
    MU.with_in_file fn (fun input ->
      MU.unfold_exc (fun () ->
        read_one_molecule nb_molecules input
      )
    )
  in
  Log.info "%d molecule(s) in %s" !nb_molecules fn;
  molecules
