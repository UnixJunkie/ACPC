
(* operations on the Tripos Mol2 File Format:
   - extract molecule names and
     pqr information (atom positions and partial charges)
   Reference document:
   Tripos Mol2 File Format, SYBYL 7.1 (Mid-2005),
   http://www.tripos.com/data/support/mol2.pdf *)

open Batteries

module A  = Array
module At = Atom
module L  = List
module Log = Dolog.Log
module F  = Filename
module MU = My_utils

type molecule_name = string

type parser_state =
  | Looking_for_molecule
  | Reading_atoms of (molecule_name * int * Atom.atom list)

let mol_start_str   = "@<TRIPOS>MOLECULE"
let mol_start_tag   = Str.regexp_string mol_start_str
let atoms_start_tag = Str.regexp_string "@<TRIPOS>ATOM"
let atoms_end_tag   = Str.regexp_string "@<TRIPOS>BOND"

(* read all molecules from fn
   returns: [(name, [atoms]); ...]
 *)
let read_molecules fn =
  if not (F.check_suffix fn ".mol2") then (
    Log.warn "%s not a .mol2 file" fn
  );
  let nb_molecules = ref 0 in
  let res = ref [] in
  let state = ref Looking_for_molecule in
  let input = open_in fn in
  (try while true do
      let line = input_line input in
      (match !state with
        | Looking_for_molecule ->
          if Str.string_match mol_start_tag line 0 then (
            (* read molecule's name *)
            let mol_name = input_line input in
            state := Reading_atoms (mol_name, !nb_molecules, []);
            incr nb_molecules;
            (* skip until list of atoms *)
            while not (Str.string_match atoms_start_tag (input_line input) 0)
            do () done);
        | Reading_atoms (mol_name, index, atoms) ->
          (* detect end of atoms list *)
          if Str.string_match atoms_end_tag line 0 then (
            (* preserve atoms' order from the molecule *)
            res := (mol_name, index, L.rev atoms) :: !res;
            state := Looking_for_molecule)
          else (* read atom *) (
            let atom = At.atom_of_mol2_line line in
            state := Reading_atoms (mol_name, index, atom :: atoms)));
    done with _ -> close_in input);
  Log.info "%d molecule(s) in %s" !nb_molecules fn;
  (* preserve molecules' order from the input file *)
  L.rev !res

(* store all lines of all molecules, for later reference by index *)
let read_raw_molecules fn =
  Log.info "indexing molecules from %s..." fn;
  if not (F.check_suffix fn ".mol2") then (
    Log.warn "%s not a .mol2 file" fn
  );
  let all_lines = MU.string_list_of_file fn in
  let all_molecules =
    L.tl (* nsplit has a strange semantic: the first list it returns
            is either empty or uninteresting for us on MOL2 files *)
      (L.nsplit
         (fun l -> Str.string_match mol_start_tag l 0)
         all_lines)
  in
  let nb_molecules = L.length all_molecules in
  (* put back the header of each molecule, removed by BatList.nsplit *)
  let res = A.of_list (L.map (fun m -> mol_start_str :: m) all_molecules) in
  Log.info "%d molecule(s) in %s" nb_molecules fn;
  res

exception Break

(* will throw End_of_file once there is no more to read *)
let read_one_molecule input =
  let res = ref [] in
  (try
     while true do (* skip lines until we start to read a new molecule *)
       let line = input_line input in
       if Str.string_match mol_start_tag line 0 then raise Break
     done
   with Break -> ());
  (* the line after the mol_start_tag is the molecule name in MOL2 format *)
  let molecule_name = input_line input in
  (try
     while true do (* skip lines until we start to read atoms *)
       let line = input_line input in
       if Str.string_match atoms_start_tag line 0
       then raise Break
     done
   with Break -> ());
  (try
     while true do (* read all atoms *)
       let line = input_line input in
       if Str.string_match atoms_end_tag line 0
       then raise Break
       else
         let atom = At.atom_of_mol2_line line in
         res := atom :: !res
     done
   with Break -> ());
  (molecule_name, L.rev !res)
