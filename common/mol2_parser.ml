
(* operations on the Tripos Mol2 File Format:
   - extract molecule names and
     pqr information (atom positions and partial charges)
   Reference document:
   Tripos Mol2 File Format, SYBYL 7.1 (Mid-2005),
   http://www.tripos.com/data/support/mol2.pdf *)

open Batteries

module A    = Array
module At   = Atom
module Dacc = Lacc.Dacc
module L    = List
module F    = Filename
module Mol  = Molecule
module MU   = My_utils

open Printf

let mol_start_str   = "@<TRIPOS>MOLECULE"
let mol_start_tag   = Str.regexp_string mol_start_str
let atoms_start_tag = Str.regexp_string "@<TRIPOS>ATOM"
let atoms_end_tag   = Str.regexp_string "@<TRIPOS>BOND"
let bonds_end_tag   = Str.regexp_string "@<TRIPOS>SUBSTRUCTURE"

let starts_with prfx_regexp str =
  Str.string_match prfx_regexp str 0

exception Break

(* throw End_of_file once there is no more to read *)
let read_one_molecule counter input =
  (* skip lines until we start to read a new molecule *)
  while not (starts_with mol_start_tag (Legacy.input_line input))
  do () done;
  (* the line after the mol_start_tag is the molecule name in MOL2 format *)
  let molecule_name = Legacy.input_line input in
  (* skip lines until we start to read atoms *)
  while not (starts_with atoms_start_tag (Legacy.input_line input))
  do () done;
  (* read all atoms *)
  let atoms, _exn =
    MU.unfold_exc
      (fun () -> let line = Legacy.input_line input in
                 if starts_with atoms_end_tag line
                 then raise Break
                 else At.of_mol2_line line
      )
  in
  let mol = Mol.create molecule_name !counter atoms in
  incr counter;
  mol

let read_one_molecule_with_bonds counter input =
  let m = read_one_molecule counter input in
  (* read all bonds *)
  let bonds, _exn =
    MU.unfold_exc
      (fun () -> let line = Legacy.input_line input in
                 if starts_with bonds_end_tag line
                 then raise Break
                 else At.bond_of_mol2_line line
      )
  in
  (m, bonds)

(* read all molecules from fn *)
let read_molecules fn =
  let read_one =
    if F.check_suffix fn ".bin" then
      (* read from a previously parsed file *)
      Mol.from_chan
    else
      (* read from a regular MOL2 file *)
      let _ =
        if not (F.check_suffix fn ".mol2")
        then Log.warn "%s not a .mol2 file" fn;
      in
      read_one_molecule
  in
  MU.with_in_file fn (fun input ->
    let nb_molecules = ref 0 in
    let molecules, exn =
      MU.unfold_exc
        (fun () -> read_one nb_molecules input)
    in
    assert (exn = End_of_file);
    Log.info "%d molecule(s) in %s" !nb_molecules fn;
    molecules
  )

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
         (fun l -> starts_with mol_start_tag l)
         all_lines)
  in
  let nb_molecules = L.length all_molecules in
  (* put back the header of each molecule, removed by BatList.nsplit *)
  let res = A.of_list (L.map (fun m -> mol_start_str :: m) all_molecules) in
  Log.info "%d molecule(s) in %s" nb_molecules fn;
  res
