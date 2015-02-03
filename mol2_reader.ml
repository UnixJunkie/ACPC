
(* test reading of mol2 files with partial charges added to the atoms *)

open Batteries
open Printf
open Molecule

module A    = Array
module Feat = Feature
module HT   = Hashtbl
module L    = List
module MU   = My_utils

let rec get_molecule_names acc = function
  | x :: ((y :: ys) as rest) ->
    if Str.string_match Mol2_parser.mol_start_tag x 0 then
      get_molecule_names (y :: acc) ys
    else
      get_molecule_names acc rest
  | _ ->
    Log.info "%d molecule(s) found" (L.length acc);
    L.rev acc

let get_molecules_with_names lines =
  let molecules_without_start_tags =
    L.filter
      ((<>) []) (* rm border effect of nsplit *)
      (L.nsplit
         (fun l -> Str.string_match Mol2_parser.mol_start_tag l 0)
         lines)
  in
  Log.info "%d molecule(s) found" (L.length molecules_without_start_tags);
  L.map
    (fun l ->
      (* extract the name *)
      MU.hd l ~err:"Mol2_reader: get_molecules_with_names: molecule name",
      (* put back the start tag *)
      "@<TRIPOS>MOLECULE" :: l)
    molecules_without_start_tags

let main () =

  Log.set_log_level Log.INFO;
  Log.color_on();

  let argc = A.length Sys.argv in
  if argc <> 2 && argc <> 3 then (
    Log.fatal "\nusage: %s input.mol2 [-uniq|-stddev]" Sys.argv.(0);
    exit 1
  );

  let input = Sys.argv.(1) in
  let args = Array.to_list Sys.argv in
  let filter = L.mem "-uniq"   args in
  let stddev = L.mem "-stddev" args in
  if stddev then
    let mols = Mol2_parser.read_molecules input in
    L.iter
      (fun mol ->
         let charges = L.map (fun a -> Atom.charge Feat.Charge a) mol.atoms in
         let pc_mean, pc_stddev = MU.mean_stddev charges in
         printf "%s %f %f\n" mol.name pc_mean pc_stddev;
      )
      mols
  else if filter then
    (* keep only the first molecule encountered with a given SMILE string *)
    let smiles = Hashtbl.create 10000 in
    let lines = MU.string_list_of_file input in
    let named_molecules = get_molecules_with_names lines in
    L.iter
      (fun (mol_name, mol_lines) ->
         let smile = snd (String.split mol_name ~by:" ") in
         if not (Hashtbl.mem smiles smile) then
           let _ = Hashtbl.add smiles smile () in
           L.iter (printf "%s\n") mol_lines
         else
           Log.warn "skipped %s" mol_name
      )
      named_molecules;
    Log.info "%d molecule(s) written out" (HT.length smiles)
  else
    let mols = Mol2_parser.read_molecules input in
    L.iter
      (fun mol ->
        printf "%s %d\n" mol.name mol.index;
        L.iter
          (fun a -> printf "%s\n" (Atom.to_mol2_string a))
          mol.atoms
      )
      mols
;;

main()
