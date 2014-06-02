
(* output protein atoms that are not too far away from any ligand atom *)

open Batteries

module A   = Array
module AA  = Amino_acid
module At  = Atom
module FN  = Filename
module F   = File
module L   = List
module Pdb = Pdb_parser
module S   = String

(* like Python's readlines() *)
let string_list_of_file f =
  L.of_enum (F.lines_of f)

let atof = float_of_string

let enforce_file_extension filename ext =
  if not (FN.check_suffix filename ext) then
    failwith ("file " ^ filename ^ " not a " ^ ext)
  else
    ()

let xyz_of_pqr_line l =
  let atom = At.of_pqr_line l in
  (At.to_xyz atom, l)

let squared_dist (x1, y1, z1) (x2, y2, z2) =
  (x1 -. x2) ** 2. +.
  (y1 -. y2) ** 2. +.
  (z1 -. z2) ** 2.

let is_atom_or_heteroatom l =
  S.starts_with l "ATOM" ||
  S.starts_with l "HETATM"

let fail_on_empty_list l f =
  match l with
      [] -> Printf.fprintf stderr "no lines selected from %s\n" f;
            exit 1;
    | _  -> ()
;;

type 'a line_type = Water | Something of 'a | Unknown

let main () =
  let argc = A.length Sys.argv in
  if argc != 4 then begin
    Printf.fprintf stderr "%s"
      (" ERROR: incorrect number of parameters\n" ^
          " USAGE: " ^ Sys.argv.(0) ^ " protein_only_pqr ligand_only_pqr radius\n" ^
          (*           0                1                2               3      4 *)
          " warning: ALL the ATOM and HETATM lines of your files will be\n" ^
          "          considered. Edit files before if you are not happy\n" ^
          "          with this.\n");
    exit 1;
  end else begin
    let protein_file, ligand_file = Sys.argv.(1), Sys.argv.(2) in
    enforce_file_extension protein_file ".pqr";
    enforce_file_extension ligand_file  ".pqr";
    let acceptance_radius = (atof (Sys.argv.(3))) ** 2. in
    let protein_atoms     =
      L.map xyz_of_pqr_line
        (L.filter is_atom_or_heteroatom (string_list_of_file protein_file)) in
    fail_on_empty_list protein_atoms protein_file;
    let ligand_atoms      =
      L.map xyz_of_pqr_line
        (L.filter is_atom_or_heteroatom (string_list_of_file ligand_file))  in
    fail_on_empty_list ligand_atoms ligand_file;
    let selected_protein_lines =
      L.filter
        (fun (xyz_p, _) ->
          let min_dist = L.min
            (L.map (fun (xyz_l, _) -> squared_dist xyz_p xyz_l) ligand_atoms) in
          min_dist < acceptance_radius)
        protein_atoms in
    let classes =
      L.map
        (fun (_, l) ->
          try
            let res = AA.aa_three_of_string (Pdb.res_name_of_pdb_line l) in
            let atom = At.atom_type_of_string (Pdb.atom_name_of_pdb_line l) in
            Printf.printf "%s\n" l; (* print out selected lines *)
            (* Some (AA.classify res) (\* per amino acid classification *\) *)
            (* per amino_acid and atom classification *)
            Something (At.classify_atom (res, atom))
          with AA.Unknown_aa_three exotic_aa ->
            if exotic_aa = "WAT" then Water
            else
              let _ = Printf.eprintf "Exotic AA: %s\n" exotic_aa in
              Unknown
        )
        selected_protein_lines
    in
    let water_count       = ref 0 in
    let unknown_count     = ref 0 in
    let hydrophobic_count = ref 0 in
    let aromatic_count    = ref 0 in
    let cationic_count    = ref 0 in
    let anionic_count     = ref 0 in
    let hb_don_count      = ref 0 in
    let hb_acc_count      = ref 0 in
    L.iter At.(function
      | Water   -> incr water_count
      | Unknown -> incr unknown_count
      | Something (props, _r) ->
        L.iter
          (function
            | Hydrophobic -> incr hydrophobic_count
            | Aromatic    -> incr aromatic_count
            | Cationic    -> incr cationic_count
            | Anionic     -> incr anionic_count
            | HbD         -> incr hb_don_count
            | HbA         -> incr hb_acc_count
          )
          props
    )
    classes;
    let total =
      float_of_int ((L.length classes) - !water_count) in
    Printf.eprintf "#Hydrophobic Aromatic Cationic Anionic HbD HbA Unknown\n";
    Printf.eprintf "%.2f "  (float_of_int !hydrophobic_count /. total);
    Printf.eprintf "%.2f "  (float_of_int !aromatic_count    /. total);
    Printf.eprintf "%.2f "  (float_of_int !cationic_count    /. total);
    Printf.eprintf "%.2f "  (float_of_int !anionic_count     /. total);
    Printf.eprintf "%.2f "  (float_of_int !hb_don_count      /. total);
    Printf.eprintf "%.2f "  (float_of_int !hb_acc_count      /. total);
    Printf.eprintf "%.2f\n" (float_of_int !unknown_count     /. total);
  end
;;

main()
