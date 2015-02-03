
open Batteries
open Molecule
open Printf

module A    = Array
module At   = Atom
module Fn   = Filename
module Feat = Feature
module Itv  = Interval_tree.Interval
module Itvt = Interval_tree
module L    = List
module MU   = My_utils
module Mol  = Molecule
module Mol2 = Mol2_parser
module Pl   = Pl_parser
module Pqr  = Pqr_parser
module S    = String
module V3   = Vector3

let read_query_molecule feature query_fn =
  let read_molecules = match feature with
    | Feat.Charge -> Mol2.read_molecules
    | Feat.Radius -> Pqr.read_molecules
    | Feat.LogP   -> Pl.read_molecules
  in
  let query = read_molecules query_fn in
  match query with
    | []  -> failwith ("read_query_molecule: no molecule in query file: "
                       ^ query_fn)
    | [q] -> q
    | _ -> failwith ("read_query_molecule: more than one molecule in: "
                     ^ query_fn)

(* get the reading function to read a database of molecules *)
let get_reader db_fn =
  if Fn.check_suffix db_fn ".bin" then
    (* read from a binary file *)
    Mol.from_chan
  else
    let feature = Feat.of_filename db_fn in
    match feature with
      | Feat.Charge -> (* a regular MOL2 file *)
        Mol2.read_one_molecule
      | Feat.Radius -> (* a regular PQR file *)
        Pqr.read_one_molecule
      | Feat.LogP -> (* a regular PL file *)
        Pl.read_one_molecule

(* try to count molecules
   returns: Some count | None *)
let try_counting_molecules db_fn =
  if Fn.check_suffix db_fn ".bin" then
    None
  else
    let feature = Feat.of_filename db_fn in
    let cmd_str =
      match feature with
      | Feat.Charge -> (* a MOL2 file *)
        "egrep -c '^@<TRIPOS>MOLECULE$' " ^ db_fn
      | Feat.Radius -> (* a PQR file *)
        "egrep -c '^COMPND ' " ^ db_fn
      | Feat.LogP -> (* a PL file *)
        "egrep -c '^MOL ' " ^ db_fn
    in
    Log.info "running:\n%s" cmd_str;
    let cmd_out = S.strip (MU.get_command_output cmd_str) in
    Some (Int.of_string cmd_out)

(* delete each atom in turn from the input molecule in 'fn' and
   return the list of all such molecules *)
let scan_atoms feature fn =
  let mol = read_query_molecule feature fn in
  let scanned = L.mapi (fun i _ -> MU.remove_at i mol.atoms) mol.atoms in
  L.mapi (fun i atoms -> Mol.create mol.name i atoms) scanned

let output_optimized_query feature out_fn optimized_query_molecule =
  let write_out_fun = match feature with
    | Feat.Charge -> Mol.pseudo_mol2_file
    | Feat.LogP   -> Mol.pl_file
    | Feat.Radius -> Mol.pseudo_pqr_file
  in
  MU.string_list_to_file out_fn (write_out_fun optimized_query_molecule)

let atom_scan feature initial_query_auc (query_fn, query_mol) scan_out_fn query_fun =
  let to_scan = scan_atoms feature query_fn in
  let bild_fn = query_fn ^ ".bld" in
  let aucs = L.map query_fun to_scan in
  (* do_query feature cmp_f false (query_fn, query_mol)
     !dx db_kdes !post_filter !db_file false true false 0 "" *)
  let deltas = L.map (fun auc -> auc -. initial_query_auc) aucs in
  Log.info "to look at the result type:\n\
            chimera %s %s" query_fn bild_fn;
  let atom_and_deltas = L.combine query_mol.atoms deltas in
  let interesting_atoms =
    L.map
      fst
      (* only keep atoms which gave a better AUC
         when they were taken into account *)
      (L.filter (fun (_atom, dAUC) -> dAUC < 0.0) atom_and_deltas)
  in
  let color_str, radius, transparency =
    Feat.chimera_sphere_properties feature in
  MU.string_list_to_file bild_fn
    (((sprintf ".transparency %f" transparency) ::
         (".color " ^ color_str) ::
         (* atoms that passed the filter *)
         (L.map
            (fun a ->
              sprintf ".sphere %s %f" (V3.to_string (At.position a)) radius)
            interesting_atoms))
     @ (* all atoms with their deltaAUC as comments *)
       (L.map
          (fun (a, dAUC) ->
            sprintf ".comment a: %s dAUC: %f" (At.to_string feature a) dAUC)
          atom_and_deltas)
    );
  output_optimized_query
    feature scan_out_fn (Mol.create query_mol.name 0 interesting_atoms);
  L.length aucs

let auto_correlation on atoms =
  (* CACPC doesn't sort values by increasing distance because
     we don't need it anymore (ACPC was enforcing this) *)
  let neg_acs = ref [] in
  let pos_acs = ref [] in
  let rec loop = function
    | [] -> ()
    | a1 :: others ->
      let _ =
        L.iter
          (fun a2 ->
            let (_d_ij, q_ij) as ac_ij = At.process_atoms on a1 a2 in
            if q_ij < 0.0 then
              neg_acs := ac_ij :: !neg_acs
            else (* q_ij >= 0.0 *)
              pos_acs := ac_ij :: !pos_acs
          )
          others
      in
      loop others
  in
  loop atoms;
  (!neg_acs, !pos_acs)

(* store auto correlation values into an interval tree, for triangular kernel
   of half base 1/a *)
let index_ac a ac_values =
  let dx = 1.0 /. a in
  let triplets =
    L.map
      (fun ((d_ij, _q_ij) as ac) -> Itv.create (d_ij -. dx) (d_ij +. dx) ac)
      ac_values
  in
  Itvt.create triplets

(* triangular kernel:
   1 at x = 0
   linearly decreasing to 0 until abs(x) = 1/a
   0 everywhere else *)
let triangle a x =
  let xabs = abs_float x in
  if xabs > 1.0 /. a then
    0.0
  else
    1.0 -. a *. xabs

(* Inspired by Eq. 1 from
   J. Chem. Inf. Comput. Sci. 1996, 36, 118-127
   "Chemical Similarity Using Physiochemical Property Descriptors"
   Simon K. Kearsley, Susan Sallamack, Eugene M. Fluder, Joseph D. Andose,
   Ralph T. Mosley, and Robert P. Sheridan *)
let kernel_overlap a ac1_idx ac2 =
  L.fold_left
    (fun acc (dj, qj) ->
       let require_kernel_eval = Itvt.query ac1_idx dj in
       L.fold_left
         (fun acc2 itv ->
            let di, qi = Itv.(itv.value) in (* query molecule contribution *)
            let w_ij = triangle a (di -. dj) in (* current contrib. weight *)
            acc2 +. (w_ij *. qi *. qj)
         )
         acc
         require_kernel_eval
    )
    0.0
    ac2

(* overlap measure based on a given kernel *)
let continuous_overlap on a (query_neg_ac_idx, query_pos_ac_idx) candidate =
  let neg_ac, pos_ac = auto_correlation on candidate.atoms in
  kernel_overlap a query_neg_ac_idx neg_ac +.
  kernel_overlap a query_pos_ac_idx pos_ac

(* cross correlation for linearly binned auto correlation vectors *)
let cross_correlation a1 a2 =
  let n = min (A.length a1) (A.length a2) in
  let res = ref 0. in
  for i = 0 to n - 1 do
    res := !res +. ((A.unsafe_get a1 i) *. (A.unsafe_get a2 i))
  done;
  !res

let correlate_linbin_autocorrs
    (_i1, _name1, (n1, p1))
    (_i2, _name2, (n2, p2)) =
  0.5 *. ((cross_correlation n1 n2) +. (cross_correlation p1 p2))

let correlate_linbin_autocorrs'
    (_name1, (n1, p1))
    (_name2, (n2, p2)) =
  0.5 *. ((cross_correlation n1 n2) +. (cross_correlation p1 p2))

let score f a1 a2 =
  let n = min (A.length a1) (A.length a2) in
  let l1 = A.to_list (A.sub a1 0 n) in
  let l2 = A.to_list (A.sub a2 0 n) in
  f l1 l2

let score_linbin_autocorrs f
    ((_i1 : int), (_name1 : string), (n1, p1))
    ((_i2 : int), (_name2 : string), (n2, p2)) =
  0.5 *. ((f n1 n2) +. (f p1 p2))

let tanimoto_linbin_autocorrs    = score_linbin_autocorrs (score MU.tanimoto_coeff)
let tversky_ref_linbin_autocorrs = score_linbin_autocorrs (score MU.tversky_ref)
let tversky_db_linbin_autocorrs  = score_linbin_autocorrs (score MU.tversky_db)
