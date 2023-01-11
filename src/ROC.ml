
(* functions for ROC analysis *)

open Batteries
open Printf

module F  = File
module L  = List
module HT = Hashtbl
module MU = My_utils
module S  = String

(* group molecules by name and keep only the best score for each name
   THIS IS NEEDED BY MULTIPLE COMFORMERS DATASETS *)
let unique dup_scores =
  let map = HT.create 10000 in
  let removed = ref 0 in
  (* keep only best score for all molecules with same name *)
  L.iter
    (fun (i, mol_name, curr_score) ->
       (* i is used to preserve the order from the input file to ease
          software validation *)
       try
         let _j, prev_score = HT.find map mol_name in
         incr removed;
         if curr_score > prev_score
         then HT.replace map mol_name (i, curr_score);
       with Not_found ->
         HT.add map mol_name (i, curr_score)
    )
    dup_scores;
  if !removed > 0 then
  Log.warn "removed %d molecules (duplicated names)" !removed;
  (* return them in the same order that they were found *)
  let key_values = MU.key_values map in
  let to_sort =
    L.map
      (fun (mol_name, (i, score)) -> (i, mol_name, score))
      key_values
  in
  L.fast_sort
    (fun (i1, _n1, _s1) (i2, _n2, _s2) -> Int.compare i1 i2)
    to_sort

let dump_scored_labels fn scores =
  Log.info "writing scores and labels in %s" fn;
  F.with_file_out fn (fun roc_out ->
    L.iter
      (fun (_i, mol_name, score) ->
        let is_active =
          if S.starts_with mol_name "active"
          then 1
          else 0
        (* this is the format the CROC Python module reads in *)
        in fprintf roc_out "%f %d\n" score is_active
      )
      scores
  )

(* dump ranks and return the top_n scoring molecules *)
let dump_ranks fn scores maybe_list_htq top_n =
  Log.info "writing names, scores and ranks in %s" fn;
  let highest_scores_first =
    L.fast_sort (fun (_i1, _n1, x) (_i2, _n2, y) -> Float.compare y x) scores
  in
  F.with_file_out fn (fun out ->
    match maybe_list_htq with
      | None ->
        L.iter
          (fun (i, mol_name, score) ->
            fprintf out "%s %f %d\n" mol_name score i)
          highest_scores_first
      | Some s ->
        L.iter
          (fun (i, mol_name, score) ->
            if score > s
            then Log.info "HTQ: %s %f %d" mol_name score i;
            fprintf out "%s %f %d\n" mol_name score i)
          highest_scores_first
  );
  L.take top_n highest_scores_first

let read_AUC_from_string s =
  try
    Scanf.sscanf s "Area Under Curve = %f" identity
  with Scanf.Scan_failure _ ->
    Log.error "ROC.ml: read_AUC_from_string: \
               no AUC could be read";
    0.0
