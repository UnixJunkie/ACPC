
(* functions for ROC analysis *)

open Batteries
open Printf

module F  = File
module HT = My_utils.Hashtbl
module L  = List
module MU = My_utils
module S  = String

(* group molecules by name and keep only the best score for each name
   THIS IS NEEDED BY MULTIPLE COMFORMERS DATASETS *)
let unique dup_scores =
  let map = HT.create 10_000 in
  let removed = ref 0 in
  (* keep only best score of all molecules with a same name *)
  L.iter
    (fun ((mol_name, curr_score, _index, _label) as score_label) ->
       (* preserve order by indices, to ease software validation *)
       try
         let _name, prev_score, _j, _label = HT.find map mol_name in
         incr removed;
         if curr_score > prev_score then
           HT.replace map mol_name score_label
         ;
       with Not_found ->
         HT.add map mol_name score_label
    )
    dup_scores;
  if !removed > 0 then
  Log.warnf "removed %d molecules (duplicated names)" !removed;
  (* return them in the same order that they were found *)
  let values = HT.values map in
  L.fast_sort
    (fun (_n1, _s1, i1, _l1) (_n2, _s2, i2, _l2) -> Int.compare i1 i2)
    values

let dump_scored_labels fn scores =
  Log.info (lazy (sprintf "writing scores and labels in %s" fn));
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
  Log.info (lazy (sprintf "writing names, scores and ranks in %s" fn));
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
            then Log.infof "HTQ: %s %f %d" mol_name score i;
            fprintf out "%s %f %d\n" mol_name score i)
          highest_scores_first
  );
  L.take top_n highest_scores_first

let read_AUC_from_string s =
  try
    Scanf.sscanf s "Area Under Curve = %f" identity
  with Scanf.Scan_failure _ ->
    Log.error (lazy (sprintf "ROC.ml: read_AUC_from_string: \
                              no AUC could be read"));
    0.0

let trapezoid_surface x1 x2 y1 y2 =
  let base = abs_float (x1 -. x2) in
  let height = 0.5 *. (y1 +. y2) in
  base *. height

(* area under the ROC curve given an unsorted list of score-labels
   TP cases have the label set to true
   FP cases have the label unset *)
let auc score_labels =
  let high_scores_first =
    L.sort
      (fun (_n1, s1, _i1, _l1) (_n2, s2, _i2, _l2) -> Float.compare s2 s1)
      score_labels
  in
  let fp, tp, fp_prev, tp_prev, a, _p_prev =
    L.fold_left
      (fun (fp, tp, fp_prev, tp_prev, a, p_prev) (_ni, si, _ii, li) ->
        let new_a, new_p_prev, new_fp_prev, new_tp_prev =
          if si <> p_prev then
            a +. trapezoid_surface fp fp_prev tp tp_prev,
            si,
            fp,
            tp
          else
            a,
            p_prev,
            fp_prev,
            tp_prev
        in
        let new_tp, new_fp =
          if li then
            tp +. 1., fp
          else
            tp, fp +. 1.
        in
        (new_fp, new_tp, new_fp_prev, new_tp_prev, new_a, new_p_prev)
      )
      (0., 0., 0., 0., 0., neg_infinity)
      high_scores_first
  in
  (a +. trapezoid_surface fp fp_prev tp tp_prev) /. (fp *. tp)

(* proportion of actives given an unsorted list of score-labels
   TP cases have the label set to true
   FP cases have the label unset
   returns: (nb_molecules, actives_rate) *)
let actives_rate score_labels =
  let tp_count, fp_count =
    L.fold_left
      (fun (tp_c, fp_c) (_name, _score, _index, label) ->
         if label then
           (tp_c + 1, fp_c)
         else
           (tp_c, fp_c + 1)
      )
      (0, 0)
      score_labels
  in
  let itof = float_of_int in
  let nb_molecules = tp_count + fp_count in
  (nb_molecules,
   (itof tp_count) /. (itof nb_molecules))

(* enrichment rate at x (e.g. x = 0.01 --> ER @ 1%) given a list
   of unsorted score-labels
   returns: (top_n, top_actives_rate, rand_actives_rate, enr_rate) *)
let enr_rate p score_labels =
  let nb_molecules, rand_actives_rate = actives_rate score_labels in
  let top_n = Float.round_to_int (p *. (float_of_int nb_molecules)) in
  let top_p_percent_molecules =
    L.take top_n
      (L.sort
         (fun (_n1, s1, _i1, _l1) (_n2, s2, _i2, _l2) ->
            (* descending sort based on scores *)
            Float.compare s2 s1)
         score_labels)
  in
  let _, top_actives_rate = actives_rate top_p_percent_molecules in
  let enr_rate = top_actives_rate /. rand_actives_rate in
  (top_n, top_actives_rate, rand_actives_rate, enr_rate)

let only_ER p score_labels =
  let _, _, _, er = enr_rate p score_labels in
  er
