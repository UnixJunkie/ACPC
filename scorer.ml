
(* tool to combine .scores files using a weight for each feature space.
   Experimental. Scores should be trimmed probably: ACPC(query, query) should be the
   max score allowed against a given query.
*)

open Batteries
open Legacy.Printf

module A   = Array
module HT  = Hashtbl
module Sco = Scores
module S   = String
module L   = List
module MU  = My_utils

(* remove the mean then divide by the stddev each score *)
let center_then_scale name_score_index_labels =
  let scores_only =
    L.map
      (fun (_name, score, _index, _label) -> score)
      name_score_index_labels
  in
  let mean, stddev = MU.mean_stddev scores_only in
  L.map
    (fun (name, score, index, label) ->
       (name, (score -. mean) /. stddev, index, label)
    )
    name_score_index_labels

(* load a .scores file into a list of (name, score, index, label)
   for use by ROC.auc / ROC.enr_rate.
   SCORES ARE CENTERED THEN SCALED
*)
let read_scores_file fn =
  let lines = MU.string_list_of_file fn in
  let score_labels =
    L.map
      (fun line ->
        try
          Scanf.sscanf line "%s %f %d"
            (fun name score index ->
              let label = S.starts_with name "active" in
              (name, score, index, label)
           )
        with exn ->
          Log.fatal "scorer.ml: read_scores_file %s: couldn't parse %s"
            fn line;
          raise exn
      )
      lines
  in
  center_then_scale score_labels

let rescore mol2_scores pqr_scores pl_scores w_mol2 w_pqr w_pl out_fn =
  let res = HT.create 10000 in
  (* retain the mol2 score of each name-index-label triplet *)
  L.iter
    (fun (name, score, index, label) ->
       let r = Sco.create () in
       HT.add res name (Sco.set_mol2_score r score, index, label)
    )
    mol2_scores
  ;
  (* retain the pqr score of each name-index-label triplet *)
  L.iter
    (fun (name, score, _index, _label) ->
       let r, index, label = HT.find res name in
       HT.replace res name (Sco.set_pqr_score r score, index, label)
    )
    pqr_scores
  ;
  (* retain the pl score of each name-index-label triplet *)
  L.iter
    (fun (name, score, _index, _label) ->
       let r, index, label = HT.find res name in
       HT.replace res name (Sco.set_pl_score r score, index, label)
    )
    pl_scores
  ;
  (* rescore *)
  let new_scores = ref [] in
  MU.with_out_file out_fn (fun out ->
    fprintf out "#mol2, pqr and pl scores were centered then scaled\n";
    fprintf out "#w_mol2: %f w_pqr: %f w_pl: %f\n" w_mol2 w_pqr w_pl;
    fprintf out "#mol_name mol2_score pqr_score pl_score combined_score\n";
    HT.iter
      (fun name (scores, index, label) ->
        let new_score = Sco.rescore scores w_mol2 w_pqr w_pl in
        let mol2_score, pqr_score, pl_score = Sco.to_triplet scores in
        fprintf out "%s %f %f %f %f\n"
          name mol2_score pqr_score pl_score new_score;
        new_scores := (name, new_score, index, label) :: !new_scores
      )
      res;
  );
  !new_scores

let main () =

  Log.set_log_level Log.INFO;
  Log.color_on();

  let argc = A.length Sys.argv in

  if argc <> 8 then (
    (*        0  1           2          3         4      5     6    7 *)
    Log.fatal
      "usage: %s mol2_scores pqr_scores pl_scores w_mol2 w_pqr w_pl scores_out_file"
      Sys.argv.(0)
    ;
    exit 1
  );

  let mol2_scores_fn = Sys.argv.(1) in
  let pqr_scores_fn  = Sys.argv.(2) in
  let pl_scores_fn   = Sys.argv.(3) in
  let w_mol2         = Float.of_string Sys.argv.(4) in
  let w_pqr          = Float.of_string Sys.argv.(5) in
  let w_pl           = Float.of_string Sys.argv.(6) in
  let out_fn         = Sys.argv.(7) in
  assert(w_mol2 >= 0.0 && w_mol2 <= 1.0);
  assert(w_pqr  >= 0.0 && w_pqr  <= 1.0);
  assert(w_pl   >= 0.0 && w_pl   <= 1.0);
  let mol2_scores = read_scores_file mol2_scores_fn in
  let pqr_scores  = read_scores_file pqr_scores_fn in
  let pl_scores   = read_scores_file pl_scores_fn in
  let score_labels =
    rescore mol2_scores pqr_scores pl_scores w_mol2 w_pqr w_pl out_fn in
  let auc = ROC.auc score_labels in
  let er = ROC.only_ER 0.01 score_labels in
  printf "%f %f %f %.2f %.2f\n" w_mol2 w_pqr w_pl auc er
;;

main()
