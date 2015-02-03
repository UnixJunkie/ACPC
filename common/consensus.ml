
open Batteries
open Printf

module L  = List
module MU = My_utils
module S  = String

(* What we will do with the ranks *)
type consensus_policy = Min

let policy_of_string = function
  | "Min" -> Min
  | s -> failwith ("cons_rank.ml: policy_of_string: unimplemented policy: " ^ s)

let parse_rank_line l =
  try Scanf.sscanf l "%s %f %d" (fun name score _index -> (name, score))
  with _ -> failwith ("cons_rank.ml: parse_rank_line: cannot parse: " ^ l)

(* compute ranks based on scores (highest score become lowest rank,
   i.e. top of the ranked list of compounds)
   then return names and ranks, sorted by molecule names *)
let name_ranks_of_file fn =
  let lines        = MU.string_list_of_file fn   in
  let parsed       = L.map parse_rank_line lines in
  let top_scores_first =
    L.fast_sort
      (fun (_n1, s1) (_n2, s2) -> Float.compare s2 s1)
      parsed
  in
  let ranked_by_score =
    L.mapi
      (fun i (n, _s) -> (n, [i]))
      top_scores_first
  in
  L.fast_sort
    (fun (n1, _r1) (n2, _r2) -> S.compare n1 n2)
    ranked_by_score

let do_consensus rank_files policy_str verbose ratio =
  let policy = policy_of_string policy_str in
  match rank_files with
    | [] -> 
      Log.fatal "Consensus.do_consensus: rank_files = []";
      exit 1
    | rf1 :: rfs ->
      let name_ranks_1 = name_ranks_of_file rf1 in
      let max_rank = float_of_int (L.length name_ranks_1) in
      let for_consensus =
        L.fold_left
          (fun acc rfi ->
            let name_ranks_i = name_ranks_of_file rfi in
            L.map2
              (fun (n1, ranks) (n2, r2) ->
                if not (S.equal n1 n2) then
                  failwith "sort by molecule names went bad";
                (* sort by names went OK and working on same molecules sets *)
                match r2 with
                  | [rank2] -> (n1, rank2 :: ranks)
                  | _ -> failwith
                    "only the acc can have more than one score per molecule"
              )
              acc name_ranks_i
          )
          name_ranks_1
          rfs
      in
      let reduce l =
        match policy with
          | Min -> float_of_int (L.reduce Int.min l)
      in
      let consensus =
        let unused_index = -1 in
        L.map
          (fun (n, ranks) -> (unused_index, n, reduce ranks))
          for_consensus
      in
      (if verbose then
          L.iter
            (fun (_i, n, r) -> Printf.printf "%s %f\n" n r)
            consensus
      );
      (* best molecule must have highest score:
         we convert ranks into score-labels here *)
      let score_labels =
        L.map
          (fun (index, name, cons_rank) ->
             let score = max_rank -. cons_rank in
             let label = S.starts_with name "active" in
             (name, score, index, label)
          )
          consensus
      in
      let auc = ROC.auc score_labels in
      let _top_n, _top_actives_rate, _rand_actives_rate, enr_rate =
        ROC.enr_rate ratio score_labels
      in
      printf "%.2f %.2f\n" auc enr_rate
