
(* AUC and Enrichment Rate calculation tool for .score-label files *)

open Batteries
open Legacy.Printf

module L     = List
module MU    = My_utils
module S     = String
module Scanf = Legacy.Scanf

let score_label_sscanf nb_molecules l =
  Scanf.sscanf l "%f %d"
    (fun score int_label ->
      let label = (int_label = 1) in 
      incr nb_molecules;
      ("", score, !nb_molecules, label)
    )

let scores_sscanf nb_molecules l =
  Scanf.sscanf l "%s@ %f %d"
    (fun name score _index ->
       let label = S.starts_with name "active" in 
       incr nb_molecules;
       ("", score, !nb_molecules, label)
    )

let main () =

  Log.set_log_level Log.INFO;
  Log.color_on();

  (* default options *)
  let input_fn = ref "" in
  let ratio = ref 0.0 in
  let cmd_line = Arg.align [
    "-i", Arg.Set_string input_fn, "in.scores input file";
    "-p", Arg.Set_float ratio, "x a float such that 0.0 < x <= 1.0 ";
  ] in
  Arg.parse cmd_line ignore
    (sprintf "Example:\n%s -i query.mol2.scores" Sys.argv.(0));
  (* check options *)
  if !input_fn = "" then (
    Log.fatal "-i is mandatory";
    exit 1
  );
  if (!ratio <= 0.0) || (!ratio > 1.0) then (
    Log.fatal "option -p x: x is out of range";
    exit 1
  );
  let sscanf =
    if S.ends_with !input_fn ".score-label" ||
      S.ends_with !input_fn ".scored-label"
    then
      score_label_sscanf
    else (
      if S.ends_with !input_fn ".scores" then
        scores_sscanf
      else (
        Log.fatal "auc_tool.ml: %s neither a .scores nor a .score-label"
          !input_fn;
        exit 1
      )
    )
  in
  (* load all the score labels *)
  let nb_molecules = ref 0 in
  let score_labels, exn =
    MU.with_in_file !input_fn (fun input ->
      MU.unfold_exc
        (fun () -> let l = Legacy.input_line input in
                   try sscanf nb_molecules l
                   with _ ->
                     failwith ("auc_tool.ml: could not parse: " ^ l)
        )
    )
  in
  if exn <> End_of_file then
    raise exn
  ;
  let auc = ROC.auc score_labels in
  let _top_n, _top_actives_rate, _rand_actives_rate, enr_rate =
    ROC.enr_rate !ratio score_labels
  in
  printf "%.2f %.2f\n" auc enr_rate
;;

main()
