
(* enrichment rate calculation tool *)

open Batteries
open Legacy.Printf

module L  = List
module MU = My_utils
module S  = String

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
  MU.enforce_any_file_extension !input_fn ["scores"];
  (* load all the score labels *)
  let nb_molecules = ref 0 in
  let score_labels, exn =
    MU.with_in_file !input_fn (fun input ->
      MU.unfold_exc
        (fun () ->
           let l = Legacy.input_line input in
           try Legacy.Scanf.sscanf l "%s %f %d"
                 (fun name score index ->
                   let label = S.starts_with name "active" in 
                   incr nb_molecules;
                   (name, score, index, label)
                 )
           with _ ->
             failwith ("er_tool.ml: could not parse: " ^ l)
        )
    )
  in
  assert(exn = End_of_file);
  let top_n, top_actives_rate, rand_actives_rate, enr_rate =
    ROC.enr_rate !ratio score_labels
  in
  Log.info "%s p: %.2f n: %d topAR: %.2f randAR: %.2f ER: %.2f"
    !input_fn !ratio top_n top_actives_rate rand_actives_rate enr_rate
;;

main()
