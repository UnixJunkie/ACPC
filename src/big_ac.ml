
open Batteries
open Printf

module A    = Array
module AC   = Autocorr
module S    = String
module F    = Filename
module HT   = Hashtbl
module L    = List
module Log  = Dolog.Log
module MU   = My_utils
module Mol2 = Mol2_parser
module V3   = Vector3

let do_query query_file db_file dx =
  let nb_molecules = ref 0 in
  let res = ref [] in
  let query = match Mol2.read_molecules query_file with
    | [m] -> m
    | []  -> failwith ("no molecule in query file: " ^ query_file)
    | _ -> failwith ("more than one molecule in query file: " ^ query_file)
  in
  let query_basename = F.chop_extension query_file in
  let scores_fn = query_basename ^ ".scores" in
  let query_desc = Kde.linbin_autocorr None dx query in
  let db_file = open_in db_file in
  (try while true do
      let name, atoms = Mol2.read_one_molecule db_file in
      let candidate_ac =
        Kde.linbin_autocorr None dx (name, !nb_molecules, atoms) in
      let score = AC.correlate_linbin_autocorrs' query_desc candidate_ac in
      res := (name, score) :: !res;
      incr nb_molecules
  done with End_of_file -> ());
  close_in db_file;
  let grouped = L.group (fun (m1, _) (m2, _) -> S.compare m1 m2) !res in
  let scores_file = open_out scores_fn in
  L.iter
    (fun group -> let scores = L.map snd group in
                  let best = L.max scores in
                  let name =
                    fst (MU.hd group ~err:"big_ac.ml: do_query: name") in
                  fprintf scores_file "%s %f\n" name best
    )
    grouped;
  close_out scores_file;
  !nb_molecules

let main () =
  let start = Unix.gettimeofday() in
  Log.set_log_level Log.INFO;
  Log.color_on();
  Log.info "\n\nCopyright (C) 2014, Zhang Initiative Research Unit,\n\
            Institute Laboratories, RIKEN\n\
            2-1 Hirosawa, Wako, Saitama 351-0198, Japan\n";
  (* default option values *)
  let query_file  = ref ""    in
  let db_file     = ref ""    in
  let dx          = ref 0.005 in
  let cmd_line = Arg.align [
    "-q"   , Arg.Set_string query_file, "query.mol2 query";
    "-db"  , Arg.Set_string db_file   , "db.mol2 database";
    "-dx"  , Arg.Set_float dx         ,
    (sprintf "float X axis discretization (default: %f)" !dx)
  ] in
  Arg.parse cmd_line ignore
    (sprintf "Example: %s -q query.mol2 -db database.mol2\n" Sys.argv.(0));
  if !query_file = "" then begin
    Log.fatal "-q query_file.mol2 is mandatory";
    exit 1
  end;
  let nb_molecules_in_db = do_query !query_file !db_file !dx in
  let stop = Unix.gettimeofday() in
  let elapsed = stop -. start in
  Log.info "speed: %.2f molecules/s"
    ((float_of_int nb_molecules_in_db) /. elapsed)
;;

main()
