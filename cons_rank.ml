
(* tool to make a consensus out of several .scores files *)

open Batteries
open Printf

module Cons = Consensus
module S    = String
module L    = List
module Log  = Log
module MU   = My_utils

let main () =

  Log.set_log_level Log.INFO;
  Log.color_on();

  (* default options *)
  let rank_files    = ref ""  in
  let rank_files_fn = ref ""  in
  let ratio         = ref 0.0 in
  let verbose       = ref false in
  let cmd_line = Arg.align [
    "-if", Arg.Set_string rank_files_fn,
    "scores_files a file containing all .scores files";
    "-rf", Arg.Set_string rank_files,
    "scores_file1,scores_file2,...,scores_fileN a comma separated list \
     of .scores files";
    "-p", Arg.Set_float ratio, "x a float such that 0.0 < x <= 1.0 ";
    "-v", Arg.Set verbose, " output intermediate results to stdout"
  ] in
  Arg.parse cmd_line ignore
    (sprintf "Example: %s -if rank_files\n" Sys.argv.(0));
  if !rank_files_fn = "" && !rank_files = "" then
    failwith "-if or -rf is mandatory"
  ;
  if !rank_files_fn <> "" && !rank_files <> "" then
    failwith "use -if XOR -rf"
  ;
  if !ratio <= 0.0 || !ratio > 1.0 then (
    Log.fatalf "option -p x: x is out of range";
    exit 1
  );
  let files =
    if !rank_files_fn <> "" then
      let lines = MU.string_list_of_file !rank_files_fn in
      match lines with
        | []  -> failwith ("no rank files found in: " ^ !rank_files_fn)
        | [_] -> failwith ("only one rank file in: " ^ !rank_files_fn)
        | _ -> lines
    else
      let lines = S.nsplit !rank_files ~by:"," in
      match lines with
        | []  -> failwith ("no rank files found in: " ^ !rank_files)
        | [_] -> failwith ("only one rank file in: " ^ !rank_files)
        | _ -> lines
  in
  Cons.do_consensus files "Min" !verbose !ratio
;;

main()
