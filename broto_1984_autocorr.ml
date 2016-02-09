
(* trial at reproducing method from  *)
(* "Molecular structures: perception, autocorrelation descriptor and sar studies" *)
(* Pierre Broto, Gilles Moreau and Corinne Vandycke. *)
(* Eur. J. Med. Chem. 1984-19 N1 pp 66-70 *)

(* we can reproduce easily their contribution to LogP and steric *)
(* autocorrelation vectors *)
(* they use histograms without normalization with discretization step = 0.2A *)

open Printf
open Batteries

module Mol2 = Mol2_parser

let auto_correlation (atoms: Atom.atom list): (float * float) list =
  (* sort values by increasing distance: this is assumed later on *)
  let sort =
    List.fast_sort (fun (d1, _) (d2, _) -> BatFloat.compare d1 d2)
  in
  let autocorr l =
    let rec loop acc = function
      | [] -> acc
      | a1 :: others ->
        let res = List.map (Autocorr.process_atoms a1) others in
        loop (List.rev_append res acc) others
    in
    loop [] l
  in
  sort (autocorr atoms)

(* take while 'p' is true in 'l' but also return the rest *)
let take_while p l =
  let rec loop acc to_process =
    match to_process with
    | [] -> (List.rev acc, [])
    | x :: xs ->
      if p x then loop (x :: acc) xs
      else (List.rev acc, to_process)
  in
  loop [] l

let fsum l =
  let rec loop acc = function
    | [] -> acc
    | (_d, v) :: xs -> loop (acc +. v) xs
  in
  loop 0.0 l

let histo dx pairs =
  let rec loop curr_dx acc to_process =
    match to_process with
    | [] -> List.rev acc
    | _ ->
      let values, rest = take_while (fun (d, _) -> d <= curr_dx) to_process in
      let sum = fsum values in
      loop (curr_dx +. dx) ((curr_dx, sum) :: acc) rest
  in
  loop dx [] pairs

let b84_histo (atoms: Atom.atom list): (float * float) list =
  histo 0.2 (auto_correlation atoms)

let print_histo (hist: (float * float) list): unit =
  List.iter (fun (d, x) -> printf "%f %f\n" d x) hist;
  printf "\n"

let main () =
  Log.set_log_level Log.INFO;
  Log.color_on();
  (* default option values *)
  let query_file = ref "" in
  let db_file    = ref "" in
  let cmd_line = Arg.align [
      "-q" , Arg.Set_string query_file, "query.mol2 query molecule";
      "-db", Arg.Set_string db_file   , "db.mol2 database";
    ]
  in
  Arg.parse cmd_line ignore
    (sprintf "Example: %s -q query.mol2 -db database.mol2\n" Sys.argv.(0));
  if !query_file = "" && !db_file = "" then begin
    Log.fatal "-q and -db are mandatory";
    exit 1
  end;
  let db_molecules = Mol2.read_molecules !db_file in
  List.iter (fun (_name, _i, atoms) ->
      let histo = b84_histo atoms in
      print_histo histo
    ) db_molecules

let () = main ()
