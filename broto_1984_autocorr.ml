
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
module MU = My_utils

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

let print_histo (hist: (float * float) list): unit =
  List.iter (fun (d, x) -> printf "%f %f\n" d x) hist;
  printf "\n"

let print_autocorr =
  print_histo

let b84_histo (atoms: Atom.atom list): (float * float) list =
  let ac = auto_correlation atoms in
  (* print_autocorr ac; *)
  histo 0.2 ac

let truncate_histo min_len curr_len histo =
  if curr_len > min_len then
    List.take min_len histo
  else
    histo

let main () =
  Log.set_log_level Log.INFO;
  Log.color_on();
  (* default option values *)
  let query_file = ref "" in
  let db_file    = ref "" in
  let scores_file = ref "" in
  let cmd_line = Arg.align
      ["-q" , Arg.Set_string query_file , "query.mol2 query molecule";
       "-db", Arg.Set_string db_file    , "db.mol2 database";
       "-o" , Arg.Set_string scores_file, "out.scores score-labels file"]
  in
  Arg.parse cmd_line ignore
    (sprintf "Example: %s -q query.mol2 -db database.mol2\n" Sys.argv.(0));
  if !query_file = "" || !db_file = "" || !scores_file = "" then begin
    Log.fatal "-q, -db and -o are all mandatory";
    exit 1
  end;
  let _, _, query_atoms = 
    match Mol2.read_molecules !query_file with
    | [one] -> one
    | _ -> failwith (sprintf "not one molecule in %s" !query_file)
  in
  let query_histo = b84_histo query_atoms |> List.map snd in (* discard distances *)
  let query_histo_len = List.length query_histo in
  (* print_histo query_histo; *)
  let db_molecules = Mol2.read_molecules !db_file in
  MU.with_out_file !scores_file (fun out ->
      List.iter (fun (name, _i, atoms) ->
          let histo = b84_histo atoms |> List.map snd in
          let histo_len = List.length histo in
          let min_len = min query_histo_len histo_len in
          let query_histo = truncate_histo min_len query_histo_len query_histo in
          let histo = truncate_histo min_len histo_len histo in
          let score = MU.tversky_ref query_histo histo in
          let is_active =
            if BatString.starts_with name "active"
            then 1
            else 0
          in
          (* this is the format the CROC Python module reads in *)
          fprintf out "%f %d\n" score is_active
        ) db_molecules
    )

let () = main ()
