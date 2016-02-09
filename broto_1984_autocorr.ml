
(* trial at reproducing methods from: *)

(* "Molecular structures: perception, autocorrelation descriptor and sar studies" *)
(* Pierre Broto, Gilles Moreau and Corinne Vandycke. *)
(* Eur. J. Med. Chem. 1984-19 N1 pp 66-70. *)
(* We can reproduce easily their contribution to LogP and steric *)
(* autocorrelation vectors. *)
(* Albeit they were using their own contribution to LogP model. *)
(* They use histograms (without normalization) and discretization step = 0.2A *)

(* there is also code to reproduce the autocorrelation vector encoding from *)
(* "Autocorrelation of Molecular Electrostatic Potential Surface Properties *)
(* Combined with Partial Least Squares Analysis as Alternative Attractive *)
(* Tool to Generate Ligand-Based 3D-QSArs" *)
(* Stefano Moro et. al. *)
(* Current Drug Discovery Technologies 2005. *)
(* Albeit they were using EP at the molecular surface while we only use *)
(* partial charges at atom centers. *)
(* They use normalized histograms, *)
(* from 1 to 13A with 12bins (1A discretization step) *)

open Printf
open Batteries

module MU = My_utils
module PQR = Pqr_parser
module PL = Pl_parser

let auto_correlation feature_space atoms: (float * float) list =
  (* sort values by increasing distance: this is assumed later on *)
  let sort =
    List.fast_sort (fun (d1, _) (d2, _) -> BatFloat.compare d1 d2)
  in
  let autocorr l =
    let rec loop acc = function
      | [] -> acc
      | a1 :: others ->
        let res = List.map (Atom.process_atoms feature_space a1) others in
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

let b84_histo feature_space atoms: (float * float) list =
  let ac = auto_correlation feature_space atoms in
  (* print_autocorr ac; *)
  histo 0.2 ac

let truncate_histo min_len curr_len histo =
  if curr_len > min_len then
    List.take min_len histo
  else
    histo

let triple_of_molecule m =
  Molecule.(m.name, m.index, m.atoms)

let read_molecules query_file db_file =
  let fmt_q = FilePath.get_extension query_file in
  let fmt_db = FilePath.get_extension db_file in
  if fmt_q <> fmt_db then begin
    Log.fatal "read_molecules: fmt_q <> fmt_db: %s %s" fmt_q fmt_db;
    exit 1
  end;
  match fmt_q with
  | "pqr" ->
    let query_molecule =
      match PQR.read_molecules query_file with
      | [molec] -> triple_of_molecule molec
      | _ -> failwith (sprintf "read_molecules: not one molecule in %s" query_file)
    in
    let db_molecules = PQR.read_molecules db_file in
    (Feature.Radius, query_molecule, db_molecules)
  | "pl" ->
    let query_molecule =
      match PL.read_molecules query_file with
      | [molec] -> triple_of_molecule molec
      | _ -> failwith (sprintf "read_molecules: not one molecule in %s" query_file)
    in
    let db_molecules = PL.read_molecules db_file in
    (Feature.LogP, query_molecule, db_molecules)
  | fmt ->
    let _ = Log.fatal "read_molecules: unsupported format: %s" fmt in
    exit 1

let main () =
  Log.set_log_level Log.INFO;
  Log.color_on();
  (* default option values *)
  let query_file = ref "" in
  let db_file = ref "" in
  let scores_file = ref "" in
  let method_str = ref "" in
  let cmd_line = Arg.align
      ["-q" , Arg.Set_string query_file , "query.{pqr|pl} query molecule";
       "-db", Arg.Set_string db_file    , "db.{pqr|pl} database";
       "-o" , Arg.Set_string scores_file, "out.scores score-labels file";
       "-m" , Arg.Set_string method_str , "{b84|m2005}"]
  in
  Arg.parse cmd_line ignore
    (sprintf "Example: %s -q query.mol2 -db database.mol2\n" Sys.argv.(0));
  if !query_file = "" || !db_file = "" || !scores_file = "" then begin
    Log.fatal "-q, -db and -o are all mandatory";
    exit 1
  end;
  let auc_file = !query_file ^ ".b84.auc" in
  let feature_space, (_, _, query_atoms), db_molecules =
    read_molecules !query_file !db_file
  in
  let query_histo = b84_histo feature_space query_atoms |> List.map snd in (* discard distances *)
  let query_histo_len = List.length query_histo in
  (* print_histo query_histo; *)
  MU.with_out_file !scores_file (fun out ->
      List.iter (fun molec ->
          let name, _i, atoms = triple_of_molecule molec in
          let histo = b84_histo feature_space atoms |> List.map snd in
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
    );
  let auc =
    ROC.read_AUC_from_string
      (MU.get_command_output
         (sprintf "croc-curve < %s 1> /dev/null 2> %s; cat %s"
            !scores_file auc_file auc_file))
  in
  Log.info "q: %s AUC: %f" !query_file auc

let () = main ()
