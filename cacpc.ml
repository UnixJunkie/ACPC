
(* continuous version of ACPC, previous one was discrete *)

open Batteries
open Legacy.Printf
open Molecule

module AC    = Autocorr
module Color = Chimera_colors
module Feat  = Feature
module Fn    = Filename
module L     = List
module MU    = My_utils
module Mol   = Molecule
module Mol2  = Mol2_parser
module MolG  = Molecular_graph
module Q     = Query

let main () =
  let start_time = Unix.gettimeofday() in
  Log.set_log_level Log.INFO;
  Log.color_on();
  (* default option values *)
  let not_set_float = 0.0 in
  let not_set_int = 0 in
  let top_N = ref not_set_int in
  let a = ref not_set_float in
  let er_ratio = ref not_set_float in
  let database_fn = ref "" in
  let actives_fn = ref "" in
  let dump_scores = ref false in
  let break_RBNs = ref false in
  let post_filter = ref true in
  let query_fn = ref "" in
  let scan_out_fn = ref "" in
  let out_fn = ref "" in
  let no_hydrogens = ref false in
  let cmd_line = Arg.align [
    "-a", Arg.Set_float a, "float kernel parameter (default depends \
                            on considered feature space and was optimized \
                            on DUD-E)";
    "-actives", Arg.Set_string actives_fn, "some_file name of active molecules \
                                            (one per line)";
    "-brbn", Arg.Set break_RBNs,
    " break the query molecule along its rotatable bonds";
    "-db", Arg.Set_string database_fn, "db.{mol2|pqr|pl}[.bin] database";
    "-er", Arg.Set_float er_ratio,
    "x enrichment rate parameter; e.g. 0.01 --> ER_1%";
    "-nopp", Arg.Clear post_filter,
    " don't rm duplicate molecules (based on names)";
    "-noH", Arg.Set no_hydrogens, " ignore hydrogen atoms in mol2 file";
    "-q", Arg.Set_string query_fn, "query.{mol2|pqr|pl} query";
    "-scan", Arg.Set_string scan_out_fn,
    "out_file scan delta AUC per atom in \
      the query and create a new query molecule";
    "-top", Arg.Set_int top_N, "int number of top scoring molecules to output \
                                (only works for MOL2 files)";
    "-o", Arg.Set_string out_fn, "out_file output file for -top";
    "-v", Arg.Set dump_scores, " dump scores out"
  ] in
  Arg.parse cmd_line ignore
    (sprintf "Example: %s -q query.mol2 -db database.mol2[.bin]\n"
       Sys.argv.(0));
  if (!query_fn = "" || !database_fn = "") then (
    Log.fatal "-q and -db are all mandatory";
    exit 1
  );
  if !top_N <> not_set_int then (
    if !out_fn = "" then (
      Log.fatal "-top implies -o";
      exit 1
    );
  );
  if !out_fn <> "" then (
    if !top_N = not_set_int then (
      Log.fatal "-o but no -top";
      exit 1
    );
  );
  let feature = Feat.of_filename !query_fn in
  let db_feature = Feat.of_filename !database_fn in
  if feature <> db_feature then (
    Log.fatal "query and database feature don't match: %s %s"
      !query_fn !database_fn;
    exit 1
  );
  if !a = not_set_float then
    a := Feat.best_param feature
  ;
  let query_molecules =
    if not !break_RBNs then
      [AC.read_query_molecule feature !query_fn]
    else (
      let rot_bonds = MolG.get_rot_bonds !query_fn in
      Log.info "found %d RBNs" (L.length rot_bonds);
      let rot_bond_atoms, fragments = MolG.cut_rot_bonds !query_fn in
      Log.info "found %d fragments" (L.length fragments);
      (* dump fragments and rot. bonds out for verification *)
      let rot_bonds_bild_fn =
        MU.filename_with_different_extension !query_fn ".mol2" "_bonds.bld"
      in
      Mol.rot_bonds_to_bild_file rot_bond_atoms rot_bonds_bild_fn;
      Log.info "rot. bonds: %s" rot_bonds_bild_fn;
      let current_color = ref Color.Red in
      L.iteri
        (fun i m ->
          let bild_file =
            Fn.concat (Fn.dirname !query_fn) (sprintf "%s.bld" m.name) in
          Log.info "fragment %d: %s" i bild_file;
          Mol.to_bild_file m !current_color bild_file;
          current_color := Color.next !current_color
        )
        fragments;
      fragments
    )
  in
  let read_one_db_molecule = AC.get_reader !database_fn in
  let nb_molecules = ref 0 in
  let maybe_actives_fn = match !actives_fn with
    | "" -> None
    | fn -> Some fn
  in
  let score_labels =
    let scores =
      Q.do_query
        maybe_actives_fn
        feature !a !database_fn read_one_db_molecule nb_molecules
        query_molecules
    in
    if !post_filter then
      ROC.unique scores
    else (
      Log.warn "not removing duplicated molecules";
      scores
    )
  in
  let stop_time = Unix.gettimeofday() in
  let elapsed_time = stop_time -. start_time in
    (* we output speed before doing the AUC calculation
       because AUC is an optional post processing step *)
  Log.info "speed: %.2f molecules/s"
    ((float_of_int !nb_molecules) /. elapsed_time);
  let initial_query_auc = ROC.auc score_labels in
  let enr_rate = ref nan in
  if !er_ratio <> not_set_float then (
    assert(!er_ratio > 0.0 && !er_ratio <= 1.0);
    let _top_n, _top_actives_rate, _rand_actives_rate, er =
      ROC.enr_rate !er_ratio score_labels
    in
    enr_rate := er
  );
  Log.info "q: %s db: %s a: %.1f AUC: %.3f ER@%.2f: %.2f"
    !query_fn
    !database_fn
    !a
    initial_query_auc
    !er_ratio
    !enr_rate
  ;
  if !dump_scores then (
    (* create .scores file for rank-based consensus queries *)
    let out_fn = !query_fn ^ ".scores" in
    MU.with_out_file out_fn (fun out ->
      L.iter
        (fun (name, score, index, _label) ->
          fprintf out "%s %f %d\n" name score index)
        score_labels
    );
  );
  if !scan_out_fn <> "" then (
    match query_molecules with
    | [query_mol] -> (* atom-grained scanning of the query molecule *)
      let nb_AUCs =
        AC.atom_scan
          feature initial_query_auc (!query_fn, query_mol) !scan_out_fn
          (fun q_mol ->
             let curr_score_labels =
               Q.do_query
                 (* (ref 0) means we ignore the number of DB molecules *)
                 maybe_actives_fn
                 feature !a !database_fn read_one_db_molecule (ref 0) [q_mol]
             in
             ROC.auc curr_score_labels
          )
      in
      assert(nb_AUCs = L.length query_mol.atoms)
    | _ ->
      Log.fatal "you are not supposed to scan fragments";
      exit 1
  );
  if !top_N <> not_set_int && !out_fn <> "" then (
    Log.info "writing top %d molecules to %s..." !top_N !out_fn;
    (* put all molecules in a table *)
    let indexed_molecules =
      (match feature with
      | Feat.Charge -> Mol2.read_raw_molecules
      | Feat.Radius
      | Feat.LogP -> failwith "read_raw_molecules not implemented")
        !database_fn
    in
    let top_scoring_molecules_first =
      L.fast_sort
        (fun
          (_name1, score1, _index1, _label1)
          (_name2, score2, _index2, _label2) ->
            (* high scores first *)
            Float.compare score2 score1
        )
        score_labels
    in
    let top_N_molecules = L.take !top_N top_scoring_molecules_first in
    (* output top N best ones *)
    MU.with_out_file !out_fn (fun out ->
      L.iter
        (fun (_name, _score, i, _label) ->
           let molecule_lines = indexed_molecules.(i) in
           L.iter (fprintf out "%s\n") molecule_lines
        )
        top_N_molecules
    );
  )
;;

main()
