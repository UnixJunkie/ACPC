
(* parallel (Parmap) version of the continuous version of ACPC *)

open Batteries
open Legacy.Printf
open Molecule

module AC   = Autocorr
module Feat = Feature
module L    = List
module MU   = My_utils
module PMU  = Common.Parmap_utils
module Q    = Query

let main () =

  let start_time = Unix.gettimeofday() in

  Log.set_log_level Log.INFO;
  Log.color_on();

  (* default options *)
  let nprocs = ref 1 in
  let not_set = 0.0 in
  let a = ref not_set in
  let database_fn = ref "" in
  let output_fn = ref "/dev/stdout" in
  let query_fn = ref "" in
  let cmd_line = Arg.align [
    "-a", Arg.Set_float a, "float kernel parameter (default was optimized on \
                            DUD-E and depends on considered feature space)";
    "-db", Arg.Set_string database_fn, "db.{mol2|pqr|pl}[.bin] molecules database";
    "-o", Arg.Set_string output_fn, "output_file to store name, index and scores (default is stdout)";
    "-q", Arg.Set_string query_fn, "query.{mol2|pqr|pl} query molecule";
    "-np", Arg.Set_int nprocs, "int max number of cores to use";
  ] in
  Arg.parse cmd_line ignore
    (sprintf "Example: %s -q query.mol2 -db database.mol2\n"
       Sys.argv.(0));
  if (!query_fn = "" || !database_fn = "") then (
    Log.fatalf "-q and -db are all mandatory";
    exit 1
  );
  Log.warnf "only scoring molecules, not removing duplicates";
  let feature = Feat.of_filename !query_fn in
  let db_feature = Feat.of_filename !database_fn in
  if feature <> db_feature then (
    Log.fatalf "query and database feature space don't match: %s %s"
      !query_fn !database_fn;
    exit 1
  );
  if !a = not_set then
    a := Feat.best_param feature
  ;
  let query_mol = AC.read_query_molecule feature !query_fn in
  let read_one_db_molecule = AC.get_reader !database_fn in
  let nb_molecules = ref 0 in
  (* Warning: this will blow up RAM in case the file is too large *)
  let all_molecules =
    MU.with_in_file !database_fn (fun input ->
      let molecules, exn =
        (* perf. hack: we don't preserve input order
           since we don't care about it *)
        MU.unfold_exc_rev (fun () ->
          read_one_db_molecule nb_molecules input
        )
      in
      assert(exn = End_of_file);
      molecules
    )
  in
  let query_ac = AC.auto_correlation feature query_mol.atoms in
  let indexed_query = Q.index_query !a query_ac in
  let map =
    if !nprocs > 1 then
      let csize = 1 in
      (* don't play with csize, it is almost useless *)
      PMU.list_parmap !nprocs csize
    else
      L.map
  in
  let name_index_scores =
    map
      (Q.query_for_parmap feature !a indexed_query)
      all_molecules
  in
  MU.with_out_file !output_fn (fun output ->
    L.iter
      (fun (name, index, score) ->
         fprintf output "%s %d %f\n" name index score)
      name_index_scores
  );
  let stop_time = Unix.gettimeofday() in
  let elapsed_time = stop_time -. start_time in
  Log.infof "speed: %.2f molecules/s"
    ((float_of_int !nb_molecules) /. elapsed_time);
;;

main()
