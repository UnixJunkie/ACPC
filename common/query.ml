
open Batteries
open Molecule

module AC = Autocorr
module HT = Hashtbl
module L  = BatList
module MU = My_utils
module S  = BatString

let query_for_parmap feature a indexed_query candidate =
  let score = AC.continuous_overlap feature a indexed_query candidate in
  (* the name is superfluous since the father process can know the mapping
     between names and indexes, however doing it here puts more work into
     the parallel section of the program so its good for performances
     in fact ! \(^v^)/ *)
  (candidate.name, candidate.index, score)

let is_active_by_name mol_name =
  S.starts_with mol_name "active"

let is_active_by_ht ht mol_name =
  HT.mem ht mol_name

let query_w_indexed_molecule
    maybe_actives_fn
    feature a database_fn read_one_db_molecule nb_molecules indexed_queries =
  let actives = match maybe_actives_fn with
    | None -> HT.create 0
    | Some actives_fn -> 
      let ht = HT.create 500 in
      MU.iter_on_lines_of_file
        (fun mol_name -> HT.add ht mol_name ())
        actives_fn;
      ht
  in
  let is_active = match maybe_actives_fn with
    | None -> is_active_by_name
    | Some _ -> is_active_by_ht actives
  in
  (* compute score-labels of all DB molecules against the query *)
  MU.with_in_file database_fn (fun input ->
    let score_labels, exn =
      L.unfold_exc
        (fun () ->
           let i = !nb_molecules in
           let candidate = read_one_db_molecule nb_molecules input in
           let name = candidate.name in
           let score =
             (* the final score is the sum of all fragments' scores (in case
                indexed_queries has more than one element) *)
             L.fsum
               (L.map
                  (fun indexed_query ->
                     AC.continuous_overlap feature a indexed_query candidate)
                  indexed_queries
               )
           in
           (name, score, i, is_active name)
        )
    in
    assert(exn = End_of_file);
    score_labels
  )

(* index the values of an autocorrelated query molecule *)
let index_query a (q_neg_ac, q_pos_ac) =
  (AC.index_ac a q_neg_ac,
   AC.index_ac a q_pos_ac)

(* query_molecules is either one molecule or the list of fragments of the
   initial query molecule *)
let do_query
    maybe_actives_fn
    feature a database_fn read_one_db_molecule nb_molecules query_molecules =
  let indexed_queries =
    L.map
      (fun m -> let query_ac = AC.auto_correlation feature m.atoms in
                index_query a query_ac
      )
      query_molecules
  in
  query_w_indexed_molecule
    maybe_actives_fn
    feature a database_fn read_one_db_molecule nb_molecules indexed_queries
