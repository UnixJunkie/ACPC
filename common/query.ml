
open Batteries
open Molecule

module AC = Autocorr
module L  = List
module MU = My_utils
module S  = BatString

let query_for_parmap feature a indexed_query candidate =
  let score = AC.continuous_overlap feature a indexed_query candidate in
  (* the name is superfluous since the father process can know the mapping
     between names and indexes, however doing it here puts more work into
     the parallel section of the program so its good for performances
     in fact ! \(^v^)/ *)
  (candidate.name, candidate.index, score)

let query_w_indexed_molecule
    feature a database_fn read_one_db_molecule nb_molecules indexed_queries =
  (* compute score-labels of all DB molecules against the query *)
  MU.with_in_file database_fn (fun input ->
    let score_labels, exn =
      MU.unfold_exc
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
           let label = S.starts_with name "active" in
           (name, score, i, label)
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
    feature a database_fn read_one_db_molecule nb_molecules query_molecules =
  let indexed_queries =
    L.map
      (fun m -> let query_ac = AC.auto_correlation feature m.atoms in
                index_query a query_ac
      )
      query_molecules
  in
  query_w_indexed_molecule
    feature a database_fn read_one_db_molecule nb_molecules indexed_queries
