(* cache of distances between autocorrelation vectors of molecules *)

module AC = Autocorr
module HT = Hashtbl

type t = (int * int, float) HT.t

let create ncores autocorr_molecules: t =
  let res = HT.create 1000 in
  if ncores = 1 then
    begin
      let rec loop i ac_molecules =
        match ac_molecules with
        | [] -> ()
        | mol0 :: xs ->
          List.iteri (fun j mol1 ->
              let k = i + j + 1 in
              let dist = 1.0 -. (AC.tanimoto mol0 mol1) in
              HT.add res (i, k) dist
            ) xs;
          loop (i + 1) xs
      in
      loop 0 autocorr_molecules;
      res
    end
  else
    begin
      let rec prepare acc i ac_molecules =
        match ac_molecules with
        | [] -> acc
        | mol0 :: xs ->
          prepare ((i, mol0, xs) :: acc) (i + 1) xs
      in
      Log.info "preparing tripels";
      let prepared = prepare [] 0 autocorr_molecules in
      Log.info "parallel map";
      assert(ncores > 1);
      let triples_ll =
        Parmap.parmap ~ncores
          (fun (i, mol0, others) ->
             List.mapi
               (fun j mol1 ->
                  let k = i + j + 1 in
                  let dist = 1.0 -. (AC.tanimoto mol0 mol1) in
                  (i, k, dist)
               ) others
          ) (Parmap.L prepared)
      in
      Log.info "sequential reduce";
      List.iter (fun triples ->
          List.iter (fun (x, y, z) ->
              HT.add res (x, y) z
            ) triples
        ) triples_ll;
      res
    end

(* return distance between molecules at indexes 'i' and 'j'
   from the distance cache 'ht' *)
let get (ht: t) i j: float =
  if i = j then 0.0
  else HT.find ht (min i j, max i j)
