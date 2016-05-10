(* cache of distances between autocorrelation vectors of molecules *)

module AC = Autocorr
module HT = Hashtbl

type t = (int * int, float) HT.t

(* FBR: TODO use parmap *)
let create autocorr_molecules: t =
  let res = HT.create 1000 in
  let rec loop i = function
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

(* return distance between molecules at indexes 'i' and 'j'
   from the distance cache 'ht' *)
let get (ht: t) i j: float =
  if i = j then 0.0
  else HT.find ht (min i j, max i j)
