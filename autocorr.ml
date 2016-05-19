
open Batteries

module A  = Array
module At = Atom
module L  = List
module MU = My_utils

let process_atoms a1 a2 =
  At.(distance a1 a2, a1.charge *. a2.charge)

let auto_correlation atoms =
  let sort =
    L.fast_sort (fun (d1, _) (d2, _) -> Float.compare d1 d2)
  in
  let autocorr l =
    let rec loop acc = function
      | [] -> acc
      | a1 :: others ->
        let res = L.map (process_atoms a1) others in
        loop (L.rev_append res acc) others
    in
    loop [] l
  in
  let ac = autocorr atoms in
  let neg_ac, pos_ac =
    L.partition
      (fun (_, mep2) -> mep2 < 0.0)
      ac
  in
  (* sort values by increasing distance: this is assumed later on *)
  (sort neg_ac,
   sort pos_ac)

(* cross correlation for linearly binned auto correlation vectors *)
let cross_correlation a1 a2 =
  let n = min (A.length a1) (A.length a2) in
  let res = ref 0. in
  for i = 0 to n - 1 do
    res := !res +. ((A.unsafe_get a1 i) *. (A.unsafe_get a2 i))
  done;
  !res

let correlate_linbin_autocorrs
    (_i1, _name1, (n1, p1))
    (_i2, _name2, (n2, p2)) =
  0.5 *. ((cross_correlation n1 n2) +. (cross_correlation p1 p2))

let correlate_linbin_autocorrs'
    (_name1, (n1, p1))
    (_name2, (n2, p2)) =
  0.5 *. ((cross_correlation n1 n2) +. (cross_correlation p1 p2))

let score f a1 a2 =
  let n = min (A.length a1) (A.length a2) in
  let l1 = A.to_list (A.sub a1 0 n) in
  let l2 = A.to_list (A.sub a2 0 n) in
  f l1 l2

let tanimoto    = score MU.tanimoto_coeff
let tversky_ref = score MU.tversky_ref
let tversky_db  = score MU.tversky_db

let score_linbin_autocorrs f
    ((_i1 : int), (_name1 : string), (n1, p1))
    ((_i2 : int), (_name2 : string), (n2, p2)) =
  0.5 *. ((f n1 n2) +. (f p1 p2))

let tanimoto_linbin_autocorrs    = score_linbin_autocorrs tanimoto
let tversky_ref_linbin_autocorrs = score_linbin_autocorrs tversky_ref
let tversky_db_linbin_autocorrs  = score_linbin_autocorrs tversky_db

type autocorrelated = int * string * (float array * float array)

(* Plain Tanimoto score (and not average of two Tanimotos like before)
   Tanimoto(vec_x, vec_y) = xy_sum /. (x2_sum +. y2_sum -. xy_sum) *)
let tanimoto
    ((_i0, _name0, (neg_ac0, pos_ac0)): autocorrelated)
    ((_i1, _name1, (neg_ac1, pos_ac1)): autocorrelated): float =
  let sum_of_squares a =
    Array.fold_left
      (fun acc x -> acc +. (x *. x))
      0.0 a
  in
  let sum_of_products a1 a2 =
    let n = min (A.length a1) (A.length a2) in
    let res = ref 0.0 in
    for i = 0 to n - 1 do
      let x = a1.(i) in
      let y = a2.(i) in
      res := !res +. x *. y
    done;
    !res
  in
  let neg_ac0_2 = sum_of_squares neg_ac0 in
  let neg_ac1_2 = sum_of_squares neg_ac1 in
  let pos_ac0_2 = sum_of_squares pos_ac0 in
  let pos_ac1_2 = sum_of_squares pos_ac1 in
  let x2_sum = neg_ac0_2 +. pos_ac0_2 in
  let y2_sum = neg_ac1_2 +. pos_ac1_2 in
  let neg_xy_sum = sum_of_products neg_ac0 neg_ac1 in
  let pos_xy_sum = sum_of_products pos_ac0 pos_ac1 in
  let xy_sum = neg_xy_sum +. pos_xy_sum in
  xy_sum /. (x2_sum +. y2_sum -. xy_sum)
