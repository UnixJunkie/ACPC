
open Batteries

module A  = Array
module AC = Autocorr
module F  = File
module HT = Hashtbl
module L  = List
module Log = Dolog.Log
module MU = My_utils

open Printf

(* table for y values covering the whole range of auto_corr's x values *)
let xrange dx auto_corr =
  let d_min = fst (MU.hd auto_corr         ~err:"Kde.xrange.d_min") in
  let d_max = fst (MU.hd (L.rev auto_corr) ~err:"Kde.xrange.d_max") in
  let span = d_max +. d_min in (* a right margin is added, since there is
                                  a left one (from 0 to d_min) *)
  let nb_elts = 1 + (int_of_float (Legacy.ceil (span /. dx))) in
  A.make nb_elts 0.

(* Cf. Wand, M. P. (1994).
   "Fast Computation of Multivariate Kernel Estimators".
   Journal of Computational and Graphical Statistics,
   Vol. 3, Num. 4, 433-445. *)
let linear_binning dx auto_corr =
  let res = xrange dx auto_corr in
  L.iter
    (fun (y, mep2) ->
       let i = int_of_float (y /. dx) in
       let j = i + 1 in
       let x = (float_of_int i) *. dx in
       let z = x +. dx in
       let lw = (z -. y) /. dx in
       let rw = 1. -. lw in
       let prev_lval = A.unsafe_get res i in
       let prev_rval = A.unsafe_get res j in
       A.unsafe_set res i (prev_lval +. (lw *. mep2));
       A.unsafe_set res j (prev_rval +. (rw *. mep2)))
    auto_corr;
  res
(*$T linear_binning
  linear_binning 0.5 [(0., 0.), (1., 0.), (2., 0.)] = A.make 4 0.
*)

let linbin_autocorr maybe_debug_fn dx (mol_name, _index, atoms) =
  let dumb_ac = [(0., 0.)] in
  let protect ac logit =
    if ac = [] then
      let _ = logit() in
      dumb_ac
    else
      ac
  in
  let neg_ac', pos_ac' = AC.auto_correlation atoms in
  let neg_ac , pos_ac =
    protect neg_ac'
      (fun _ -> Log.warn "empty neg_ac in %s" mol_name),
    protect pos_ac'
      (fun _ -> Log.warn "empty pos_ac in %s" mol_name)
  in
  let lb_neg_ac = linear_binning dx neg_ac in
  let lb_pos_ac = linear_binning dx pos_ac in
  (match maybe_debug_fn with
    | None -> ()
    | Some fn ->
      let neg_ac_fn = fn ^ ".nac" in
      let pos_ac_fn = fn ^ ".pac" in
      let neg_lb_ac_fn = fn ^ ".nlbac" in
      let pos_lb_ac_fn = fn ^ ".plbac" in
      MU.list_to_file
        neg_ac_fn
        (fun (d, mep2) -> sprintf "%f %f\n" d mep2)
        neg_ac;
      MU.list_to_file
        pos_ac_fn
        (fun (d, mep2) -> sprintf "%f %f\n" d mep2)
        pos_ac;
      F.with_file_out neg_lb_ac_fn (fun out ->
        A.iteri
          (fun i x -> fprintf out "%f %f\n" (dx *. float_of_int i) x)
          lb_neg_ac);
      F.with_file_out pos_lb_ac_fn (fun out ->
        A.iteri
          (fun i x -> fprintf out "%f %f\n" (dx *. float_of_int i) x)
          lb_pos_ac)
  );
  (mol_name,
   (lb_neg_ac, lb_pos_ac))
