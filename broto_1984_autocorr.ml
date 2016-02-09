
(* trial at reproducing method from  *)
(* "Molecular structures: perception, autocorrelation descriptor and sar studies" *)
(* Pierre Broto, Gilles Moreau and Corinne Vandycke. *)
(* Eur. J. Med. Chem. 1984-19 N1 pp 66-70 *)

(* we can reproduce easily their contribution to LogP and steric *)
(* autocorrelation vectors *)
(* they use histograms without normalization with discretization step = 0.2A *)

open Batteries

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

let histo dx pairs =
  let rec loop acc to_process =
    match to_process with
    | [] -> List.rev acc
    | _ ->
      failwith "not implemented yet"
  in
  loop [] pairs

let b84_histo (atoms: Atom.atom list): float * float list =
  let acpc_ac = auto_correlation atoms in
  failwith "not implemented yet"
