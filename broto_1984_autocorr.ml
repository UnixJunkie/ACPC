
(* trial at reproducing method from  *)
(* "Molecular structures: perception, autocorrelation descriptor and sar studies" *)
(* Pierre Broto, Gilles Moreau and Corinne Vandycke. *)
(* Eur. J. Med. Chem. 1984-19 N1 pp 66-70 *)

(* we can reproduce easily their contribution to LogP and steric *)
(* autocorrelation vectors *)
(* they use histograms without normalization with discretization step = 0.2A *)

module L = List

let auto_correlation (atoms: Atom.atom list) =
  (* sort values by increasing distance: this is assumed later on *)
  let sort =
    L.fast_sort (fun (d1, _) (d2, _) -> BatFloat.compare d1 d2)
  in
  let autocorr l =
    let rec loop acc = function
      | [] -> acc
      | a1 :: others ->
        let res = L.map (Autocorr.process_atoms a1) others in
        loop (L.rev_append res acc) others
    in
    loop [] l
  in
  sort (autocorr atoms)

let b84_histo (atoms: Atom.atom list) =
  let acpc_ac = auto_correlation atoms in
  failwith "not implemented yet"
