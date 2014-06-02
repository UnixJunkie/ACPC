
open Batteries

module S = String

type feature = Charge | Radius | LogP

let of_string = function
  | "pc"   -> Charge
  | "r"    -> Radius
  | "logP" -> LogP
  | s -> failwith ("Feature.of_string: unknown feature: " ^ s)

let to_string = function
  | Charge -> "pc"
  | Radius -> "r"
  | LogP   -> "logP"

let of_filename fn =
  if (S.ends_with fn ".mol2") ||
    (S.ends_with fn ".mol2.bin")
  then Charge
  else (
    if (S.ends_with fn ".pqr") ||
      (S.ends_with fn ".pqr.bin")
    then Radius
    else (
      if (S.ends_with fn ".pl") ||
        (S.ends_with fn ".pl.bin")
      then LogP
      else failwith ("Feature.of_filename: giving up: " ^ fn)
    )
  )

(* best default value for the 'a' parameter given a feature space
   - this was optimised on 51 (out of 102) targets of DUDE chosen randomly
     and using a single query per target
     (the first ligand in the ligands list) *)
let best_param = function
  | Charge -> 350.0
  | Radius -> 560.0
  | LogP   ->  42.0

(* chimera sphere properties used to highlight an atom in the given
   feature space as a triplet: (color, radius, transparency) *)
let chimera_sphere_properties = function
  | Charge -> "purple", 0.6, 0.5
  | Radius -> "gray"  , 0.5, 0.4
  | LogP   -> "green" , 0.4, 0.3
