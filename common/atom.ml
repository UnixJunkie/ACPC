
open Batteries

module AA   = Amino_acid
module Feat = Feature
module L    = List
module MU   = My_utils
module S    = String
module V3   = Vector3

open Printf

(* To avoid programming and manipulation errors, we consider a single
   property at a time:
   van der vals radius | partial charge | logP contribution *)
type property = Radius of float | Charge of float | LogP of float

type atom = { position: V3.t ; (* x, y, z in the PDB *)
              prop: property }

(* ACCESORS *)

let position a = a.position

let radius on a = match on with
  | Feat.Radius ->
    (match a.prop with
      | Radius r -> r
      | Charge _ -> failwith "property type mismatch: radius Vs. charge"
      | LogP _ -> failwith "property type mismatch: radius Vs. LogP"
    )
  | Feat.Charge -> failwith "'Atom.radius on a' called with on = Charge"
  | Feat.LogP -> failwith "'Atom.radius on a' called with on = LogP"

let charge on a = match on with
  | Feat.Radius -> failwith "'Atom.charge on a' called with on = Radius"
  | Feat.Charge ->
    (match a.prop with
      | Radius _ -> failwith "property type mismatch: charge Vs. radius"
      | Charge q -> q
      | LogP _ -> failwith "property type mismatch: charge Vs. LogP"
    )
  | Feat.LogP -> failwith "'Atom.charge on a' called with on = LogP"

let logP on a = match on with
  | Feat.Radius -> failwith "'Atom.logP on a' called with on = Radius"
  | Feat.Charge -> failwith "'Atom.logP on a' called with on = Charge"
  | Feat.LogP -> (match a.prop with
      | Radius _ -> failwith "property type mismatch: logP Vs. radius"
      | Charge _ -> failwith "property type mismatch: logP Vs. charge"
      | LogP lp -> lp
  )

let create x y z prop =
  { position = V3.make x y z ;
    prop }

let one_or_more_spaces = Str.regexp "[ ]+"

let a2f = float_of_string

let of_pl_line l =
  let tokens = Str.split one_or_more_spaces l in
  let res =
    match tokens with
        [ _atom ; xs ; ys ; zs ; logPs ] ->
          (try Some { position = V3.make (a2f xs) (a2f ys) (a2f zs) ;
                      prop     = LogP (a2f logPs)                   }
           with _ -> None)
      | _ -> None
  in
  match res with
      Some a -> a
    | None -> failwith ("atom_of_pl_line: could not parse: " ^ l)

let of_pqr_line l =
  (* end of a pqr line from Open Babel 2.3.0: x y z q r element_name *)
  let tokens = L.rev (Str.split one_or_more_spaces l) in
  let res =
    match tokens with
        _element :: rs :: _qs :: zs :: ys :: xs :: _ ->
          (try Some { position = V3.make (a2f xs) (a2f ys) (a2f zs) ;
                      prop     = Radius (a2f rs)                    }
           with _ -> None)
      | _ -> None
  in
  match res with
      Some a -> a
    | None -> failwith ("atom_of_pqr_line: could not parse: " ^ l)

(* example line:
"      1 C1          9.2324   -6.5321    2.1574 C.3       1 <0>        -0.0750"
*)
let of_mol2_line l =
  try Scanf.sscanf l " %d %s %f %f %f %s %d %s %f"
        (fun _id _name x y z _type _subst_id _subst_name q ->
          { position = V3.make x y z ;
            prop     = Charge q      }
        )
  with _ -> failwith ("atom.ml: of_mol2_line: could not parse: " ^ l)

(* same as of_mol2_line except that hydrogen atoms are ignored *)
let of_mol2_line_no_H l =
  try Scanf.sscanf l " %d %s %f %f %f %s %d %s %f"
        (fun _id _name x y z atom_type _subst_id _subst_name q ->
           if atom_type = "H" then None
           else Some { position = V3.make x y z ;
                       prop     = Charge q      }
        )
  with _ -> failwith ("atom.ml: of_mol2_line_no_H: could not parse: " ^ l)

let bond_of_mol2_line l =
  try Scanf.sscanf l " %d %d %d %s"
        (fun bond_i origin_atom_i target_atom_i _bond_type ->
          (bond_i, origin_atom_i, target_atom_i)
        )
  with _ -> failwith ("atom.ml: bond_of_mol2_line: could not parse: " ^ l)

let to_string feature a =
  V3.(sprintf "%f %f %f %f"
        a.position.x a.position.y a.position.z
        (match feature with
          | Feat.Charge -> charge Feat.Charge a
          | Feat.Radius -> radius Feat.Radius a
          | Feat.LogP   -> logP   Feat.LogP a
        ))

let pseudo_mol2_line i a =
  (* OB's MOL2 format code:
     snprintf(buffer,BUFF_SIZE,"%7d%1s%-6s%12.4f%10.4f%10.4f%1s%-5s%4s%1s %-8s%10.4f",
     atom->GetIdx(),"",label,
     atom->GetX(),atom->GetY(),atom->GetZ(),
     "",str1.c_str(),
     rnum,"",rlabel,
     atom->GetPartialCharge());
  *)
  V3.(sprintf "%7d %-6s %12.4f%10.4f%10.4f %-5s%4s  %-8s%10.4f"
        i "C0" a.position.x a.position.y a.position.z
        "C.ar" "1" "<0>11" (charge Feat.Charge a)
  )

(* .pl file line for an atom *)
let pl_line a =
  (* eg: 'ATOM 3.906 2.7422 0.9614 0.2811' *)
  V3.(sprintf "ATOM %f %f %f %f"
        a.position.x
        a.position.y
        a.position.z
        (logP Feat.LogP a))

let pseudo_pqr_line i a =
  (* eg: "ATOM      1  C0  <0> A   1       1.236  -4.057   3.848 -0.09430000   1.700  C  " *)
  V3.(sprintf
        "ATOM    %03d  C0  <0> A   1    %8.3f%8.3f%8.3f  0.00000000%8.3f  C  "
        i
        a.position.x
        a.position.y
        a.position.z
        (radius Feat.Radius a))

let move atom displacement =
  { atom with position = V3.add atom.position displacement }

(* (\* add r to radius *\) *)
(* let dilate a r = *)
(*   { a with radius = a.radius +. r } *)

let distance atom1 atom2 =
  V3.dist atom1.position atom2.position

let distance_to_point a p =
  V3.mag (V3.diff a.position p)

let squared_distance_to_point a p =
  V3.mag2 (V3.diff a.position p)

(* let point_is_in_atom p a = *)
(*   squared_distance_to_point a p <= a.radius *. a.radius *)

(* let point_in_any_atom p atoms = *)
(*   L.exists (point_is_in_atom p) atoms *)

let point_within_d2_of_atom p d2 a =
  squared_distance_to_point a p <= d2

(* at some point there is a product of properties when autocorrelating *)
let process_atoms on a1 a2 =
  distance a1 a2,
  match on with
    | Feat.Charge -> charge on a1 *. charge on a2
    | Feat.LogP   -> logP   on a1 *. logP   on a2
    | Feat.Radius -> radius on a1 *. radius on a2

(* (\* do they overlap? *\) *)
(* let clash atom1 atom2 = *)
(*   let r1 = atom1.radius in *)
(*   let r2 = atom2.radius in *)
(*   distance atom1 atom2 <= r1 +. r2 *)

(* (\* return xmin, xmax, ymin, ymax, zmin, zmax *\) *)
(* let bounding_box a = *)
(*   let x, y, z = V3.to_triplet a.position in *)
(*   let r = a.radius in *)
(*   (x -. r, x +. r, *)
(*    y -. r, y +. r, *)
(*    z -. r, z +. r) *)

(* (\* list surface points *\) *)
(* let dot_surface a unit_vectors = *)
(*   let scale_then_translate v = *)
(*     V3.add a.position (V3.mult v a.radius) *)
(*   in *)
(*   L.map scale_then_translate unit_vectors *)

(* (\* surface of a sphere *\) *)
(* let surface a = *)
(*   let pi = MU.pi in *)
(*   let r2 = a.radius *. a.radius in *)
(*   4. *. pi *. r2 *)

(* (\* r is in (points / A^2) *\) *)
(* let surface_at_resolution a r = *)
(*   let s = surface a in *)
(*   let nb_points = int_of_float (ceil (r *. s)) in *)
(*   dot_surface a (Spir.unit_vectors nb_points) *)

let to_xyz a =
  V3.to_triplet a.position

let to_mol2_string a =
  let x, y, z = to_xyz a in
  let q = charge Feat.Charge a in
  Printf.sprintf "%f %f %f %f"
    x
    y
    z
    q

let to_pqr_string a =
  let x, y, z = to_xyz a in
  let r = radius Feat.Radius a in
  Printf.sprintf "%f %f %f %f"
    x
    y
    z
    r

let to_pl_string a =
  let x, y, z = to_xyz a in
  let lp = logP Feat.LogP a in
  Printf.sprintf "%f %f %f %f"
    x
    y
    z
    lp

(* let of_string s = *)
(*   Scanf.sscanf s "%f %f %f %f %f" *)
(*     (fun x y z r q -> { position = V3.make x y z ; *)
(*                         radius   = r ; *)
(*                         charge   = q }) *)

(* reference:
   http://www-cryst.bioc.cam.ac.uk/~richard/piccolo/atom_table.php *)
type atom_class = Hydrophobic | Aromatic | Cationic | Anionic | HbD | HbA

(* atom types found in PDB files *)
type atom_type =
  | C
  | CA
  | CB
  | CD
  | CD1
  | CD2
  | CE
  | CE1
  | CE2
  | CE3
  | CG
  | CG1
  | CG2
  | CH2
  | CZ
  | CZ2
  | CZ3
  | N
  | ND1
  | ND2
  | NE
  | NE1
  | NE2
  | NH1
  | NH2
  | NZ
  | O
  | OD1
  | OD2
  | OE1
  | OE2
  | OG
  | OG1
  | OH
  | SD
  | SG

exception Unknown_atom_type of string

let atom_type_of_string = function
 | "C" -> C
 | "CA" -> CA
 | "CB" -> CB
 | "CD" -> CD
 | "CD1" -> CD1
 | "CD2" -> CD2
 | "CE" -> CE
 | "CE1" -> CE1
 | "CE2" -> CE2
 | "CE3" -> CE3
 | "CG" -> CG
 | "CG1" -> CG1
 | "CG2" -> CG2
 | "CH2" -> CH2
 | "CZ" -> CZ
 | "CZ2" -> CZ2
 | "CZ3" -> CZ3
 | "N" -> N
 | "ND1" -> ND1
 | "ND2" -> ND2
 | "NE" -> NE
 | "NE1" -> NE1
 | "NE2" -> NE2
 | "NH1" -> NH1
 | "NH2" -> NH2
 | "NZ" -> NZ
 | "O" -> O
 | "OD1" -> OD1
 | "OD2" -> OD2
 | "OE1" -> OE1
 | "OE2" -> OE2
 | "OG" -> OG
 | "OG1" -> OG1
 | "OH" -> OH
 | "SD" -> SD
 | "SG" -> SG
 | unk -> raise (Unknown_atom_type unk)

exception Unknown_aa_atom_pair of AA.aa_three * atom_type

(* Residue, Atom ->
   [Hydrophobic; Aromatic; Cationic; Anionic; HbD; HbA], VdWr (* Covr *)
   source: http://www-cryst.bioc.cam.ac.uk/~richard/piccolo/atom_table.php
   HID, HIE and HIP were provided by Arnout Voet *)
let classify_atom = function
 | AA.ALA, C   -> [], 1.61 (* 0.77 *)
 | AA.ALA, CA  -> [], 1.88 (* 0.77 *)
 | AA.ALA, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.ALA, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.ALA, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.ARG, C   -> [], 1.61 (* 0.77 *)
 | AA.ARG, CA  -> [], 1.88 (* 0.77 *)
 | AA.ARG, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.ARG, CD  -> [], 1.88 (* 0.77 *)
 | AA.ARG, CG  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.ARG, CZ  -> [Cationic], 1.61 (* 0.77 *)
 | AA.ARG, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.ARG, NE  -> [Cationic; HbD], 1.64 (* 0.70 *)
 | AA.ARG, NH1 -> [Cationic; HbD], 1.64 (* 0.70 *)
 | AA.ARG, NH2 -> [Cationic; HbD], 1.64 (* 0.70 *)
 | AA.ARG, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.ASN, C   -> [], 1.61 (* 0.77 *)
 | AA.ASN, CA  -> [], 1.88 (* 0.77 *)
 | AA.ASN, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.ASN, CG  -> [], 1.61 (* 0.77 *)
 | AA.ASN, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.ASN, ND2 -> [HbD], 1.64 (* 0.70 *)
 | AA.ASN, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.ASN, OD1 -> [HbA], 1.42 (* 0.66 *)
 | AA.ASP, C   -> [], 1.61 (* 0.77 *)
 | AA.ASP, CA  -> [], 1.88 (* 0.77 *)
 | AA.ASP, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.ASP, CG  -> [Anionic], 1.61 (* 0.77 *)
 | AA.ASP, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.ASP, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.ASP, OD1 -> [Anionic; HbA], 1.42 (* 0.66 *)
 | AA.ASP, OD2 -> [Anionic; HbA], 1.42 (* 0.66 *)
 | AA.CYS, C   -> [], 1.61 (* 0.77 *)
 | AA.CYS, CA  -> [], 1.88 (* 0.77 *)
 | AA.CYS, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.CYS, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.CYS, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.CYS, SG  -> [HbD; HbA], 1.77 (* 1.04 *)
 | AA.GLN, C   -> [], 1.61 (* 0.77 *)
 | AA.GLN, CA  -> [], 1.88 (* 0.77 *)
 | AA.GLN, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.GLN, CD  -> [], 1.61 (* 0.77 *)
 | AA.GLN, CG  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.GLN, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.GLN, NE2 -> [HbD], 1.64 (* 0.70 *)
 | AA.GLN, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.GLN, OE1 -> [HbA], 1.42 (* 0.66 *)
 | AA.GLU, C   -> [], 1.61 (* 0.77 *)
 | AA.GLU, CA  -> [], 1.88 (* 0.77 *)
 | AA.GLU, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.GLU, CD  -> [Anionic], 1.61 (* 0.77 *)
 | AA.GLU, CG  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.GLU, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.GLU, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.GLU, OE1 -> [Anionic; HbA], 1.42 (* 0.66 *)
 | AA.GLU, OE2 -> [Anionic; HbA], 1.42 (* 0.66 *)
 | AA.GLY, C   -> [], 1.61 (* 0.77 *)
 | AA.GLY, CA  -> [], 1.88 (* 0.77 *)
 | AA.GLY, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.GLY, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.HIS, C   -> [], 1.61 (* 0.77 *)
 | AA.HIS, CA  -> [], 1.88 (* 0.77 *)
 | AA.HIS, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.HIS, CD2 -> [Aromatic; Cationic], 1.76 (* 0.77 *)
 | AA.HIS, CE1 -> [Aromatic; Cationic], 1.76 (* 0.77 *)
 | AA.HIS, CG  -> [Aromatic; Cationic], 1.61 (* 0.77 *)
 | AA.HIS, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.HIS, ND1 -> [Aromatic; Cationic; HbD; HbA], 1.64 (* 0.70 *)
 | AA.HIS, NE2 -> [Aromatic; Cationic; HbD; HbA], 1.64 (* 0.70 *)
 | AA.HIS, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.HID, C   -> [], 1.61 (* 0.77 *)
 | AA.HID, CA  -> [], 1.88 (* 0.77 *)
 | AA.HID, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.HID, CD2 -> [Aromatic; Cationic], 1.76 (* 0.77 *)
 | AA.HID, CE1 -> [Aromatic; Cationic], 1.76 (* 0.77 *)
 | AA.HID, CG  -> [Aromatic; Cationic], 1.61 (* 0.77 *)
 | AA.HID, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.HID, ND1 -> [Aromatic; Cationic; HbD], 1.64 (* 0.70 *)
 | AA.HID, NE2 -> [Aromatic; HbA], 1.64 (* 0.70 *)
 | AA.HID, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.HIE, C   -> [], 1.61 (* 0.77 *)
 | AA.HIE, CA  -> [], 1.88 (* 0.77 *)
 | AA.HIE, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.HIE, CD2 -> [Aromatic; Cationic], 1.76 (* 0.77 *)
 | AA.HIE, CE1 -> [Aromatic; Cationic], 1.76 (* 0.77 *)
 | AA.HIE, CG  -> [Aromatic; Cationic], 1.61 (* 0.77 *)
 | AA.HIE, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.HIE, ND1 -> [Aromatic; HbA], 1.64 (* 0.70 *)
 | AA.HIE, NE2 -> [Aromatic; Cationic; HbD], 1.64 (* 0.70 *)
 | AA.HIE, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.HIP, C   -> [], 1.61 (* 0.77 *)
 | AA.HIP, CA  -> [], 1.88 (* 0.77 *)
 | AA.HIP, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.HIP, CD2 -> [Aromatic; Cationic], 1.76 (* 0.77 *)
 | AA.HIP, CE1 -> [Aromatic; Cationic], 1.76 (* 0.77 *)
 | AA.HIP, CG  -> [Aromatic; Cationic], 1.61 (* 0.77 *)
 | AA.HIP, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.HIP, ND1 -> [Aromatic; Cationic; HbD], 1.64 (* 0.70 *)
 | AA.HIP, NE2 -> [Aromatic; Cationic; HbD], 1.64 (* 0.70 *)
 | AA.HIP, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.ILE, C   -> [], 1.61 (* 0.77 *)
 | AA.ILE, CA  -> [], 1.88 (* 0.77 *)
 | AA.ILE, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.ILE, CD1 -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.ILE, CG1 -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.ILE, CG2 -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.ILE, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.ILE, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.LEU, C   -> [], 1.61 (* 0.77 *)
 | AA.LEU, CA  -> [], 1.88 (* 0.77 *)
 | AA.LEU, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.LEU, CD1 -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.LEU, CD2 -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.LEU, CG  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.LEU, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.LEU, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.LYS, C   -> [], 1.61 (* 0.77 *)
 | AA.LYS, CA  -> [], 1.88 (* 0.77 *)
 | AA.LYS, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.LYS, CD  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.LYS, CE  -> [], 1.88 (* 0.77 *)
 | AA.LYS, CG  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.LYS, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.LYS, NZ  -> [Cationic; HbD], 1.64 (* 0.70 *)
 | AA.LYS, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.MET, C   -> [], 1.61 (* 0.77 *)
 | AA.MET, CA  -> [], 1.88 (* 0.77 *)
 | AA.MET, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.MET, CE  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.MET, CG  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.MET, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.MET, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.MET, SD  -> [Hydrophobic; HbA], 1.77 (* 1.04 *)
 | AA.PHE, C   -> [], 1.61 (* 0.77 *)
 | AA.PHE, CA  -> [], 1.88 (* 0.77 *)
 | AA.PHE, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.PHE, CD1 -> [Hydrophobic; Aromatic], 1.76 (* 0.77 *)
 | AA.PHE, CD2 -> [Hydrophobic; Aromatic], 1.76 (* 0.77 *)
 | AA.PHE, CE1 -> [Hydrophobic; Aromatic], 1.76 (* 0.77 *)
 | AA.PHE, CE2 -> [Hydrophobic; Aromatic], 1.76 (* 0.77 *)
 | AA.PHE, CG  -> [Hydrophobic; Aromatic], 1.61 (* 0.77 *)
 | AA.PHE, CZ  -> [Hydrophobic; Aromatic], 1.76 (* 0.77 *)
 | AA.PHE, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.PHE, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.PRO, C   -> [], 1.61 (* 0.77 *)
 | AA.PRO, CA  -> [], 1.88 (* 0.77 *)
 | AA.PRO, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.PRO, CD  -> [], 1.88 (* 0.77 *)
 | AA.PRO, CG  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.PRO, N   -> [], 1.64 (* 0.70 *)
 | AA.PRO, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.SER, C   -> [], 1.61 (* 0.77 *)
 | AA.SER, CA  -> [], 1.88 (* 0.77 *)
 | AA.SER, CB  -> [], 1.88 (* 0.77 *)
 | AA.SER, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.SER, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.SER, OG  -> [HbD; HbA], 1.46 (* 0.66 *)
 | AA.THR, C   -> [], 1.61 (* 0.77 *)
 | AA.THR, CA  -> [], 1.88 (* 0.77 *)
 | AA.THR, CB  -> [], 1.88 (* 0.77 *)
 | AA.THR, CG2 -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.THR, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.THR, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.THR, OG1 -> [HbD; HbA], 1.46 (* 0.66 *)
 | AA.TRP, C   -> [], 1.61 (* 0.77 *)
 | AA.TRP, CA  -> [], 1.88 (* 0.77 *)
 | AA.TRP, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.TRP, CD1 -> [Aromatic], 1.76 (* 0.77 *)
 | AA.TRP, CD2 -> [Hydrophobic; Aromatic], 1.61 (* 0.77 *)
 | AA.TRP, CE2 -> [Aromatic], 1.61 (* 0.77 *)
 | AA.TRP, CE3 -> [Hydrophobic; Aromatic], 1.76 (* 0.77 *)
 | AA.TRP, CG  -> [Hydrophobic; Aromatic], 1.61 (* 0.77 *)
 | AA.TRP, CH2 -> [Hydrophobic; Aromatic], 1.76 (* 0.77 *)
 | AA.TRP, CZ2 -> [Hydrophobic; Aromatic], 1.76 (* 0.77 *)
 | AA.TRP, CZ3 -> [Hydrophobic; Aromatic], 1.76 (* 0.77 *)
 | AA.TRP, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.TRP, NE1 -> [Aromatic; HbD], 1.64 (* 0.70 *)
 | AA.TRP, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.TYR, C   -> [], 1.61 (* 0.77 *)
 | AA.TYR, CA  -> [], 1.88 (* 0.77 *)
 | AA.TYR, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.TYR, CD1 -> [Hydrophobic; Aromatic], 1.76 (* 0.77 *)
 | AA.TYR, CD2 -> [Hydrophobic; Aromatic], 1.76 (* 0.77 *)
 | AA.TYR, CE1 -> [Hydrophobic; Aromatic], 1.76 (* 0.77 *)
 | AA.TYR, CE2 -> [Hydrophobic; Aromatic], 1.76 (* 0.77 *)
 | AA.TYR, CG  -> [Hydrophobic; Aromatic], 1.61 (* 0.77 *)
 | AA.TYR, CZ  -> [Aromatic], 1.61 (* 0.77 *)
 | AA.TYR, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.TYR, O   -> [HbA], 1.42 (* 0.66 *)
 | AA.TYR, OH  -> [HbD], 1.46 (* 0.66 *)
 | AA.VAL, C   -> [], 1.61 (* 0.77 *)
 | AA.VAL, CA  -> [], 1.88 (* 0.77 *)
 | AA.VAL, CB  -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.VAL, CG1 -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.VAL, CG2 -> [Hydrophobic], 1.88 (* 0.77 *)
 | AA.VAL, N   -> [HbD], 1.64 (* 0.70 *)
 | AA.VAL, O   -> [HbA], 1.42 (* 0.66 *)
 | aa, atom    -> raise (Unknown_aa_atom_pair (aa, atom))
