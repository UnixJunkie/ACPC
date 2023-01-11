
open Batteries

type aa_three =
  | ALA
  | ARG
  | ASN
  | ASP
  | CYS
  | GLU
  | GLN
  | GLY
  | HID
  | HIP
  | HIE
  | HIS
  | ILE
  | LEU
  | LYS
  | MET
  | PHE
  | PRO
  | SER
  | THR
  | TRP
  | TYR
  | VAL

type aa_one =
  | A
  | R
  | N
  | D
  | C
  | E
  | Q
  | G
  | H
  | I
  | L
  | K
  | M
  | F
  | P
  | S
  | T
  | W
  | Y
  | V

let three_to_one = function
  | ALA -> A
  | ARG -> R
  | ASN -> N
  | ASP -> D
  | CYS -> C
  | GLU -> E
  | GLN -> Q
  | GLY -> G
  | HIS -> H
  | HID -> H
  | HIP -> H
  | HIE -> H
  | ILE -> I
  | LEU -> L
  | LYS -> K
  | MET -> M
  | PHE -> F
  | PRO -> P
  | SER -> S
  | THR -> T
  | TRP -> W
  | TYR -> Y
  | VAL -> V

let one_to_three = function
  | A -> ALA
  | R -> ARG
  | N -> ASN
  | D -> ASP
  | C -> CYS
  | E -> GLU
  | Q -> GLN
  | G -> GLY
  | H -> HIS
  | I -> ILE
  | L -> LEU
  | K -> LYS
  | M -> MET
  | F -> PHE
  | P -> PRO
  | S -> SER
  | T -> THR
  | W -> TRP
  | Y -> TYR
  | V -> VAL

exception Unknown_aa_three of string

let aa_three_of_string = function
  | "ALA" -> ALA
  | "ARG" -> ARG
  | "ASN" -> ASN
  | "ASP" -> ASP
  | "CYS" -> CYS
  | "GLU" -> GLU
  | "GLN" -> GLN
  | "GLY" -> GLY
  | "HIS" -> HIS
  | "HID" -> HID
  | "HIP" -> HIP
  | "HIE" -> HIE
  | "ILE" -> ILE
  | "LEU" -> LEU
  | "LYS" -> LYS
  | "MET" -> MET
  | "PHE" -> PHE
  | "PRO" -> PRO
  | "SER" -> SER
  | "THR" -> THR
  | "TRP" -> TRP
  | "TYR" -> TYR
  | "VAL" -> VAL
  | "A" -> ALA
  | "R" -> ARG
  | "N" -> ASN
  | "D" -> ASP
  | "C" -> CYS
  | "E" -> GLU
  | "Q" -> GLN
  | "G" -> GLY
  | "H" -> HIS
  | "I" -> ILE
  | "L" -> LEU
  | "K" -> LYS
  | "M" -> MET
  | "F" -> PHE
  | "P" -> PRO
  | "S" -> SER
  | "T" -> THR
  | "W" -> TRP
  | "Y" -> TYR
  | "V" -> VAL
  | unk   -> raise (Unknown_aa_three unk)

exception Unknown_aa_one of string

let aa_one_of_string = function
  | "ALA" -> A
  | "ARG" -> R
  | "ASN" -> N
  | "ASP" -> D
  | "CYS" -> C
  | "GLU" -> E
  | "GLN" -> Q
  | "GLY" -> G
  | "HIS" -> H
  | "HID" -> H
  | "HIP" -> H
  | "HIE" -> H
  | "ILE" -> I
  | "LEU" -> L
  | "LYS" -> K
  | "MET" -> M
  | "PHE" -> F
  | "PRO" -> P
  | "SER" -> S
  | "THR" -> T
  | "TRP" -> W
  | "TYR" -> Y
  | "VAL" -> V
  | "A" -> A
  | "R" -> R
  | "N" -> N
  | "D" -> D
  | "C" -> C
  | "E" -> E
  | "Q" -> Q
  | "G" -> G
  | "H" -> H
  | "I" -> I
  | "L" -> L
  | "K" -> K
  | "M" -> M
  | "F" -> F
  | "P" -> P
  | "S" -> S
  | "T" -> T
  | "W" -> W
  | "Y" -> Y
  | "V" -> V
  | unk -> raise (Unknown_aa_one unk)

let string_of_aa_three = function
  | ALA -> "ALA"
  | ARG -> "ARG"
  | ASN -> "ASN"
  | ASP -> "ASP"
  | CYS -> "CYS"
  | GLU -> "GLU"
  | GLN -> "GLN"
  | GLY -> "GLY"
  | HIS -> "HIS"
  | HID -> "HID"
  | HIP -> "HIP"
  | HIE -> "HIE"
  | ILE -> "ILE"
  | LEU -> "LEU"
  | LYS -> "LYS"
  | MET -> "MET"
  | PHE -> "PHE"
  | PRO -> "PRO"
  | SER -> "SER"
  | THR -> "THR"
  | TRP -> "TRP"
  | TYR -> "TYR"
  | VAL -> "VAL"

let string_of_aa_one = function
  | A -> "A"
  | R -> "R"
  | N -> "N"
  | D -> "D"
  | C -> "C"
  | E -> "E"
  | Q -> "Q"
  | G -> "G"
  | H -> "H"
  | I -> "I"
  | L -> "L"
  | K -> "K"
  | M -> "M"
  | F -> "F"
  | P -> "P"
  | S -> "S"
  | T -> "T"
  | W -> "W"
  | Y -> "Y"
  | V -> "V"

let char_of_aa_one = function
  | A -> 'A'
  | R -> 'R'
  | N -> 'N'
  | D -> 'D'
  | C -> 'C'
  | E -> 'E'
  | Q -> 'Q'
  | G -> 'G'
  | H -> 'H'
  | I -> 'I'
  | L -> 'L'
  | K -> 'K'
  | M -> 'M'
  | F -> 'F'
  | P -> 'P'
  | S -> 'S'
  | T -> 'T'
  | W -> 'W'
  | Y -> 'Y'
  | V -> 'V'

type hydropathy = Hydrophobic | Hydrophilic | Neutral

type charge = Positive | Negative | Uncharged

type polarity = Polar | Nonpolar

(*
  data comes from:
  http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/\
    IMGTclasses.html#table

  paper:
  Pommié, C. et al., J. Mol. Recognit., 17, 17-32 (2004) PMID: 14872534,
  LIGM:284
*)

let classify = function
 | ALA -> (Hydrophobic, Uncharged, Nonpolar)
 | ARG -> (Hydrophilic, Positive , Polar   )
 | ASN -> (Hydrophilic, Uncharged, Polar   )
 | ASP -> (Hydrophilic, Negative , Polar   )
 | CYS -> (Hydrophobic, Uncharged, Nonpolar)
 | GLN -> (Hydrophilic, Uncharged, Polar   )
 | GLU -> (Hydrophilic, Negative , Polar   )
 | GLY -> (Neutral    , Uncharged, Nonpolar)
 | HIS -> (Neutral    , Positive , Polar   )
 | HID -> failwith "HID not implemented yet"
 | HIP -> failwith "HIP not implemented yet"
 | HIE -> failwith "HIE not implemented yet"
 | ILE -> (Hydrophobic, Uncharged, Nonpolar)
 | LEU -> (Hydrophobic, Uncharged, Nonpolar)
 | LYS -> (Hydrophilic, Positive , Polar   )
 | MET -> (Hydrophobic, Uncharged, Nonpolar)
 | PHE -> (Hydrophobic, Uncharged, Nonpolar)
 | PRO -> (Neutral    , Uncharged, Nonpolar)
 | SER -> (Neutral    , Uncharged, Polar   )
 | THR -> (Neutral    , Uncharged, Polar   )
 | TRP -> (Hydrophobic, Uncharged, Nonpolar)
 | TYR -> (Neutral    , Uncharged, Polar   )
 | VAL -> (Hydrophobic, Uncharged, Nonpolar)

let string_of_classes (h, c, p) =
  let h_str = match h with
    | Hydrophobic -> "Hydrophobic"
    | Hydrophilic -> "Hydrophilic"
    | Neutral     -> "Neutral"
  in
  let c_str = match c with
    | Positive  -> "Positive"
    | Negative  -> "Negative"
    | Uncharged -> "Uncharged"
  in
  let p_str = match p with
    | Polar    -> "Polar"
    | Nonpolar -> "Nonpolar"
  in
  (h_str ^ " " ^ c_str ^ " " ^ p_str)
