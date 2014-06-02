
(* module to store the results found in several .scores files *)

type t = { mol2_score : float ;
           pqr_score  : float ;
           pl_score   : float ;
           nb_scores  : int   }

let create () =
  { mol2_score = 0. ;
    pqr_score  = 0. ;
    pl_score   = 0. ;
    nb_scores  = 0  }

let set_mol2_score r s =
  { mol2_score = s               ;
    pqr_score  = r.pqr_score     ;
    pl_score   = r.pl_score      ;
    nb_scores  = 1 + r.nb_scores }

let set_pqr_score r s =
  { mol2_score = r.mol2_score    ;
    pqr_score  = s               ;
    pl_score   = r.pl_score      ;
    nb_scores  = 1 + r.nb_scores }

let set_pl_score r s =
  { mol2_score = r.mol2_score    ;
    pqr_score  = r.pqr_score     ;
    pl_score   = s               ;
    nb_scores  = 1 + r.nb_scores }

let rescore r w_mol2 w_pqr w_pl =
  assert(r.nb_scores = 3);
  (r.mol2_score *. w_mol2 +.
   r.pqr_score  *. w_pqr  +.
   r.pl_score   *. w_pl)

let to_triplet r =
  (r.mol2_score,
   r.pqr_score,
   r.pl_score)
