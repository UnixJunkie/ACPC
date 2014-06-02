
(* to factorize code *)
let list_parmap ncores chunksize f l =
  Parmap.parmap ~ncores ~chunksize f (Parmap.L l)

(* to factorize code *)
let list_pariter ncores chunksize f l =
  Parmap.pariter ~ncores ~chunksize f (Parmap.L l)
