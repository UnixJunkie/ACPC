
(* reader of PQR (.pqr) files *)

open Batteries
open Legacy.Printf
open Molecule

module A   = Array
module At  = Atom
module L   = List
module Pqr = Pqr_parser

let main () =

  Log.set_log_level Log.INFO;
  Log.color_on();

  if A.length Sys.argv <> 2 then begin
    Log.fatal "incorrect number of parameters\n\
               usage: %s pqr_file" Sys.argv.(0);
    (*       0  1 *)
    exit 1;
  end else
    let molecules = Pqr.read_molecules Sys.argv.(1) in
    L.iter
      (fun mol ->
         printf "%s\n" mol.name;
         printf "%d\n" mol.index;
         L.iter
           (fun a -> printf "%s\n" (At.to_pqr_string a))
           mol.atoms
      )
      molecules
;;

main()
