
(* reader of PL (.pl) files *)

open Batteries
open Printf
open Molecule

module A  = Array
module At = Atom
module L  = List
module Pl = Pl_parser

let main () =

  Log.set_log_level Log.INFO;
  Log.color_on();

  if A.length Sys.argv <> 2 then begin
    Log.fatal "incorrect number of parameters\n\
               usage: %s pl_file" Sys.argv.(0);
    (*       0  1 *)
    exit 1;
  end else
    let molecules = Pl.read_molecules Sys.argv.(1) in
    L.iter
      (fun mol ->
         printf "%s\n" mol.name;
         printf "%d\n" mol.index;
         L.iter
           (fun a -> printf "%s\n" (At.to_pl_string a))
           mol.atoms
      )
      molecules
;;

main()
