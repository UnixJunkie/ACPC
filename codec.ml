
(* encoder (pre-parsing) and decoder of {mol2|pqr|pl}.bin files 
   .bin files are faster to read from than ASCII-based formats such as MOL2, PQR and PL *)

open Batteries
open Printf

module AC   = Autocorr
module MU   = My_utils
module Mol  = Molecule

let main () =

  let start = Unix.gettimeofday() in
  Log.set_log_level Log.INFO;
  Log.color_on();

  (* default options *)
  let input_fn = ref "" in
  let output_fn = ref "" in
  let cmd_line = Arg.align [
    "-i", Arg.Set_string input_fn, "in.{mol2|pqr|pl}[.bin] input file";
    "-o", Arg.Set_string output_fn, "out_fn output file"
  ] in
  Arg.parse cmd_line ignore
    (sprintf "Examples:\n\
              encoding: %s -i database.mol2 -o database.mol2.bin\n\            
              decoding: %s -i database.{mol2|pqr|pl}.bin"
       Sys.argv.(0) Sys.argv.(0));
  if !input_fn = "" then (
    Log.fatal (lazy ("-i is mandatory"));
    exit 1
  );
  let in_counter = ref 0 in
  let read_one_molecule = AC.get_reader !input_fn in
  if !output_fn = "" then (
    (* test reading speed *)
    MU.with_in_file !input_fn (fun input ->
      try
        while true do
          ignore(read_one_molecule in_counter input);
        done
      with End_of_file -> ()
    );
  ) else (
    (* encode to .bin file *)
    MU.enforce_any_file_extension !output_fn ["bin"];
    let out_counter = ref 0 in
    MU.with_out_file !output_fn (fun output ->
      MU.with_in_file !input_fn (fun input ->
        try
          while true do
            let read_in = read_one_molecule in_counter input in
            Mol.to_chan out_counter output read_in;
          done
        with End_of_file -> ()
      );
    );
    assert(!in_counter = !out_counter);
  );
  let stop = Unix.gettimeofday() in
  let elapsed = stop -. start in
  Log.infof "molecules: %d reading speed: %.2f molecules/s"
    !in_counter ((float_of_int !in_counter) /. elapsed)
;;

main()
