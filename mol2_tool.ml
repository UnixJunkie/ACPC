
(* basic operations on the Tripos Mol2 File Format:
   - extract multiple molecules
   - correct the atom names to be unique if necessary
   Reference document: Tripos Mol2 File Format, SYBYL 7.1 (Mid-2005),
   http://www.tripos.com/data/support/mol2.pdf *)

(* This tool reads the whole MOL2 file in memory at once and hence
   will blow up on too big files *)

open Batteries
open Printf

module A    = Array
module Fn   = Filename
module L    = List
module MU   = My_utils
module Mol  = Molecule
module Mol2 = Mol2_parser
module S    = String

type parser_state = Reading_atoms | Other

let process_mol2_line state_r count_r l acc =
  if (S.starts_with l "@<TRIPOS>MOLECULE") then
    match acc with
        [] -> (None            , [l]) (* first molecule in the file *)
      | _  -> (Some (L.rev acc), [l])
  else begin
    let l' =
      if (S.starts_with l " ") && (!state_r = Reading_atoms) then
        match Str.split (Str.regexp "[ ]+") l with
            _id::name::x::_y::_z::_type::_optional_fields -> (
              (* atom_name fields from the original mol2 file are made unique
                 so that pdb2pqr will not choke on them later on *)
              let name_pos    = S.find l name                             in
              let x_pos       = S.find l x                                in
              let prfx        = String.sub l 0 name_pos                   in
              let sufx        = String.sub l x_pos ((S.length l) - x_pos) in
              let stripped    =
                Str.global_replace (Str.regexp "[0-9]*$") "" name         in
              let unique_name = stripped ^ (string_of_int !count_r)       in
              incr count_r;
              let prev_len    = x_pos - name_pos                          in
              let new_name    =
                unique_name ^
                  (S.make (prev_len - (S.length unique_name)) ' ')        in
              sprintf "%s%s %s" prfx new_name sufx)
          | _ -> failwith ("mol2_parser.ml: process_mol2_line: \
                            invalid atom line: " ^ l)
      else l
    in
    if (S.starts_with l "@<TRIPOS>ATOM") then
      (state_r := Reading_atoms; count_r := 0);
    if (S.starts_with l "@<TRIPOS>BOND") then
      state_r := Other;
    (None, l' :: acc)
  end

let molecule_name m = match m with
    []  -> failwith "mol2_parser.ml: molecule_name: empty molecule"
  | [_] -> failwith "mol2_parser.ml: molecule_name: molecule without a name"
  | _ :: name :: _ -> name

(* remove all lines before the first '@<TRIPOS>MOLECULE' one *)
let rec skip_header mol2_lines =
  match mol2_lines with
      []      -> []
    | x :: xs ->
      if S.starts_with x "@<TRIPOS>MOLECULE"
      then mol2_lines
      else skip_header xs

(* get the list of molecules *)
let read_mol2_file f =
  MU.enforce_any_file_extension f [".mol2"];
  let state = ref Other                              in
  let count = ref 0                                  in
  let lines = skip_header (MU.string_list_of_file f) in
  let molecules, molecule =
    L.fold_left
      (fun (global_acc, local_acc) l ->
        match process_mol2_line state count l local_acc with
          | Some molecule, new_local_acc ->
            (molecule :: global_acc, new_local_acc)
          | None         , new_local_acc ->
            (global_acc            , new_local_acc))
      ([], []) lines in
  let all_molecules = (L.rev molecule) :: molecules in
  L.rev_map
    (* the "@<TRIPOS>SUBSTRUCTURE" line is required at the end of
       each molecule by pdb2pqr.py *)
    (fun m ->
      let corrected_m =
        if L.mem "@<TRIPOS>SUBSTRUCTURE" m then
          m
        else
          L.append m ["@<TRIPOS>SUBSTRUCTURE"] in
      (corrected_m, molecule_name m))
    all_molecules

(* explode a multiple molecules mol2 file: each molecule will be stored under
   its own directory in a single molecule mol2 file
   - return the list of files created *)
let explode input_f =
  let molecules = read_mol2_file    input_f in
  let out_dir   = Fn.chop_extension input_f in
  Log.info (lazy (sprintf "molecules in %s: %d" input_f (L.length molecules)));
  L.mapi
    (fun j (m, m_name) ->
      let new_dir  = sprintf "%s_%07d" out_dir (j + 1) in
      (* a directory needs at least rwx for the user so that we can
         go and write into it later on *)
      (try Unix.mkdir new_dir 0o700 with exn ->
        (match exn with
          | Unix.Unix_error (err_code, f_name, param) ->
            (match err_code with
              | Unix.EEXIST -> ()
            (* dir already exist, we don't care so much in fact *)
            (*
              printf
              "mol2_tool.ml: warning: directory already exist: %s\n"
              new_dir
            *)
              | _ -> failwith (sprintf
                                 "mol2_tool.ml: %s: %s %s"
                                 (Unix.error_message err_code) f_name param))
          | _ -> raise exn));
      let new_file = sprintf "%s/%07d.mol2" new_dir (j + 1) in
      MU.string_list_to_file new_file m;
      (new_file, m_name)
    )
    molecules

let main () =

  Log.set_log_level Log.INFO;
  Log.color_on();

  (* default options *)
  let input_fn = ref "" in
  let cmd_line = Arg.align [
    "-i", Arg.Set_string input_fn, "in.mol2 input file";
  ] in
  Arg.parse cmd_line ignore
    (sprintf "Example:\n%s -i database.mol2" Sys.argv.(0));
  if !input_fn = "" then (
    Log.fatal (lazy ("-i some_file.[mol2|bin] is mandatory"));
    exit 1
  );
  MU.enforce_any_file_extension !input_fn ["mol2"];
  (* "explode" the file *)
  let new_files = explode !input_fn in
  Log.info (lazy "Files created:");
  L.iter
    (fun (fn, _m_name) -> Log.info (lazy fn))
    new_files
;;

main()
