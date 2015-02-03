
open Batteries

open Molecule
open Legacy.Printf

module A    = Array
module L    = List
module MU   = My_utils
module Mol2 = Mol2_parser

(* imperative (i.e. modifiable) concrete unlabeled undirected graph *)
module G = struct

  module Int = struct
    type t = int
    let compare = Pervasives.compare
    let hash = Hashtbl.hash
    let equal = (=)
    let default = 0
  end

  include Graph.Imperative.Graph.Concrete(Int)

end

(* get the list of rotatable bonds of the molecule in the given MOL2 file *)
let get_rot_bonds fn =
  MU.run_command ("findrotatable.pl " ^ fn ^ " 1> /dev/null 2> /dev/null");
  let aux_fn = MU.filename_with_different_extension fn "mol2" "aux" in
  let aux_lines = MU.string_list_of_file aux_fn in
  match aux_lines with
  | _nb_rot_bonds :: bonds ->
    L.map
      (fun line ->
        try
          Scanf.sscanf line "ROT_BOND %d %d" (fun i j -> (i, j))
        with exn ->
          Log.error "in file %s, could not parse line %s" aux_fn line;
          raise exn
      )
      bonds
  | _ ->
    Log.fatal "could not read rotatable bonds from %s" aux_fn;
    exit 1

(* returns (molecule, molecular_graph) of the molecule in the given MOL2 file
   whose rotatable bonds have been removed *)
let process fn =
  MU.enforce_any_file_extension fn ["mol2"];
  let mol, all_bonds =
    MU.with_in_file fn
      (fun input -> Mol2.read_one_molecule_with_bonds (ref 0) input)
  in
  let rot_bonds = get_rot_bonds fn in
  (* extremities of the rotatable bonds; i.e. atom pairs *)
  let rot_bond_atoms =
    L.map
      (* MOL2 index to atom array index *)
      (fun (i, j) -> L.at mol.atoms (i - 1), L.at mol.atoms (j - 1))
      rot_bonds
  in
  let graph = G.create () in
  (* all atoms are vertices of the graph *)
  L.iteri
    (* in a MOL2 file, atom indexes start at 1; hence the (i + 1) *)
    (fun i _atom -> G.add_vertex graph (i + 1))
    mol.atoms;
  (* all bonds are edges of the graph *)
  L.iter
    (fun (_bond_i, bond_start, bond_end) -> G.add_edge graph bond_start bond_end)
    all_bonds;
  (* remove edges which are rotatable bonds *)
  L.iter
    (fun (bond_start, bond_end) -> G.remove_edge graph bond_start bond_end)
    rot_bonds;
  (mol, graph, rot_bond_atoms)

(* to find the strongly connected components of a graph *)
module Scc = Graph.Components.Make(G)

(* cut the molecule in MOL2 file 'fn' into its conformer-invariant
   sub parts (fragments)
   returns: a list of molecules (each one is a fragment) *)
let cut_rot_bonds fn =
  let molecule, graph, rot_bond_atoms = process fn in
  let atoms_array = A.of_list molecule.atoms in
  let strongly_connected_components = Scc.scc_list graph in
  (rot_bond_atoms,
   L.mapi
     (fun fragment_index atom_indexes ->
        let fragment_atoms = 
          L.map
           (* indexes start at 1 in the MOL2 file but not in atoms_array *)
            (fun i -> atoms_array.(i - 1))
            atom_indexes
        in
        Molecule.create
          (sprintf "%s_%d" molecule.name fragment_index)
          molecule.index (* keep same index as the input molecule *)
          fragment_atoms
     )
     strongly_connected_components)
