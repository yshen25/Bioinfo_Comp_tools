#! usr/env/bin python3
import argparse
import numpy as np
from Bio import PDB
from scipy.spatial import cKDTree
import string

def load_structure(pdb_file, structure_id):
    """Load a PDB file and return the structure."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(structure_id, pdb_file)
    return structure

def get_contact_residues(structure, contact_residues=None):
    """Get specified contact residues or all residues if None specified."""
    if contact_residues is None:
        residues = list(structure.get_residues())
    else:
        # Convert contact_residues to a set for faster lookup
        contact_set = set(contact_residues)
        residues = [res for res in structure.get_residues() if res.id[1] in contact_set]
    if not residues:
        raise ValueError(f"No residues found in the structure: {structure.id}")
    return residues

def get_center_of_mass(residues):
    """Calculate the center of mass of given residues (excluding hydrogens)."""
    coord_sum = np.zeros(3)
    total_atoms = 0
    for res in residues:
        for atom in res:
            if atom.element != 'H':
                coord_sum += atom.coord
                total_atoms += 1  # Assuming all heavy atoms have equal mass
    if total_atoms == 0:
        raise ValueError("No heavy atoms found in the provided residues for center of mass calculation.")
    return coord_sum / total_atoms

def get_heavy_atom_coords(residues):
    """Get the coordinates of all heavy atoms in the given residues."""
    coords = np.array([atom.coord for res in residues for atom in res if atom.element != 'H'])
    if len(coords) == 0:
        raise ValueError("No heavy atoms found in the provided residues.")
    return coords

def move_structure(structure, translation_vector):
    """Apply a translation vector to all atoms in the structure."""
    for atom in structure.get_atoms():
        atom.coord += translation_vector

def align_structures(antibody, antigen, antibody_residues, antigen_residues):
    """Align the centers of mass of the antibody and antigen along the x-axis."""
    # Calculate centers of mass based on contact residues
    antibody_com = get_center_of_mass(antibody_residues)
    antigen_com = get_center_of_mass(antigen_residues)

    # Align antigen COM to origin
    antigen_translation = -antigen_com
    move_structure(antigen, antigen_translation)

    # Move antibody far away along the x-axis and align COM
    antibody_translation = -antibody_com
    move_structure(antibody, antibody_translation)
    antibody_far_distance = 100.0  # Move antibody 100 Å along x-axis
    move_structure(antibody, np.array([antibody_far_distance, 0, 0]))

def adjust_distance(antibody, antigen, antibody_residues, antigen_residues, desired_distance):
    """Move the antibody closer to the antigen to achieve the desired minimal distance between contact residues."""
    # Get heavy atom coordinates from contact residues
    antibody_coords = get_heavy_atom_coords(antibody_residues)
    antigen_coords = get_heavy_atom_coords(antigen_residues)

    # Calculate current minimal distance
    current_min_distance = calculate_min_distance(antibody_coords, antigen_coords)
    print(f"Current minimal heavy atom distance: {current_min_distance:.2f} Å")

    # Calculate how much to move the antibody along x-axis
    distance_difference = current_min_distance - desired_distance
    move_structure(antibody, np.array([-distance_difference, 0, 0]))

    # Recalculate minimal distance after moving antibody
    antibody_coords = get_heavy_atom_coords(antibody_residues)
    new_min_distance = calculate_min_distance(antibody_coords, antigen_coords)
    print(f"New minimal heavy atom distance: {new_min_distance:.2f} Å")

def calculate_min_distance(coords1, coords2):
    """Calculate the minimal distance between two sets of atom coordinates."""
    tree1 = cKDTree(coords1)
    distances, _ = tree1.query(coords2)
    min_distance = distances.min()
    return min_distance

def combine_structures(antibody, antigen):
    """Combine the antibody and antigen structures into a single complex."""
    # Create a new structure for the complex
    complex_structure = PDB.Structure.Structure("complex")
    complex_model = PDB.Model.Model(0)
    complex_structure.add(complex_model)

    # Add antigen to the complex
    for chain in antigen[0]:
        complex_model.add(chain.copy())

    # Collect chain IDs to ensure uniqueness
    antigen_chain_ids = set(chain.id for chain in antigen[0])

    # Add antibody to the complex, ensuring chain IDs are unique
    available_chain_ids = [cid for cid in string.ascii_letters + string.digits if cid not in antigen_chain_ids]
    for chain in antibody[0]:
        original_chain_id = chain.id
        if original_chain_id in antigen_chain_ids:
            if not available_chain_ids:
                raise ValueError("No available chain IDs to assign to antibody chain.")
            new_chain_id = available_chain_ids.pop(0)
            print(f"Assigning new chain ID '{new_chain_id}' to antibody chain '{original_chain_id}'.")
            chain.id = new_chain_id
            antigen_chain_ids.add(new_chain_id)
        complex_model.add(chain.copy())

    # Renumber atom serial numbers
    atom_number = 1
    for atom in complex_structure.get_atoms():
        atom.serial_number = atom_number
        atom_number += 1

    return complex_structure

def save_structure(structure, file_name):
    """Save a structure to a PDB file."""
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(file_name)
    print(f"Complex structure saved as '{file_name}'")

def generate_complex(antibody_file, antigen_file, complex_file, antibody_contacts=None, antigen_contacts=None, distance=4.0):
    """Generate the antibody-antigen complex with simplified algorithm and contact residues."""
    # Load structures
    antibody = load_structure(antibody_file, "antibody")
    antigen = load_structure(antigen_file, "antigen")

    # Get contact residues
    antibody_residues = get_contact_residues(antibody, antibody_contacts)
    antigen_residues = get_contact_residues(antigen, antigen_contacts)

    # Align structures
    align_structures(antibody, antigen, antibody_residues, antigen_residues)

    # Adjust distance
    adjust_distance(antibody, antigen, antibody_residues, antigen_residues, distance)

    # Combine structures
    complex_structure = combine_structures(antibody, antigen)

    # Save complex structure
    save_structure(complex_structure, complex_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a complex structure of antibody and antigen.")
    parser.add_argument("antibody", help="Antibody PDB file")
    parser.add_argument("antigen", help="Antigen PDB file")
    parser.add_argument("complex", help="Path of complex pdb output file")
    parser.add_argument("--antibody_contacts", nargs="+", type=int, help="Antibody contact residue numbers")
    parser.add_argument("--antigen_contacts", nargs="+", type=int, help="Antigen contact residue numbers")
    parser.add_argument("--distance", type=float, default=10, help="Distance between contact surfaces")

    args = parser.parse_args()

    generate_complex(args.antibody, args.antigen, args.complex, args.antibody_contacts, args.antigen_contacts, args.distance)