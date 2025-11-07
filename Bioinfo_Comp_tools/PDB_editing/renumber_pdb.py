from Bio.PDB import PDBParser, PDBIO, Select
from typing import Dict, List, Tuple, Optional
import warnings
import os
from pathlib import Path


class ChainSelect(Select):
    """Select specific chains for output"""
    def __init__(self, chain_ids):
        self.chain_ids = set(chain_ids) if chain_ids else None
    
    def accept_chain(self, chain):
        if self.chain_ids is None:
            return True
        return chain.id in self.chain_ids


class RenumberPDB:
    """
    A class for renumbering residues in PDB files.
    
    Supports flexible renumbering schemes including:
    - Start-based renumbering
    - End-based renumbering
    - Full range renumbering
    - Gap preservation or removal
    """
    
    def __init__(self, in_pdb: str):
        """
        Initialize with a PDB file.
        
        Args:
            in_pdb: Path to input PDB file
        """
        parser = PDBParser(PERMISSIVE=1, QUIET=True)
        self.structure = parser.get_structure("input", in_pdb)
        self.chain_ids = {chain.id for chain in self.structure.get_chains()}
        self.renumber_dict = {}
        self.input_pdb_path = in_pdb
        
    def get_chain_info(self, chain_id: str) -> Dict:
        """
        Get information about a specific chain.
        
        Args:
            chain_id: Chain identifier
            
        Returns:
            Dictionary with chain information
        """
        if chain_id not in self.chain_ids:
            raise ValueError(f'Chain {chain_id} not in structure. Available: {self.chain_ids}')
        
        chain = None
        for c in self.structure.get_chains():
            if c.id == chain_id:
                chain = c
                break
        
        residues = list(chain.get_residues())
        res_ids = [res.id[1] for res in residues]
        
        return {
            'chain_id': chain_id,
            'num_residues': len(residues),
            'first_resid': res_ids[0],
            'last_resid': res_ids[-1],
            'residue_ids': res_ids,
            'has_gaps': any(res_ids[i+1] - res_ids[i] > 1 for i in range(len(res_ids)-1))
        }
    
    def gen_replace_dict(self, chain_id: str, start_index: Optional[int] = None, 
                        end_index: Optional[int] = None, keep_gap: bool = True) -> Dict[int, int]:
        """
        Generate a dictionary mapping old residue indices to new ones.
        
        Args:
            chain_id: Chain to renumber
            start_index: New starting index
            end_index: New ending index
            keep_gap: If True, preserve gaps in numbering; if False, sequential numbering
            
        Returns:
            Dictionary mapping old_resid -> new_resid
        """
        info = self.get_chain_info(chain_id)
        old_ids = info['residue_ids']
        num_residues = len(old_ids)
        
        # Validate inputs
        if start_index is None and end_index is None:
            raise ValueError('Must provide at least one of start_index or end_index')
        
        replace_dict = {}
        
        # Mode 1: Only start_index given
        if start_index is not None and end_index is None:
            offset = start_index - old_ids[0]
            if keep_gap:
                replace_dict = {old_id: old_id + offset for old_id in old_ids}
            else:
                replace_dict = {old_id: start_index + i for i, old_id in enumerate(old_ids)}
        
        # Mode 2: Only end_index given
        elif start_index is None and end_index is not None:
            offset = end_index - old_ids[-1]
            if keep_gap:
                replace_dict = {old_id: old_id + offset for old_id in old_ids}
            else:
                replace_dict = {old_id: end_index - num_residues + i + 1 
                               for i, old_id in enumerate(old_ids)}
        
        # Mode 3: Both indices given
        else:
            expected_length = end_index - start_index + 1
            if not keep_gap and expected_length != num_residues:
                warnings.warn(
                    f'Range length ({expected_length}) differs from number of residues '
                    f'({num_residues}). Using sequential numbering.'
                )
            
            if keep_gap:
                # Scale the existing gaps to fit the new range
                old_range = old_ids[-1] - old_ids[0]
                new_range = end_index - start_index
                if old_range == 0:
                    replace_dict = {old_id: start_index for old_id in old_ids}
                else:
                    scale = new_range / old_range
                    replace_dict = {
                        old_id: start_index + int((old_id - old_ids[0]) * scale)
                        for old_id in old_ids
                    }
            else:
                replace_dict = {old_id: start_index + i for i, old_id in enumerate(old_ids)}
        
        return replace_dict
    
    def renumber(self, chain_id: str, start_index: Optional[int] = None, 
                end_index: Optional[int] = None, keep_gap: bool = True):
        """
        Renumber one chain at a time.
        
        Args:
            chain_id: Chain to modify
            start_index: New start index
            end_index: New end index
            keep_gap: If True, preserve gaps; if False, sequential numbering
        """
        # Generate replacement dictionary
        replace_dict = self.gen_replace_dict(chain_id, start_index, end_index, keep_gap)
        
        # Store for later reference
        self.renumber_dict[chain_id] = replace_dict
        
        # Apply renumbering
        self._apply_renumbering(chain_id, replace_dict)
        
        return replace_dict
    
    def _apply_renumbering(self, chain_id: str, replace_dict: Dict[int, int]):
        """
        Apply the renumbering to the structure.
        
        Args:
            chain_id: Chain to renumber
            replace_dict: Mapping of old -> new residue IDs
        """
        # Step 1: Renumber to temporary negative values to avoid conflicts
        target_chain = None
        for chain in self.structure.get_chains():
            if chain.id == chain_id:
                target_chain = chain
                break
        
        if target_chain is None:
            return
        
        temp_mapping = {}
        for i, residue in enumerate(target_chain.get_residues()):
            old_id = residue.id[1]
            temp_id = -(i + 1)
            residue.id = (' ', temp_id, ' ')
            temp_mapping[temp_id] = old_id
        
        # Step 2: Apply final numbering
        for residue in target_chain.get_residues():
            temp_id = residue.id[1]
            old_id = temp_mapping[temp_id]
            new_id = replace_dict[old_id]
            residue.id = (' ', new_id, ' ')
    
    def save_pdb(self, out_pdb: Optional[str] = None, out_path: Optional[str] = None, 
                 chains: Optional[List[str]] = None):
        """
        Save the renumbered structure to a PDB file.
        
        Args:
            out_pdb: Output PDB file name (optional)
            out_path: Output directory path (optional)
            chains: List of chain IDs to save (None = save all)
            
        Behavior:
            - Only out_pdb: Save to out_pdb in current directory
            - Only out_path: Save to out_path with original name + '_renumber' suffix
            - Both: Save to out_path/out_pdb
            - Neither: Save to current directory with original name + '_renumber' suffix
        """
        # Determine output file path
        if out_pdb is None and out_path is None:
            # Neither given: use input name with _renumber suffix
            input_path = Path(self.input_pdb_path)
            out_file = input_path.stem + '_renumber' + input_path.suffix
            output_path = out_file
        elif out_pdb is not None and out_path is None:
            # Only out_pdb given: use as-is
            output_path = out_pdb
        elif out_pdb is None and out_path is not None:
            # Only out_path given: use input name with _renumber suffix in out_path
            input_path = Path(self.input_pdb_path)
            out_file = input_path.stem + '_renumber' + input_path.suffix
            output_path = os.path.join(out_path, out_file)
            # Create directory if it doesn't exist
            os.makedirs(out_path, exist_ok=True)
        else:
            # Both given: combine them
            output_path = os.path.join(out_path, out_pdb)
            # Create directory if it doesn't exist
            os.makedirs(out_path, exist_ok=True)
        
        # Save the structure
        io = PDBIO()
        io.set_structure(self.structure)
        
        if chains:
            # Validate chain IDs
            invalid_chains = set(chains) - self.chain_ids
            if invalid_chains:
                raise ValueError(f'Invalid chain IDs: {invalid_chains}')
            selector = ChainSelect(chains)
        else:
            selector = ChainSelect(None)
        
        io.save(output_path, selector)
        print(f"Saved renumbered structure to: {output_path}")
    
    def get_renumber_summary(self) -> str:
        """
        Get a summary of all renumbering operations performed.
        
        Returns:
            Formatted string summary
        """
        if not self.renumber_dict:
            return "No renumbering operations performed yet."
        
        summary = ["Renumbering Summary:", "=" * 50]
        for chain_id, mapping in self.renumber_dict.items():
            old_ids = sorted(mapping.keys())
            new_ids = [mapping[k] for k in old_ids]
            summary.append(f"\nChain {chain_id}:")
            summary.append(f"  Original range: {old_ids[0]} - {old_ids[-1]}")
            summary.append(f"  New range:      {new_ids[0]} - {new_ids[-1]}")
            summary.append(f"  Total residues: {len(mapping)}")
        
        return "\n".join(summary)


def parse_chain_index(chain_index_str: str) -> Dict[str, Tuple[Optional[int], Optional[int]]]:
    """
    Parse chain-index string into a dictionary of chain renumbering instructions.
    
    Format: A:1-30;B:5;C:-50;D:10-100
    Where:
        - A:1-30  means chain A from 1 to 30
        - B:5     means chain B starting from 5
        - C:-50   means chain C ending at 50
        - D:10-100 means chain D from 10 to 100
    
    Args:
        chain_index_str: Semicolon-separated chain specifications
        
    Returns:
        Dictionary mapping chain_id -> (start_index, end_index)
    """
    chain_specs = {}
    
    # Split by semicolon to get individual chain specifications
    for spec in chain_index_str.split(';'):
        spec = spec.strip()
        if not spec:
            continue
            
        # Split chain ID from index specification
        if ':' not in spec:
            raise ValueError(f"Invalid format '{spec}'. Expected format: CHAIN:INDEX or CHAIN:START-END")
        
        chain_id, index_part = spec.split(':', 1)
        chain_id = chain_id.strip()
        index_part = index_part.strip()
        
        # Parse index specification
        if '-' in index_part:
            # Range specification: START-END
            parts = index_part.split('-')
            
            # Handle cases like ":-50" (no start) or "5-" (no end)
            if len(parts) == 2:
                start_str, end_str = parts
                start_index = int(start_str) if start_str else None
                end_index = int(end_str) if end_str else None
            else:
                raise ValueError(f"Invalid range format in '{spec}'. Use START-END, -END, or START-")
        else:
            # Single number: start index only
            start_index = int(index_part)
            end_index = None
        
        chain_specs[chain_id] = (start_index, end_index)
    
    return chain_specs


def main():
    """CLI interface for PDB renumbering."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Renumber residues in PDB files with flexible multi-chain support.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Renumber chain A starting from 1
  %(prog)s input.pdb -c A:1
  
  # Renumber chain A from 1-100, chain B from 1-50
  %(prog)s input.pdb -c "A:1-100;B:1-50"
  
  # Renumber chain A ending at 200, chain B starting at 5
  %(prog)s input.pdb -c "A:-200;B:5"
  
  # Specify output file and path
  %(prog)s input.pdb -c A:1 -o output.pdb -p ./results
  
  # Keep gaps in numbering
  %(prog)s input.pdb -c A:1-100 --keep-gap
        """
    )
    
    parser.add_argument('input', 
                       help='Input PDB file path')
    parser.add_argument('-c', '--chain-index', 
                       required=True,
                       help='Chain renumbering specification (e.g., "A:1-30;B:5;C:-50")')
    parser.add_argument('-o', '--output', 
                       default=None,
                       help='Output PDB file name (optional)')
    parser.add_argument('-p', '--path', 
                       default=None,
                       help='Output directory path (optional)')
    parser.add_argument('--keep-gap', 
                       action='store_true',
                       help='Preserve gaps in residue numbering (default: False)')
    parser.add_argument('--save-chains',
                       default=None,
                       help='Comma-separated list of chains to save (default: save all)')
    parser.add_argument('-v', '--verbose',
                       action='store_true',
                       help='Print detailed information')
    
    args = parser.parse_args()
    
    try:
        # Parse chain-index specification
        chain_specs = parse_chain_index(args.chain_index)
        
        if args.verbose:
            print(f"Input file: {args.input}")
            print(f"Chain specifications: {chain_specs}")
            print(f"Keep gaps: {args.keep_gap}")
            print("-" * 50)
        
        # Initialize renumberer
        renumberer = RenumberPDB(args.input)
        
        if args.verbose:
            print(f"Available chains: {renumberer.chain_ids}\n")
        
        # Process each chain
        for chain_id, (start_idx, end_idx) in chain_specs.items():
            if chain_id not in renumberer.chain_ids:
                print(f"Warning: Chain '{chain_id}' not found in PDB. Skipping.")
                continue
            
            if args.verbose:
                info = renumberer.get_chain_info(chain_id)
                print(f"Chain {chain_id} - Original: {info['first_resid']} to {info['last_resid']} "
                      f"({info['num_residues']} residues)")
            
            # Renumber the chain
            renumberer.renumber(chain_id, 
                              start_index=start_idx, 
                              end_index=end_idx, 
                              keep_gap=args.keep_gap)
            
            if args.verbose:
                mapping = renumberer.renumber_dict[chain_id]
                old_ids = sorted(mapping.keys())
                new_ids = [mapping[k] for k in old_ids]
                print(f"Chain {chain_id} - New: {new_ids[0]} to {new_ids[-1]}\n")
        
        # Print summary
        if args.verbose or len(chain_specs) > 1:
            print(renumberer.get_renumber_summary())
            print()
        
        # Parse save_chains if provided
        save_chains = None
        if args.save_chains:
            save_chains = [c.strip() for c in args.save_chains.split(',')]
        
        # Save the renumbered structure
        renumberer.save_pdb(out_pdb=args.output, 
                          out_path=args.path, 
                          chains=save_chains)
        
        print("✓ Renumbering completed successfully!")
        
    except Exception as e:
        print(f"Error: {e}")
        import sys
        sys.exit(1)


# Example usage
if __name__ == "__main__":
    main()