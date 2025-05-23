#!/usr/bin/python
import socket
import sys
import os
import subprocess
import argparse

###################################################################################################################################################################
# Important: The variables below are set to values that will make the run complete faster (as a tutorial example), but will not give scientifically valid results.
#            Please change them to the "normal" default values before a real run.
###################################################################################################################################################################

rosetta_scripts_path = os.path.expanduser("/home/ec2-user/local/rosetta.source.release-371/main/source/bin/rosetta_scripts.mpi.linuxgccrelease")
nstruct = 35 # Normally 35
max_minimization_iter = 5000 # Normally 5000
abs_score_convergence_thresh = 1.0 # Normally 1.0
number_backrub_trials = 35000 # Normally 35000
backrub_trajectory_stride = 35000 # Can be whatever you want, if you would like to see results from earlier time points in the backrub trajectory. 7000 is a reasonable number, to give you three checkpoints for a 35000 step run, but you could also set it to 35000 for quickest run time (as the final minimization and packing steps will only need to be run one time).
path_to_script = 'ddG-backrub.xml'

if not os.path.isfile(rosetta_scripts_path):
    print('ERROR: "rosetta_scripts_path" variable must be set to the location of the "rosetta_scripts" binary executable')
    print('This file might look something like: "rosetta_scripts.linuxgccrelease"')
    raise Exception('Rosetta scripts missing')

def run_flex_ddg( name, input_path, input_pdb_path, chains_to_move, nstruct_i ):
    output_directory = os.path.join( 'output', os.path.join( name, '%02d' % nstruct_i ) )
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    flex_ddg_args = [
        "mpirun",
        "--use-hwthread-cpus",
        "-n",
        "16",
        os.path.abspath(rosetta_scripts_path),
        "-s", f"{os.path.abspath(input_pdb_path)}",
        '-parser:protocol', os.path.abspath(path_to_script),
        '-parser:script_vars',
        'chainstomove=' + chains_to_move,
        'mutate_resfile_relpath=' + os.path.abspath( os.path.join( input_path, 'single_mut.resfile' ) ),
        'number_backrub_trials=%d' % number_backrub_trials,
        'max_minimization_iter=%d' % max_minimization_iter,
        'abs_score_convergence_thresh=%.1f' % abs_score_convergence_thresh,
        'backrub_trajectory_stride=%d' % backrub_trajectory_stride ,
        '-restore_talaris_behavior',
        '-in:file:fullatom',
        '-ignore_unrecognized_res',
        '-ignore_zero_occupancy false',
        '-ex1',
        '-ex2',
    ]

    log_path = os.path.join(output_directory, 'rosetta.out')

    print( 'Running Rosetta with args:' )
    print( ' '.join(flex_ddg_args) )
    print( 'Output logged to:', os.path.abspath(log_path) )
    print()

    outfile = open(log_path, 'w')
    process = subprocess.Popen(flex_ddg_args, stdout=outfile, stderr=subprocess.STDOUT, close_fds = True, cwd = output_directory)
    returncode = process.wait()
    outfile.close()

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser().parse_args()
    parser.add_argument('--job_name', type=str, required=True)
    parser.add_argument('--input_pdb_path', type=str, required=True)
    parser.add_argument('--chains_to_move', type=str, required=True, help='Chains to mutate, if mutations are one defferent chains, specify chain identifiers for each corresponding mutation and seperate by +, e.g. I+L')
    parser.add_argument('--mutation', type=str, required=True, help='Index and resname to mutate, connect multiple mutations with +, e.g. G12A+R13S')
    parser.add_argument('--output_path', type=str, required=False, default=None, help='Path to output directory, default is output under the input pdb directory')
    parser.add_argument('--Ncpu', '-N', type=int, required=False, default=16)
    args = parser.parse_args()
    
    cases = []
    for nstruct_i in range(1, nstruct + 1 ):
        # for case_name in os.listdir('inputs'):
        #     case_path = os.path.join( 'inputs', case_name )
        #     for f in os.listdir(case_path):
        #         if f.endswith('.pdb'):
        #             input_pdb_path = os.path.join( case_path, f )
        #             break

        chains_to_move = "I"
        case_name = "IL22-IL22RA"
        case_path = "."
        input_pdb_path = "/home/ec2-user/Ab_design/Interface/3dlq_IR_fix.pdb"

        cases.append( (args.job_name, case_path, input_pdb_path, chains_to_move, nstruct_i) )

    for args in cases:
        run_flex_ddg( *args )