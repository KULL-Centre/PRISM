#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# @file: compute_ddG.py
##
# @brief: 	 Compute ddGs of mutation
# @details: Use the Rosetta membrane framework to compute the ddG of unfolding of
# a membrane protein in Rosetta (uses packer, mutate.py from Evan Baugh)
##
# @author: Rebecca F. Alford (rfalford12@gmail.com)
# @author: JKLeman (julia.koehler1982@gmail.com)


# Adapted by Johanna KS Tiemann (johanna.tiemann@gmail.com)


# Tools
import sys
import os
import numpy as np
##import commands
##import random
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname(os.path.realpath(__file__))

# Rosetta-specific imports
from rosetta import *
# pyrosetta.init()

import pyrosetta.rosetta.protocols.membrane
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta import create_score_function
from pyrosetta.rosetta.core.pack.task import TaskFactory
from rosetta.utility import vector1_bool
from pyrosetta.rosetta.core.chemical import aa_from_oneletter_code
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pose import PDBInfo
from pyrosetta.rosetta.core.chemical import VariantType

###############################################################################

# @brief Main - Add Membrane to Pose, Compute ddG


"""
ToDo:
        # Create a membrane energy function enabled by pH mode
        # Includes two terms not standard in the smoothed energy function: pH energy
        # and fa_elec mpframework_pHmode_fa_2014
"""


def main(args):

    parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)

    # input options
    parser.add_option('--in_pdb', '-p',
                      action="store",
                      help="Input PDB file.", )

    parser.add_option('--in_span', '-s',
                      action="store",
                      help="Input spanfile.", )

    parser.add_option('--out', '-o',
                      action="store", default='ddG.out',
                      help="Output filename with pose residue numbering. Default: 'ddG.out'", )

    parser.add_option('--out_add',
                      action="store", default='ddG_additional.out',
                      help="Output filename of additional information with pose residue numbering. Default: 'ddG_additional.out'", )

    parser.add_option('--res', '-r',
                      action="store",
                      help="Pose residue number to mutate.", )

    parser.add_option('--mut', '-m',
                      action="store",
                      help="One-letter code of residue identity of the mutant. Example: A181F would be 'F'", )

    parser.add_option('--repack_radius', '-a',
                      action="store", default=8.0,
                      type=float,
                      help="Repack the residues within this radius",)

    parser.add_option('--output_breakdown', '-b',
                      action="store", default="scores.sc",
                      help="Output mutant and native score breakdown by weighted energy term into a scorefile", )

    parser.add_option('--include_pH', '-z',
                      action="store", default=0,
                      help="Include pH energy terms: pH_energy and fa_elec. Default false.", )

    parser.add_option('--pH_value', '-q',
                      action="store", default=7,
                      type=float,
                      help="Predict ddG and specified pH value. Default 7. Will not work if include pH is not passed", )

    parser.add_option('--repeats',
                      action="store", default=3,
                      type=int,
                      help="Number of ddG calculations. Default 3.", )

    parser.add_option('--lipids', '-l',
                      action="store", default='DLPC',
                      choices=['DLPC', ],
                      help="Specify lipids from options.", )

    parser.add_option('--temperature', '-t',
                      action="store", default=37.0,
                      type=float,
                      help="Experimental temperature.", )

    # parse options
    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    # Check the required inputs (PDB file, spanfile) are present
    if (not Options.in_pdb or not Options.in_span or not Options.res):
        sys.exit("Must provide flags '-in_pdb', '-in_span', and '-res'! Exiting...")

    if (float(Options.pH_value) < 0 or float(Options.pH_value) > 14):
        sys.exit("Specified pH value must be between 0-14: Exiting...")

    # Initialize Rosetta options from user options. Enable pH mode if
    # applicable
    rosetta_options = (f'-mp:setup:spanfiles {Options.in_span}'
                       f' -mp:lipids:temperature {Options.temperature}'
                       f' -mp:lipids:composition {Options.lipids}'
                       f' -mp:lipids:has_pore false'
                       f' -run:constant_seed'
                       f' -in:ignore_unrecognized_res'
                       f' -pH_mode -value_pH {Options.pH_value}')

    # Initialize Rosetta based on user inputs
    pyrosetta.init(extra_options=rosetta_options)



    # Create an energy function
    sfxn = pyrosetta.rosetta.core.scoring.ScoreFunction()
    # Create a smoothed membrane full atom energy function (pH 7 calculations)
    sfxn = create_score_function("franklin2019")


    # Load Pose, & turn on the membrane
    pose = pyrosetta.rosetta.core.import_pose.pose_from_file(Options.in_pdb)

    # Add Membrane to Pose
    add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( Options.in_span)
    add_memb.apply(pose)

    # Setup in a topology based membrane
    init_mem_pos = pyrosetta.rosetta.protocols.membrane.MembranePositionFromTopologyMover()
    init_mem_pos.apply(pose)




    # Repack the native rotamer and residues within the repack radius
    native_res = pose.residue(int(Options.res)).name1()
    repacked_native = mutate_residue(
        pose, int(Options.res), native_res, Options.repack_radius, sfxn)

    # to output score breakdown, start by printing the score labels in
    # the top of the file
    print_score_labels_to_file(repacked_native, sfxn, Options.output_breakdown)

    # Compute mutations
    if (Options.mut):
        AAs = list(Options.mut)
    else:
        AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    for aa in AAs:
        with open(Options.out, 'a') as f, open(Options.out_add, 'a') as fp2:
            result_ddG = []
            for rep in range(Options.repeats):
                ddGs = compute_ddG(repacked_native, sfxn, int(
                    Options.res), aa, Options.repack_radius, Options.output_breakdown)
                res_str = str(native_res) + str(Options.res) + str(ddGs[0])
                result_ddG.append(ddGs[3])
                fp2.write(res_str + "," + str(ddGs[3]) + "," +
                          str(ddGs[1]) + "," + str(ddGs[2]) + "\n")
            result_ddG = np.array(result_ddG)
            f.write(f'{res_str},{round(np.mean(result_ddG),3)},{round(np.std(result_ddG),3)}\n')

###############################################################################

# @brief Compute ddG of mutation in a protein at specified residue and AA position


def compute_ddG(pose, sfxn, resnum, aa, repack_radius, sc_file):

    # Score Native Pose
    native_score = sfxn(pose)

    # Perform Mutation at residue <resnum> to amino acid <aa>
    mutated_pose = mutate_residue(pose, resnum, aa, repack_radius, sfxn)

    # Score Mutated Pose
    mutant_score = sfxn(mutated_pose)

    # If specified the user, print the breakdown of ddG values into a file
    print_ddG_breakdown(pose, mutated_pose, sfxn, resnum, aa, sc_file)

    # return scores
    return aa, round(mutant_score, 3), round(native_score, 3), round(mutant_score - native_score, 3)

###############################################################################

# @brief Replace the residue at <resid> in <pose> with <new_res> and allows
# repacking within a given <pack_radius>


def mutate_residue(pose, mutant_position, mutant_aa, pack_radius, pack_scorefxn):

    if pose.is_fullatom() == False:
        IOError('mutate_residue only works with fullatom poses')

    test_pose = Pose()
    test_pose.assign(pose)

    # Create a packer task (standard)
    task = TaskFactory.create_packer_task(test_pose)

    # the Vector1 of booleans (a specific object) is needed for specifying the
    #    mutation, this demonstrates another more direct method of setting
    #    PackerTask options for design
    aa_bool = vector1_bool()

    # PyRosetta uses several ways of tracking amino acids (ResidueTypes)
    # the numbers 1-20 correspond individually to the 20 proteogenic amino acids
    # aa_from_oneletter returns the integer representation of an amino acid
    #    from its one letter code
    # convert mutant_aa to its integer representation
    mutant_aa = aa_from_oneletter_code(mutant_aa)

    # mutation is performed by using a PackerTask with only the mutant
    #    amino acid available during design
    # to do this, construct a Vector1 of booleans indicating which amino acid
    #    (by its numerical designation, see above) to allow
    for i in range(1, 21):
        # in Python, logical expression are evaluated with priority, thus the
        #    line below appends to aa_bool the truth (True or False) of the
        #    statement i == mutant_aa
        aa_bool.append(i == mutant_aa)

    # modify the mutating residue's assignment in the PackerTask using the
    #    Vector1 of booleans across the proteogenic amino acids
    task.nonconst_residue_task(
        mutant_position).restrict_absent_canonical_aas(aa_bool)

    # prevent residues from packing by setting the per-residue "options" of
    #    the PackerTask
    center = pose.residue(mutant_position).nbr_atom_xyz()
    for i in range(1, pose.total_residue() + 1):
        dist = center.distance_squared(test_pose.residue(i).nbr_atom_xyz())
        # only pack the mutating residue and any within the pack_radius
        print(i, pack_radius, dist,  pow(float(pack_radius), 2))
        print('##################################################')
        if i != mutant_position and dist > pow(float(pack_radius), 2):
            task.nonconst_residue_task(i).prevent_repacking()

    # apply the mutation and pack nearby residues
    packer = PackRotamersMover(pack_scorefxn, task)
    packer.apply(test_pose)

    return test_pose

###############################################################################
#@brief Print ddG breakdown from the pose
# Extract weighted energies from the native and mutated pose. Calculate the ddG
# of each and print the component-wise ddG vlaues


def print_ddG_breakdown(native_pose, mutated_pose, sfxn, resnum, aa, fn):

    # Extract scores
    tmp_native = native_pose.energies().total_energies(
    ).weighted_string_of(sfxn.weights())
    tmp_mutant = mutated_pose.energies().total_energies(
    ).weighted_string_of(sfxn.weights())

    # Parse out scores
    array_native = list(filter(None, tmp_native.split(' ')))
    array_mutant = list(filter(None, tmp_mutant.split(' ')))

    # Pull out only the scores from these arrays
    native_scores = []
    for i in range(len(array_native)):
        if (i % 2 != 0):
            native_scores.append(float(array_native[i]))

    mutant_scores = []
    for i in range(len(array_mutant)):
        if (i % 2 != 0):
            mutant_scores.append(float(array_mutant[i]))

    # Make a label for the mutation
    native_res = native_pose.residue(int(resnum)).name1()
    mut_label = native_res + str(resnum) + aa

    # Calculate ddG of individual components
    ddGs = []
    ddGs.append(mut_label)
    for i in range(len(mutant_scores)):
        ddG_component = mutant_scores[i] - native_scores[i]
        ddGs.append(round(ddG_component, 3))

    ddGs_str = convert_array_to_str(ddGs)
    with open(fn, 'a') as f:
        f.write(ddGs_str + "\n")

###############################################################################
#@brief Get header for ddG breakdown output
# Save the score labels, to be printed at the top of the output breakdown file


def print_score_labels_to_file(native_pose, sfxn, fn):

    tmp_native = native_pose.energies().total_energies(
    ).weighted_string_of(sfxn.weights())
    array_native = list(filter(None, tmp_native.split(' ')))
    print(array_native)
    labels = []
    labels.append('mutation ')  # Append field for mutation label
    for i in range(len(array_native)):
        if (i % 2 == 0):
            labels.append(array_native[i].translate(':'))

    labels_str = convert_array_to_str(labels)
    print(labels_str)
    with open(fn, 'a') as f:
        f.write(labels_str + "\n")


###############################################################################
#@brief Convert an array to a space deliminted string
# Save the score labels, to be printed at the top of the output breakdown file
def convert_array_to_str(array):

    linestr = ""
    for elem in array:
        if (linestr == ""):
            linestr = linestr + str(elem)
        else:
            linestr = linestr + " " + str(elem)

    return linestr


if __name__ == "__main__":
    main(sys.argv)
