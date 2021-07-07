"""pyrosetta_ddG.py does perform the ddG calculations for this pipeline.

Author: Anders Frederiksen
Contribution: Johanna K.S. Tiemann

Date of last major changes: 2020-04-15

"""

# Standard library imports
from argparse import ArgumentParser
import json
import logging  as log
import os
import sys

# Third party imports
import numpy as np
import pandas as pd
import pyrosetta.rosetta.protocols.membrane
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta import create_score_function
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.toolbox.mutants import mutate_residue as pyrosetta_toolbox_mutants_mutate_residue
from pyrosetta.rosetta.core.scoring import all_scatom_rmsd_nosuper
from pyrosetta.rosetta.core.scoring import all_atom_rmsd_nosuper
from pyrosetta.rosetta.core.scoring import rms_at_corresponding_heavy_atoms
# Local application imports



# Local application imports
log_message="verbose"

if log_message=="verbose":
    log.basicConfig(
        format='%(levelname)s:%(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=log.INFO
    )
elif log_message=="debug":
    log.basicConfig(
        format='%(levelname)s:%(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=log.WARNING
    )
else:
    log.basicConfig(
        format='%(levelname)s:%(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=log.ERROR
    )

logger = log.getLogger(__name__)


def parse_args():
    """
    Argument parser function
    """
    parser = ArgumentParser( description="" )

    # Initiating/setup command line arguments
    parser.add_argument('--in_pdb', '-p',
        type=str,
        help="Input PDB file.", )

    parser.add_argument('--in_span', '-s',
        type=str,
        help="Input spanfile.", )

    parser.add_argument('--outdir', '-o',
        type=str,
        default='.',
        help="Output directory. Files in this directory will be appended. If not present, the directory will be created.", )

    parser.add_argument('--outfile',
        type=str,
        default='ddG.out',
        help="Output filename with pose residue numbering. Default: 'ddG.out'", )

    parser.add_argument('--out_add',
        type=str,
        default='ddG_additional.out',
        help="Output filename of additional information with pose residue numbering. Default: 'ddG_additional.out'", )

    parser.add_argument('--mutfile', '-m',
        type=str,
        default='False',
        help="Rosetta mutfile containing the mutations. If empty, all residues will be mutated into everything.", )

    parser.add_argument('--repack_radius', '-a',
        default=8.0,
        type=float,
        help="Repack the residues within this radius",)

    parser.add_argument('--rescale',
        default=1.0,
        type=float,
        help="Rescale REU to kcal/mol",)

    parser.add_argument('--lowest',
        default=1,
        type=int,
        choices=[1,0],
        help="calculate the mean ddG from the lowest half (min: lowest 3), default: True",)

    parser.add_argument('--include_pH', '-z',
        default=0,
        help="Include pH energy terms: pH_energy and fa_elec. Default false.", )

    parser.add_argument('--pH_value', '-q',
        default=7,
        type=float,
        help="Predict ddG and specified pH value. Default 7. Will not work if include pH is not passed", )

    parser.add_argument('--repeats',
        default=5,
        type=int,
        help="Number of ddG calculations. Default 5.", )

    parser.add_argument('--lipids', '-l',
        default='DLPC',
        choices=['DLPC', 'DMPC', 'DOPC', 'DPPC', 'POPC', 'DLPE', 'DMPE', 
            'DOPE', 'DPPE', 'POPE', 'DLPG', 'DMPG', 'DOPG', 'DPPG', 'POPG'],
        help="Specify lipids from options.", )

    parser.add_argument('--temperature', '-t',
        default=20.0,
        type=float,
        help="Experimental temperature.", )

    parser.add_argument('--lip_has_pore',
        type=str,
        help="File containing the dictionary if residue faces lipids, set to false.", )

    parser.add_argument('--score_function',
        default='franklin2019',
        type=str,
        help="ddG score function. Default=franklin2019, other options are 'ref15', 'mpframework_smooth_fa_2012', 'ref15_memb' ", )

    parser.add_argument('--dump_pdb',
        default=0,
        type=int,
        choices=[1,0],
        help="dumps out the mutated pdb files", )

    parser.add_argument('--repack_protocol',
        default='MP_repack',
        choices=['MP_repack', 'MP_flex_relax_ddG'],
        help="MP repacking algorithm (mainly for benchmarking). Default=MP_repack, other options are 'MP_flex_relax_ddG' ", )

    # parse options
    args = parser.parse_args()

    # Check the required inputs (PDB file, spanfile) are present
    if (not args.in_pdb or not args.in_span):
        sys.exit("Must provide flags '-in_pdb' and '-in_span'! Exiting...")

    if (float(args.pH_value) < 0 or float(args.pH_value) > 14):
        sys.exit("Specified pH value must be between 0-14: Exiting...")

    args.outdir = os.path.abspath(args.outdir)
    os.makedirs(args.outdir, exist_ok=True)
    args.outfile = os.path.join(args.outdir, os.path.basename(args.outfile))
    args.out_add = os.path.join(args.outdir, os.path.basename(args.out_add))
    if args.lowest == 1:
        args.lowest = True
    else:
        args.lowest = False
    if args.dump_pdb == 1:
        args.dump_pdb = True
        args.dump_dir = os.path.join(args.outdir, 'pdbs')
        os.makedirs(args.dump_dir, exist_ok=True)
    else:
        args.dump_pdb = False
        args.dump_dir = False

    return args


def read_mutfile(mutfile):
    mutations = []
    mutate = {}
    if mutfile != 'False':
        with open(mutfile, 'r') as f:
            first_line = next(f).split()[0]
        with open(mutfile, 'r') as f:
            if first_line == 'total':
                save = False
                next(f)
                for lines in f:
                    lines = lines.split()
                    if len(lines) == 1:
                        if save:
                            if "_".join(key) in mutate.keys():
                                old_wt, old_val = mutate["_".join(key)]
                                mutate["_".join(key)] = "_".join(wt), old_val+["_".join(val)]
                            else:
                                mutate["_".join(key)] = "_".join(wt), ["_".join(val)]
                                print(mutate)
                        else:
                            save = True
                        key=[]
                        wt=[]
                        val=[]
                    else:
                        key.append(lines[1])
                        wt.append(lines[0])
                        val.append(lines[2])
        if "_".join(key) in mutate.keys():
            old_wt, old_val = mutate["_".join(key)]
            mutate["_".join(key)] = "_".join(wt), old_val+["_".join(val)]
        else:
            mutate["_".join(key)] = "_".join(wt), ["_".join(val)]
    logger.info(f"mutate: {mutate}")
    for key in mutate.keys():
        mutate[key][1].remove( mutate[key][0])
    logger.info(f"mutate: {mutate}")

    return mutate


def read_lipids_pore(lip_has_pore_file, residue):
    with  open(lip_has_pore_file, 'r') as fp:
        lip_has_pore_dic = json.loads(fp.read())
    #take currently for multi-mutants only the one which occurs most often
    lip_has_pore = []
    for resid in residue.split('_'):
        lip_has_pore.append(lip_has_pore_dic[resid])
    lip_has_pore = max(zip((lip_has_pore.count(item) for item in set(lip_has_pore)), set(lip_has_pore)))[1]

    return lip_has_pore


def mutate_residue_repack(pose, mutant_position, mutant_aa, pack_radius, pack_scorefxn):

    if pose.is_fullatom() == False:
        IOError('mutate_residue only works with fullatom poses')

    test_pose = Pose()
    test_pose.assign(pose)

    # repacking of sidechain atoms within the desired pack radius and mutation
    # into the assigned amino acid
    mutant_position = mutant_position.split("_")
    mutant_aa = mutant_aa.split("_")
    for index, residue in enumerate(mutant_position):
        pyrosetta_toolbox_mutants_mutate_residue(test_pose, int(mutant_position[index]), mutant_aa[index], 
            pack_radius=pack_radius, pack_scorefxn=pack_scorefxn)

    return test_pose


# mutates and repacks using the cartesian algorithm
def mutate_residue_flex_relax_algo(pose, mutant_position, mutant_aa, pack_radius, pack_scorefxn, ddg_bbnbrs=1, cartesian=False, relax_iter=200):
    import time
    from pyrosetta.rosetta.core.pack.task import operation
    
    start_time = time.time()
    logger.warning("Interface mode not implemented (should be added!)")
    
    if cartesian:
        pack_scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreTypeManager.score_type_from_name('cart_bonded'), 0.5)
        #pack_scorefxn.set_weight(atom_pair_constraint, 1)#0.5
        pack_scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreTypeManager.score_type_from_name('pro_close'), 0)
        #logger.warning(pyrosetta.rosetta.basic.options.get_boolean_option('ex1'))#set_boolean_option( '-ex1', True )
        #pyrosetta.rosetta.basic.options.set_boolean_option( 'ex2', True )
    
    if pose.is_fullatom() == False:
        IOError('mutate_residue only works with fullatom poses')

    #Cloning of the pose including all settings
    working_pose = pose.clone()

    #Update residue neighbors for optimal performance
    working_pose.update_residue_neighbors()

    #Prepare mutants
    mutant_position = mutant_position.split("_")
    mutant_aa = mutant_aa.split("_")

    #Select mutant residue
    mutant_selector = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector()
    mutant_selector_list = {}
    for index, residue in enumerate(mutant_position):
        tmp_mutant_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(int(residue))
        mutant_selector_list[residue] = tmp_mutant_selector
        mutant_selector.add_residue_selector(tmp_mutant_selector)
    
    #Select all except mutant
    all_nand_mutant_selector = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector()
    all_nand_mutant_selector.set_residue_selector(mutant_selector)

    #Select neighbors with mutant
    nbr_or_mutant_selector = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector()
    for index, residue in enumerate(mutant_position):
        tmp_nbr_or_mutant_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
        tmp_nbr_or_mutant_selector.set_focus(str(residue))
        tmp_nbr_or_mutant_selector.set_distance(pack_radius)
        tmp_nbr_or_mutant_selector.set_include_focus_in_subset(True)
        nbr_or_mutant_selector.add_residue_selector(tmp_nbr_or_mutant_selector)

    #Select mutant and it's sequence neighbors
    seq_nbr_or_mutant_selector = pyrosetta.rosetta.core.select.residue_selector.PrimarySequenceNeighborhoodSelector(ddg_bbnbrs, ddg_bbnbrs, mutant_selector, False)            

    #Select mutant, it's seq neighbors and it's surrounding neighbors
    seq_nbr_or_nbr_or_mutant_selector = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector()
    seq_nbr_or_nbr_or_mutant_selector.add_residue_selector(seq_nbr_or_mutant_selector)
    seq_nbr_or_nbr_or_mutant_selector.add_residue_selector(nbr_or_mutant_selector)    

    logger.info(f'mutant_selector: {pyrosetta.rosetta.core.select.residue_selector.selection_positions(mutant_selector.apply(working_pose))}')
    logger.info(f'all_nand_mutant_selector: {pyrosetta.rosetta.core.select.residue_selector.selection_positions(all_nand_mutant_selector.apply(working_pose))}')
    logger.info(f'nbr_or_mutant_selector: {pyrosetta.rosetta.core.select.residue_selector.selection_positions(nbr_or_mutant_selector.apply(working_pose))}')
    logger.info(f'seq_nbr_or_mutant_selector: {pyrosetta.rosetta.core.select.residue_selector.selection_positions(seq_nbr_or_mutant_selector.apply(working_pose))}')
    logger.info(f'seq_nbr_or_nbr_or_mutant_selector: {pyrosetta.rosetta.core.select.residue_selector.selection_positions(seq_nbr_or_nbr_or_mutant_selector.apply(working_pose))}')
    

    #Mutate residue and pack rotamers before relax
    #generate packer task
    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.IncludeCurrent())
    
    #Set all residues except mutant to false for design and repacking
    prevent_repacking_rlt = operation.PreventRepackingRLT()
    prevent_subset_repacking = operation.OperateOnResidueSubset(prevent_repacking_rlt, all_nand_mutant_selector, False )
    tf.push_back(prevent_subset_repacking)

    #Assign mutant residue to be designed and repacked
    for index, mutate2aa in enumerate(mutant_aa):
        residue = mutant_position[index]
        single_mutant_selector = mutant_selector_list[residue]
        resfile_comm = pyrosetta.rosetta.protocols.task_operations.ResfileCommandOperation(single_mutant_selector, f"PIKAA {mutate2aa}")
        resfile_comm.set_command(f"PIKAA {mutate2aa}")
        tf.push_back(resfile_comm)

    #Apply packing of rotamers of mutant
    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover()
    packer.score_function(pack_scorefxn)
    packer.task_factory(tf)
    logger.warning(tf.create_task_and_apply_taskoperations(working_pose))
    packer.apply(working_pose)

    #allow the movement for bb for the mutant + seq. neighbors, and sc for neigbor in range, seq. neighbor and mutant
    movemap = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
    movemap.all_jumps(False)
    movemap.add_bb_action(pyrosetta.rosetta.core.select.movemap.mm_enable, seq_nbr_or_mutant_selector)
    movemap.add_chi_action(pyrosetta.rosetta.core.select.movemap.mm_enable, seq_nbr_or_nbr_or_mutant_selector)
    
    #for checking if all has been selected correctly
    logger.info(f"movemap: {movemap.create_movemap_from_pose(working_pose)}")

    #Generate a TaskFactory
    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.IncludeCurrent())
    #tf.push_back(operation.NoRepackDisulfides())

    #prevent all residues except selected from design and repacking
    prevent_repacking_rlt = operation.PreventRepackingRLT()
    prevent_subset_repacking = operation.OperateOnResidueSubset(prevent_repacking_rlt, seq_nbr_or_nbr_or_mutant_selector, True )
    tf.push_back(prevent_subset_repacking)

    # allow selected residues only repacking (=switch off design)
    restrict_repacking_rlt = operation.RestrictToRepackingRLT()
    restrict_subset_repacking = operation.OperateOnResidueSubset(restrict_repacking_rlt , seq_nbr_or_nbr_or_mutant_selector, False)
    tf.push_back(restrict_subset_repacking)


    #Perform a FastRelax
    fastrelax = pyrosetta.rosetta.protocols.relax.FastRelax()
    fastrelax.set_scorefxn(pack_scorefxn)
    logger.info(f'get_scorefxn : {fastrelax.get_scorefxn()}')
    
    if cartesian:
        fastrelax.cartesian(True)
    if relax_iter:
        fastrelax.max_iter(relax_iter)
        
    fastrelax.set_task_factory(tf)
    fastrelax.set_movemap_factory(movemap)
    fastrelax.set_movemap_disables_packing_of_fixed_chi_positions(True)
    
#    logger.info(tf.create_task_and_apply_taskoperations(working_pose))
    pre_time = time.time()
    fastrelax.apply(working_pose)
    post_time = time.time()
    logger.info(f'{post_time - start_time}, vs {pre_time - start_time}')
    return working_pose


def compute_ddG(pose, sfxn, resnum, aa, repack_radius, repack_protocol='MP_repack', add_file=None, flag='', dump_pdb=False, output_dir='.'):
    # Perform Mutation at residue <resnum> to amino acid <aa>
    if repack_protocol == 'MP_flex_relax_ddG':
        mutated_pose = mutate_residue_flex_relax_algo(
            pose, resnum, aa, repack_radius, sfxn, ddg_bbnbrs=1)
    elif repack_protocol == 'MP_ori_design':
        mutated_pose = mutate_residue_original(
            pose, resnum, aa, repack_radius, sfxn)
    elif repack_protocol == 'MP_repack':
        mutated_pose = mutate_residue_repack(
            pose, resnum, aa, repack_radius, sfxn)
    else:
        IOError("wrong repack_protocol assigned - check input")
    # Score Mutated Pose
    mutant_score = sfxn(mutated_pose)
    logger.info(f'all_scatom_rmsd_nosuper: {all_scatom_rmsd_nosuper(pose, mutated_pose)}')
    rmsd = all_atom_rmsd_nosuper(pose, mutated_pose)
    logger.info(f'all_atom_rmsd_nosuper: {rmsd}')
    logger.info(f'rms_at_corresponding_heavy_atoms: {rms_at_corresponding_heavy_atoms(pose, mutated_pose)}')

    if add_file:
        with open(add_file, 'a') as fp:
            variant_name = []
            for indi, var in enumerate(resnum.split('_')):
                variant_name.append(f"{list(pose.sequence())[int(resnum.split('_')[indi])-1]}{resnum.split('_')[indi]}{aa.split('_')[indi]}")
            variant_name = ":".join(variant_name)
            
            score_weights = mutated_pose.energies().total_energies().weighted_string_of(sfxn.weights())
            fp.write(f"{flag}_{variant_name},{round(mutant_score,5)},{rmsd},{score_weights}\n")

    if dump_pdb:
        variant_name = []
        for indi, var in enumerate(resnum.split('_')):
            variant_name.append(f"{list(pose.sequence())[int(resnum.split('_')[indi])-1]}{resnum.split('_')[indi]}{aa.split('_')[indi]}")
        variant_name = "_".join(variant_name)
        mutated_pose.dump_pdb(os.path.join(output_dir, f"{flag}_{variant_name}.pdb"))

    return [aa, mutant_score, rmsd]


def checkIfCalculated(variant, out, out_add):
    if os.path.getsize(out) > 0:
        df = pd.read_csv(out, header=None)
        df = df.rename(columns={0:'variant', 1:'ddG', 2:'std'})
    else:
        df = pd.DataFrame(columns=['variant', 'ddG', 'std'])
    
    if variant in df['variant'].values:
        mean_saved = True
    else:
        mean_saved = False
    
    if os.path.getsize(out_add) > 0:
        df_add = pd.read_csv(out_add, header=None)
        df_add = df_add.rename(columns={0:'variant_raw', 1:'dG', 2:'rmsd'})
    else:
        df_add = pd.DataFrame(columns=['variant_raw', 'dG', 'rmsd'])
    df_add['variant'] = df_add['variant_raw'].apply(lambda x: x.split('_')[1])
    df_add = df_add.loc[df_add['variant']==variant].reset_index(drop=True)
    if len(df_add['variant'])>0:
        return mean_saved, list(map(lambda x,y: [variant[-1],x,y], df_add['dG'], df_add['rmsd'])), len(df_add['variant'])
    else:
        return mean_saved, [], 0


def main():
    """
    Main function called as default at the end.
    """
    # get user input arguments
    args = parse_args()
    
    # Initialize Rosetta options from user options. Enable pH mode if applicable
    rosetta_options = (f'-mp:setup:spanfiles {args.in_span}'
                       f' -mp:lipids:temperature {args.temperature}'
                       f' -mp:lipids:composition {args.lipids}'
                  #     f' -run:constant_seed'
                       f' -in:ignore_unrecognized_res'
                       f' -pH_mode -value_pH {args.pH_value}')

    mutations_dic = read_mutfile(args.mutfile)

    for residue in mutations_dic.keys():
        lip_has_pore = read_lipids_pore(args.lip_has_pore, residue)
        rosetta_options += f' -mp:lipids:has_pore {lip_has_pore}'

        # Initialize Rosetta based on user inputs
        pyrosetta.init(extra_options=rosetta_options)

        # Create an energy function
        sfxn = pyrosetta.rosetta.core.scoring.ScoreFunction()
        # Create a smoothed membrane full atom energy function (pH 7 calculations)
        sfxn = create_score_function(args.score_function)
        if args.score_function == 'franklin2019':
            logger.info(f'get nonzero_weighted scoretypes: {sfxn.get_nonzero_weighted_scoretypes()}')
            logger.info(f'fa_water_to_bilayer before: {sfxn.get_weight(pyrosetta.rosetta.core.scoring.ScoreType.fa_water_to_bilayer)}')
            sfxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.fa_water_to_bilayer,1.5) 
            logger.info(f'fa_water_to_bilayer after: {sfxn.get_weight(pyrosetta.rosetta.core.scoring.ScoreType.fa_water_to_bilayer)}')

        # Load Pose, & turn on the membrane
        pose = pyrosetta.rosetta.core.import_pose.pose_from_file(args.in_pdb)

        # Add Membrane to Pose
        add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( args.in_span)
        add_memb.apply(pose)

        # Setup in a topology based membrane
        init_mem_pos = pyrosetta.rosetta.protocols.membrane.MembranePositionFromTopologyMover()
        init_mem_pos.apply(pose)

        #Get WT & mutants
        native_aa, mut_aas = mutations_dic[residue]

        # Calculate ddGs for WT 
        wt_variant = f"{native_aa}{residue}{native_aa}"
        mean_saved, repacked_wt_list, counts = checkIfCalculated(wt_variant, args.outfile, args.out_add)
        for repeat in range(args.repeats-counts):
            logger.warning(f'repacking protocol: {args.repack_protocol}')
            repacked_wt_list.append(compute_ddG(pose, sfxn, residue, native_aa, args.repack_radius, 
                repack_protocol=args.repack_protocol, add_file=args.out_add, flag='WT', dump_pdb=args.dump_pdb, output_dir=args.dump_dir))
        logger.info(repacked_wt_list)
        

        # Calculate ddGs for mutants
        with open(args.outfile, 'a') as fp:
            repacked_mut_list_all = []
            for index, mutant in enumerate(mut_aas):
                repacked_mut_list = []
                var_variant = f"{native_aa}{residue}{mutant}"
                mean_saved, repacked_mut_list, counts = checkIfCalculated(var_variant, args.outfile, args.out_add)
                if mean_saved:
                    repacked_mut_list_all.append([mutant, repacked_mut_list])
                else:
                    for repeat in range(args.repeats-counts):
                        logger.warning(f'repacking protocol: {args.repack_protocol}')
                        repacked_mut_list.append(compute_ddG(pose, sfxn, residue, mutant, args.repack_radius, 
                            repack_protocol=args.repack_protocol, add_file=args.out_add, flag='MUT', dump_pdb=args.dump_pdb, output_dir=args.dump_dir))
                    repacked_mut_list_all.append([mutant, repacked_mut_list])
                    
                    # get variant name
                    variant_name = []
                    logger.info(mutant)
                    for indi, var in enumerate(mutant.split('_')):
                        variant_name.append(f"{repacked_wt_list[0][0].split('_')[indi]}{residue.split('_')[indi]}{var}")
                    variant_name = ":".join(variant_name)
                    logger.info(f'variant_name: {variant_name}')

                    #get WT ddGs
                    wt_dgs = [wt[1] for wt in repacked_wt_list]
                    wt_rmsd = [wt[2] for wt in repacked_wt_list]

                    #get MT ddGs
                    mut_dgs = [mut[1] for mut in repacked_mut_list]
                    mut_rmsd = [mut[2] for mut in repacked_mut_list]
                    #calculate mean & std
                    ddg = []
                    for mut_dg in mut_dgs:
                        ddg.append((mut_dg - np.mean(wt_dgs))/args.rescale)
                    logger.info(f'wt_dgs: {wt_dgs}')
                    logger.info(f'wt_rmsd: {wt_rmsd}')
                    logger.info(f'mut_dgs: {mut_dgs}')
                    logger.info(f'mut_rmsd: {mut_rmsd}')
                    logger.info(f'ddg: {ddg}')
                    logger.info(f"{variant_name},{round(np.mean(ddg), 3)},{round(np.std(ddg), 3)}\n")


                    if args.lowest:
                        ddg = np.array(ddg)
                        small_k = int(args.repeats/2)+1
                        smallest_idx = np.argpartition(ddg, small_k)
                        ddg = ddg[smallest_idx[:small_k]]
                        logger.info(f'lowest ddg: {ddg}')

                    fp.write(f"{variant_name},{round(np.mean(ddg), 3)},{round(np.std(ddg), 3)}\n")


            logger.info(repacked_mut_list_all)

if __name__ == '__main__':
    main()