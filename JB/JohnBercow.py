#!/software/containers/john_bercow.sif

# Generates an IDT-ready .xlsx file for ordering eBlocks from a folder of PDBs and/or a FASTA file.
# Compatible overhangs for Golden Gate cloning are added automatically.
# Reverse translation is performed with Ryan's Domesticator.
# Sequences are queried against IDT (courtesy of Ryan) to ensure synthesiability.
# RECOMMENDED: check your GG assemblies at https://goldengate.neb.com/#!/
# Wondering why the script is called John Bercow? https://www.youtube.com/watch?v=VYycQTm2HrM&ab_channel=TheSun

# ============================================
# TODO
# ============================================
# Splitting beyond 2 fragments
# Automatic layout for fragmented hetero-oligomers

# ============================================
# LIBRARIES
# ============================================
import sys
import glob
import os
import re
import argparse
import datetime; date = datetime.datetime.now().strftime('%Y_%m_%d')
import numpy as np
import pandas as pd
from Bio import SeqIO, PDB, SeqUtils, Seq
import domesticator # Ryan's domesticator with minor modifications.
import idt # IDT API - from Ryan.

# ============================================
# HARD-CODED DEFINITIONS
# ============================================
# Restriction enzyme and Cterm/Nterm protein tags of the different GG vectors.
# Enzyme / Nterm-tag / Cterm-tag / 5'-sticky / 3'-sticky / description
vectors = {
    'LM0627':['BsaI','MSG','GSGSHHWGSTHHHHHH','agga','ggttcc', 'C-term SNAC-His'],
    'LM0670':['BsaI','MSG','GSHHHHHH','agga','ggttcc', 'C-term His'],
    'LM0671':['BsaI','MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFARYPDHMKQHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYKGSSG','GSHHHHHH','agga','ggttcc', 'N-term sfGFP and C-term His'],
    'LM0673':['BsaI','MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFARYPDHMKQHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYKGGSHHWSSG','GSHHHHHH','agga','ggttcc', 'N-term sfGFP-SNAC and C-term His'],
    'LM1371':['BsaI','MSHHHHHHSG','GS','agga','ggttcc', 'N-term His'],
    'LM1425':['BsaI','MSG','GSGLNDIFEAQKIEWHESHHHHHH','agga','ggttcc', 'C-term AviTag-His'],
    'BW1001':['BsaI','MSG','GSSGSGGSGGGGSGGSSSGGVTGYRLFEEILGSHHHHHH','agga','ggttcc', 'C-term smBiT-His'],
    'BW1002':['BsaI','MSG','GSSGSGGSGGGGSGGSSSGGVTGYRLFEEILGSTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTEGSHHHHHH','agga','ggttcc', 'C-term smBiT-GB1-His'],
    'BW1003':['BsaI','MSG','GSSGSGGSGGGGSGGSSSGGVFTLEDFVGDWEQTAAYNLDQVLEQGGVSSLLQNLAVSVTPIQRIVRSGENALKIDIHVIIPYEGLSADQMAQIEEVFKVVYPVDDHHFKVILPYGTLVIDGVTPNMLNYFGRPYEGIAVFDGKKITVTGTLWNGNKIIDERLITPDGSMLFRVTINSGSHHHHHH','agga','ggttcc', 'C-term lgBiT'],
    'BW1004':['BsaI','MSGHHHHHHGSVTGYRLFEEILGGSGSGGSGGGGSGGSSSGG','GS','agga','ggttcc', 'N-term His-smBiT'],
    'BW1005':['BsaI','MSGHHHHHHGSTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTEGSVTGYRLFEEILGGSGSGGSGGGGSGGSSSGG','GS','agga','ggttcc', 'N-term His-GB1-smBiT'],
    'BW1006':['BsaI','MSGHHHHHHGSVFTLEDFVGDWEQTAAYNLDQVLEQGGVSSLLQNLAVSVTPIQRIVRSGENALKIDIHVIIPYEGLSADQMAQIEEVFKVVYPVDDHHFKVILPYGTLVIDGVTPNMLNYFGRPYEGIAVFDGKKITVTGTLWNGNKIIDERLITPDGSMLFRVTINSGGSGSGGSGGGGSGGSSSGG','GS','agga','ggttcc', 'N-term His-lgBiT'],
}

vec_str = '\n'.join([f' - {k} ({v[0]}): {v[5]} | {v[1]}[...]{v[2]}' for k, v in vectors.items()]) # for argparse description

# All lab-available plasmids
gg_vectors = glob.glob('./entry_vectors/*.fa')
entry_vectors = {}
for v in gg_vectors:
    records = SeqIO.parse(v, 'fasta')
    for record in records:
        entry_vectors[record.id.split('_')[0]] = str(record.seq.lower())

# Type IIS restriction enzymes cut sites. Both the  forward and reverse complement sequences (each 5'->3') are indicated.
avoid_seq = {
    'BsaI':['GGTCTC', 'GAGACC'],
    'SapI':['GCTCTTC', 'GAAGAGC']
}

# GG adapters.
gg_adapters = {
    'BsaI':{
        '5prime':'atactacggtctca',
        '3prime':'cgagaccgtaatgc',
        },
    'SapI':{
        '5prime':'atactacgctcttcg',
        '3prime':'cgaagagcgtaatgc',
        }
}

# FW, RV, N_spacer, N_sticky -- for generating the plasmid maps
cuts = {
    'BsaI':['ggtctc', 'gagacc', 1, 4],
    'SapI':['gctcttc', 'gaagagc', 1, 3]
}

base_pairing = {'C':'G', 'G':'C', 'A':'T', 'T':'A'}
# Pool of 4bp sticky ends. Used for spliting genes into 2 fragments if necessary. ONLY FOR LM VECTORS!!!
# Overhangs (4 bp) with high T4 ligase fidelity. From Ligase Fidelity Viewer (v2) https://ggtools.neb.com/viewset/run.cgi
sticky_4bp = ['AAGG', 'ACTC', 'AGGA', 'ATCA', 'GCCG', 'CTGA', 'GCGA', 'GGAA', 'GTTT']
sticky_4bp_comp = [''.join([base_pairing[b] for b in st]) for st in sticky_4bp] # complementary sequences
sticky_4bp_revcomp = [s[::-1] for s in sticky_4bp_comp] # reverse complementary sequences
sticky_4bp = sticky_4bp + sticky_4bp_revcomp # by convention, consider 5'->3' sequences only.
sticky_4bp_already_used = ['AGGA', 'TTCC', 'TCCT', 'GGAA'] # remove these sites from the pool
sticky_4bp = list(set([st for st in sticky_4bp if st not in sticky_4bp_already_used]))

# ============================================
# ARGUMENTS
# ============================================
parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=" * Generates an IDT-ready .xlsx file for ordering eBlocks from a folder of PDBs and/or a concatenated FASTA file.\n"
                    " * Appropriate overhangs for Golden Gate cloning into entry vector(s) of interest are added automatically.\n"
                    " * RECOMMENDED: check your GG assemblies at https://goldengate.neb.com/#!/\n"
                    " * Wondering why the script is called John Bercow? https://www.youtube.com/watch?v=VYycQTm2HrM&ab_channel=TheSun\n"
                    "\n"
                    " * EXAMPLE AVAILABLE ENTRY VECTORS:\n"
                    " *** see ./entry_vectors/ for the FULL list ***\n"
                   f"{vec_str}\n"
        )
# REQUIRED
parser.add_argument(
        '--order_pdbs',
        help='path to a folder containing the PDBs you want to order. Either this option or --order_fasta needs to be specified. Both can specified.',
        action='store',
        type=str
        )
parser.add_argument(
        '--order_fasta',
        help='path to a FASTA file containing the AA sequences you want to order. Either this option or --order_pdbs needs to be specified. Both can be specified.',
        action='store',
        type=str
        )
parser.add_argument(
        '--order_name',
        help='name of the order (appended to output files). NB date, plate IDs, organism, and enzyme are added automatically.',
        action='store',
        type=str,
        required=True
        )
parser.add_argument(
        '--gg_vector',
        help='name of target vector(s) for Golden Gate cloning (determines the DNA adapters). Also determines the AA tags appended to the design in the FASTA output. Multiple vectors can be specified as space-separated values, but they need to be compatible (i.e. use the same enzyme and same overhangs).',
        action='store',
        required=True,
#        choices=vectors.keys(),
        nargs='+'
        )
parser.add_argument(
        '--species',
        help='codon optimisation will be performed for this species (e.g. e_coli, s_cerevisiae, h_sapiens, etc...)',
        action='store',
        type=str,
        required=True
        )
parser.add_argument(
        '--design_prefix',
        help='designs get IDs with this prefix (e.g. LM0001, LM0002, etc...)',
        action='store',
        type=str,
        required=True
        )
parser.add_argument(
        '--design_id',
        help='increment design indices from this number.',
        action='store',
        type=int,
        required=True
        )

# OPTIONAL
parser.add_argument(
        '--skip_idt_query',
        help="skip IDT website query that checks synthesisability of eBlocks.",
        action='store_true',
        )
parser.add_argument(
        '--idt_score',
        help="IDT complexity score threshold for accepting a reverse translated sequence (Default: 7). 0-7 is green, 7-15 is yellow, >15 is red and IDT will not make it.",
        action='store',
        type=float,
        default=7
        )
parser.add_argument(
        '--starting_kmers_weight',
        help="starting value for the kmers_weight setting of Domesticator (Default: 10). This parameter is linearly ramped (up to 100) over --n_domesticator_steps.",
        action='store',
        type=int,
        default=10
        )
parser.add_argument(
        '--n_domesticator_steps',
        help="maximum number of Domesticator steps attempted (Default: 10). The kmers_weight parameter (which increases synthesiability of repetitive sequences) is linearly ramped up to 100 over this number of steps.",
        action='store',
        type=int,
        default=10
        )
parser.add_argument(
        '--max_attempts',
        help="maximum number of reverse translation attempts at each Domesticator step (Default: 20). Since Domesticator is stochastic, re-running the optimisation problem with the same parameters can lead to different solutions.",
        action='store',
        type=int,
        default=20
        )
parser.add_argument(
        '--max_length',
        help="maximum length of eBlocks. IDT's maximum is 1500 bp, but less can be specified if sequence complexity is a issue for synthesis and you want to force the generation of smaller fragments.",
        action='store',
        type=int,
        default=1500
        )
parser.add_argument(
        '--print_hto',
        help="print heterooligomers to stdout.",
        action='store_true',
        )
parser.add_argument(
        '--no_layout',
        help="do not apply automated layout formatting.",
        action='store_true',
        )
parser.add_argument(
        '--no_plasmids',
        help="do not generate the cloned plasmid maps.",
        action='store_true',
        )
parser.add_argument(
        '--verbose',
        help="increase the verbosity of th e output (recommended).",
        action='store_true',
        )
parser.add_argument(
        '--echo',
        help="generates outputs formated as 384w plates (for ordering into ECHO-qualified plates). Currently only available in conjuction with the no_layout option.",
        action='store_true',
        )

args = parser.parse_args()

# Sanity checks.
if args.order_pdbs == None and args.order_fasta == None:
    print('Either a folder containing PDBs or a FASTA file need to be specified. System exiting...')
    sys.exit()

if args.order_pdbs == None:

    if '/' in args.order_fasta:
        output_folder = '/'.join(args.order_fasta.split('/')[:-1]) + '/'

    else:
        output_folder = './'

else:
    output_folder = args.order_pdbs

if args.max_length < 600:
    print('IDT cannot synthesise eBlocks outside of the 300-1500 bp range. Splitting a gene of less than 600 bp will generate fragments that are smaller than 300 bp each. System exiting...')
    sys.exit()

if args.idt_score > 10.0:
    print('IDT will not synthesise anything with a complexity score above 10. System exiting...')
    sys.exit()

try:
    enzyme = np.unique([vectors[gg_v][0] for gg_v in args.gg_vector])
    if len(enzyme) > 1:
        print(f'The destination vectors ({", ".join(args.gg_vector)}) are not compatible. System exiting...')
        sys.exit()

    else:
        enzyme = enzyme[0]
        sticky5prime, sticky3prime = vectors[args.gg_vector[0]][3], vectors[args.gg_vector[0]][4]

except:
    enzyme = 'BsaI'
    sticky5prime, sticky3prime = vectors['LM0627'][3], vectors['LM0627'][4]
    print('Assuming LM0627-like overhangs for the selected vector(s).')
    print('!!! DO NOT ORDER IF THIS IS NOT THE CASE AS YOU WON\'T BE ABLE TO CLONE YOUR EBLOCKS !!!')


if args.echo and not args.no_layout:
    print('The --echo  option is currently only available in conjuction with --no_layout. Switching to this...')
    args.no_layout = True

if not args.skip_idt_query: # Get/make IDT API credentials to query sequence complexity.
    idt.user_info_file, idt.token_file = idt.use_dir("~/idt_credentials")
    idt_user_info = idt.get_user_info(idt.user_info_file)

filename = f'{output_folder}{date}_{args.order_name}_{args.species}_{enzyme}'

# ============================================
# FUNCTIONS
# ============================================

def adjust_for_eblock(dna_seq, max_length, gg_int_adapters, avoid_seqs):
    '''
    Check size of DNA sequences and either pad it (if too short), or split it (if too long).
    '''
    # Check if the sequence fits the the size limits
    if len(dna_seq) <= max_length:

        if len(dna_seq) < 300:
            print(f'  [!] DNA sequence is too short to be ordered as an eBlock ({len(dna_seq)} vs. 300 bp). Adding some padding to the sequence...')
            pad_length = ((300 - (len(dna_seq))) // 2 )
            extra = 300 - len(dna_seq) - (2 * pad_length)
            if extra < 0 :
                extra = 0

            # Pad sequence
            pad_nocut = False
            while pad_nocut == False:

                pad5prime = ''.join(np.random.choice(['A','T','C','G'], size=pad_length + extra))
                pad3prime = ''.join(np.random.choice(['A','T','C','G'], size=pad_length))

                pad_nocut = True
                for a_seq in avoid_seqs:
                    # if any cut sequence is found in the padding, try again
                    if (a_seq in pad5prime) or (a_seq in pad3prime):
                        pad_nocut = False

                    else:
                        pad_nocut = True

            dna_seq =  pad5prime + dna_seq + pad3prime
            dna_fragments = {'':dna_seq}

        else:
            dna_fragments = {'':dna_seq}

    else:
        if len(dna_seq) > (2 * args.max_length - ((14+4)*2)):
            print(f'  [!] Lengths of the fragments after 2-way splitting will be too long for eBlock synthesis (>{int(len(dna_seq)/2)} bp each). Splitting into more than 2 fragments is currently not supported.')
            print(f'  [!] This design cannot currently be ordered as eBlocks. Ignoring it and moving on...')
            print(f'  ####################################################################################')
            dna_fragments = {'':None}

        else:
            print(f'  [!] DNA sequence is longer than the maximum specified ({len(dna_seq)} vs. {args.max_length} bp). Splitting the sequence...')

            middle = len(dna_seq) / 2
            ligation_sites = {}
            for st in sticky_4bp:
                matching_pos = np.array([match.start() for match in re.finditer(st, dna_seq)])

                if len(matching_pos) == 0:
                    # print(f'  No matche found for {st}')
                    pass
                else:
                    closest_to_middle = matching_pos[np.argmin(np.abs(matching_pos - middle))]
                    ligation_sites[closest_to_middle] = st
                    # print(f'  Closest ligation site to the middle for {st} found at {closest_to_middle}')

            if len(ligation_sites) == 0:
                print('  [!] No cut site found for this sequence. This eBlock will be too long for synthesis. Ignoring it and moving on...')
                print(f'  ################################################################################################################')
                dna_fragments = {'':None}

            else:
                possible_sites = np.array(list(ligation_sites.keys()))
                best = possible_sites[np.argmin(np.abs(possible_sites - middle))]
                print(f'  * Closest ligation site to the middle of the sequence is {ligation_sites[best]} at {best}')

                # Split DNA fragment and add new GG adapters
                Nterm_dna, Cterm_dna = dna_seq[:best], dna_seq[best+4:]
                dna_frag_x = Nterm_dna + ligation_sites[best] + gg_adapters[enzyme]['3prime']
                dna_frag_y = gg_adapters[enzyme]['5prime'] + ligation_sites[best] + Cterm_dna
                dna_fragments = {'x':dna_frag_x,'y':dna_frag_y}

    return dna_fragments

# ============================================
# GET AA SEQUENCES
# ============================================
aa_sequences = {}

if args.order_pdbs is not None:
    pdbs = sorted(glob.glob(f'{args.order_pdbs}*.pdb'))
    print(f'Extracting sequences from {len(pdbs)} PDBs...')

    for pdb in pdbs:
        pdb_name = pdb.split('/')[-1].replace('.pdb', '')
        parser = PDB.PDBParser(PERMISSIVE=1, QUIET=True)
        structure = parser.get_structure('design', pdb)
        chain2seq = {}

        for model in structure:
             for chain in model:
                chain2seq[chain.id] = ''

                for residue in chain:
                    chain2seq[chain.id] += SeqUtils.IUPACData.protein_letters_3to1[residue.resname.capitalize()]

        seq2chain = {v:k for k, v in chain2seq.items()} # removes duplicates for homo-oligomers

        if len(seq2chain) > 1:
            is_heterooligomer = True

        else:
            is_heterooligomer = False

        for seq, chain in seq2chain.items():
            if is_heterooligomer:
                aa_sequences[pdb_name + '_' + chain + '_isheterooligomer'] = seq  # for renaming later on.

            else:
                aa_sequences[pdb_name] = seq

if args.order_fasta is not None:
    print('Extracting sequences from FASTA file...')
    fasta_sequences = list(SeqIO.parse(args.order_fasta, 'fasta'))
    fasta_names = [fasta.id for fasta in fasta_sequences]

    for fasta in fasta_sequences:

        if fasta.id.split('_')[-1] in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ': # check if name has '_X' to indicate different chains.
            base_name = '_'.join(fasta.id.split('_')[:-1])
            num_bn = np.sum([True if base_name in fn else False for fn in fasta_names])

            if num_bn > 1: # check if base name appears more than once.
                aa_sequences[fasta.id + '_isheterooligomer'] = str(fasta.seq)

            else:
                aa_sequences[fasta.id] = str(fasta.seq)

        else:
            aa_sequences[fasta.id] = str(fasta.seq)

print(f'Number of AA sequences extracted: {len(aa_sequences)}\n')

if len(aa_sequences) == 0:
    print('No AA sequences found. System exiting...')
    sys.exit()

# Option for printing hetero-oligomeric sequences -- useful for debugging.
if args.print_hto:
    print('HETEROOLIGOMERS:')
    heterooligos = np.unique(['_'.join(n.split('_')[:-2]) for n in aa_sequences.keys() if '_isheterooligomer' in n])
    for hto in heterooligos:
        for k, v in aa_sequences.items():
            if hto in k:
                print(f'>{k}\n{v}')
        print('-----')
    print('\n')

# Check for duplicates
all_aa_seq = len(aa_sequences)
unique_aa_seq = len(np.unique(list(aa_sequences.values())))
if all_aa_seq != unique_aa_seq:
    print(f'[!] Found duplicated sequences ({unique_aa_seq} unique sequences vs. {all_aa_seq} total sequences):')

    visited = set()
    dup = [x for x in aa_sequences.values() if x in visited or (visited.add(x) or False)]
    duplicates = {d:[] for d in dup}
    for k, v in aa_sequences.items():
        if v in dup:
            duplicates[v].append(k)

    for seq, ids in duplicates.items():
        for id in ids:
            print(f'>{id}')
        print(seq + '\n-----')

    print('ERROR: Duplicates. Verify your sequences and retry. System exiting...')
    sys.exit()

# ============================================
# GENERATE eBLOCKS
# ============================================
eblocks = {
    'design_name':[],
    'Sequence':[],
    'length_eblock':[],
    'readin_order':[],
    'design_aa_seq':[],
    **{f'{gg_v}_cloned_plasmid_seq':[] for gg_v in args.gg_vector},
    **{f'ORF_from_{gg_v}':[] for gg_v in args.gg_vector},
    **{f'exp_aa_seq_from_{gg_v}':[] for gg_v in args.gg_vector},
    'idt_score':[]
    }

skipped = {}
for i, (design_name, aa_seq) in enumerate(aa_sequences.items()):

    print(f'>>> [{i+1}/{len(aa_sequences)}] {design_name.replace("_isheterooligomer", "")}: {aa_seq}')

    output = True
    n_dom_steps = 0
    kmers_weight_settings = np.linspace(args.starting_kmers_weight, 100, args.n_domesticator_steps)
    idt_score = 10.0
    while (idt_score >= args.idt_score) and (n_dom_steps < args.n_domesticator_steps) and (output == True):

        # Reverse translate
        print(f'  Reverse translating (optimising for {args.species} with kmers_weight@{kmers_weight_settings[n_dom_steps]})...')

        dna_seq = domesticator.reverse_translate(
                                    aa_seq,
                                    kmers_weight=kmers_weight_settings[n_dom_steps],
                                    cai_weight=1.0,
                                    hairpins_weight=1.0,
                                    max_tries=args.max_attempts,
                                    species=args.species,
                                    avoid=avoid_seq[enzyme]
                                    )
        n_dom_steps += 1

        # Double-check for cut sites.
        for seq in avoid_seq[enzyme]:
            if seq in dna_seq.upper():
                output = False
                skipped[design_name + f' [ERROR: {enzyme} site(s) detected]'] = aa_seq
                print(f'  {dna_seq}')
                print(f'  [!] {enzyme} site(s) detected. Problem with reverse translation!')
                print(f'  ############################################################')

            else:
                output = True

        # Add GG adapters and vector-specific sticky ends to DNA sequence.
        dna_seq = gg_adapters[enzyme]['5prime'] + sticky5prime + dna_seq + sticky3prime + gg_adapters[enzyme]['3prime']

        # Check that the sequence fits an eBlock. If not, pad or split the sequence to fall within the 300-1500 bp limits.
        dna_fragments = adjust_for_eblock(dna_seq, args.max_length, gg_adapters[enzyme], avoid_seq[enzyme])

        if None in dna_fragments.values():
            output = False
            skipped[design_name + f' [ERROR: too large to be ordered as eBlocks]'] = aa_seq

        else:

            if args.skip_idt_query:
                idt_scores = {}
                for f, dna_frag in dna_fragments.items():
                    idt_scores[f] = 0 # to escape while loop

                idt_score = np.max(list(idt_scores.values()))
                print(f'  Skipping IDT query.')

            else:
                # Check the synthesiability score(s) of the fragments.
                idt_scores = {}
                for f, dna_frag in dna_fragments.items():
                    idt_scores[f] = idt.total_score(dna_frag, idt_user_info)

                idt_score = np.max(list(idt_scores.values())) # for splitted genes, take the worst score of the two.

                if idt_score >= args.idt_score:
                    print(f'  [!] DNA sequence(s) failed the complexity test from IDT (Total Complexity Score = {idt_score:.1f}). Trying again...')

                else:
                    print(f'  DNA sequence(s) passed the complexity test from IDT (Total Complexity Score = {idt_score:.1f}).')


        if n_dom_steps >= args.n_domesticator_steps:
            output = False
            skipped[design_name + f' [ERROR: IDT cannot make it]'] = aa_seq
            print(f'  [!] Maximum number of Domesticator steps reached ({args.n_domesticator_steps}).')
            print(f"  [!] This design could not be reverse translated to IDT's specifications. Ignoring it and moving on...")
            print(f'  #####################################################################################################')

    # Print/save results.
    if output == False:
        pass

    else:
        print(f'  Final eBlock(s):')
        for f, dna_frag in dna_fragments.items():
            print(f'  - eBlock {f.upper()} ({len(dna_frag)} bp): {dna_frag}')

            
        # Save results.
        for f, dna_frag in dna_fragments.items():
            eblocks['design_name'].append(design_name)
            eblocks['Sequence'].append(dna_frag)
            eblocks['length_eblock'].append(len(dna_frag))
            eblocks['readin_order'].append(str(i) + f)
            eblocks['design_aa_seq'].append(aa_seq)
            eblocks['idt_score'].append(idt_scores[f])

            
            # Generate cloned plasmids, ORFs and expression products
            for gg_v in args.gg_vector:
                
                # Enzyme-specific cut characteristics.
                fw, rv, n_spacer, n_sticky = cuts[enzyme]

                # GG cloning.
                entry_vector = entry_vectors[gg_v]
                vector_5prime, vector_3prime = entry_vector.find(fw) \
                                                + len(fw) \
                                                + n_spacer, entry_vector.find(rv) \
                                                - n_spacer

                insert = dna_seq.lower()

                # Find eBlock section with cut sites facing in the correct directions.
                # (Necessary for cases where cut sites are accidentlly also present in the padding regions.)
                fw_locations = np.array([x.span()[0] for x in re.finditer(fw, insert)])
                rv_locations = np.array([x.span()[0] for x in re.finditer(rv, insert)])

                fw_idx = []
                rv_idx = []
                delta_bp = []
                for i, f in enumerate(fw_locations):
                    for j, r in enumerate(rv_locations):
                        delta_bp.append(r - f)
                        fw_idx.append(i)
                        rv_idx.append(j)

                delta_bp = np.array(delta_bp)
                correct_idx = np.argwhere(delta_bp==delta_bp[delta_bp>=3*len(aa_seq)].min())[0][0]

                insert_5prime = fw_locations[fw_idx[correct_idx]] \
                                + len(fw) \
                                + n_spacer \
                                + n_sticky

                insert_3prime = rv_locations[rv_idx[correct_idx]] \
                                - n_spacer \
                                - n_sticky

                assembled_plasmid = entry_vector[:vector_3prime] \
                                        + insert[insert_5prime:insert_3prime] \
                                        + entry_vector[vector_5prime:]
                
                eblocks[f'{gg_v}_cloned_plasmid_seq'].append(assembled_plasmid.lower())
                
                
                # Identify the ORF that contains the insert.
                # Search for the LONGEST START-STOP span that contains the insert sequence.
                plasmid_seq = assembled_plasmid.lower()
                starts = np.array([s.start() for s in re.finditer('atg', plasmid_seq)])
                ends = np.array(sorted([e.end() for e in re.finditer('tag', plasmid_seq)] 
                                       + [e.end() for e in re.finditer('taa', plasmid_seq)] 
                                       + [e.end() for e in re.finditer('tga', plasmid_seq)]))
                current_stop = 0

                for s in starts:
                    inframe_stops = ends[np.logical_and(ends>s, (ends-s)%3==0)]

                    if len(inframe_stops) > 0:
                        if s > current_stop:
                            current_stop = inframe_stops[0]
                            coding_seq = plasmid_seq[s:current_stop]

                            if insert[insert_5prime:insert_3prime] in coding_seq:
                                ORF = coding_seq
                                exp_product = str(Seq.Seq(ORF).translate())
                
                eblocks[f'ORF_from_{gg_v}'].append(ORF)
                eblocks[f'exp_aa_seq_from_{gg_v}'].append(exp_product)
        
                print(f'  - Expression product (if cloned into {gg_v} with {enzyme}): {exp_product}')
          
    print('\n')
    

df = pd.DataFrame.from_dict(eblocks)

# ============================================
# PLATE FORMATTING
# ============================================
# For 'no_layout' option (just fill plate from A1->H12)
w96 = [str(p) + '_' + r + str(c) for p in range(1,10) for r in 'ABCDEFGH' for c in range(1,13)]

# For 384w formatting (no layout)
w384 = [str(p) + '_' + r + str(c) for p in range(1,10) for r in 'ABCDEFGHIJKLMNOP' for c in range(1,25)]
w384_df = pd.DataFrame([r + str(c) for r in 'ABCDEFGHIJKLMNOP' for c in range(1,25)], columns=['Well Position'])

# For single fragments plate layouts -- staggered arrangement for gel loading.
w96_zigzag = []
for p in range(1, 10):
    for pair in ['AB', 'CD', 'EF', 'GH']:
        for c in range(1,13):
            for r in pair:
                w96_zigzag.append(str(p) + '_' + r + str(c))


# For ECHO transfer file generation
w96_zigzag_single = []
for pair in ['AB', 'CD', 'EF', 'GH']:
    for c in range(1,13):
        for r in pair:
            w96_zigzag_single.append(r + str(c))


# For dual fragments and hetero-dimer layouts -- two horizontal blcoks of 48.
w96_2blocks = []
for p in range(1, 10):
    for pair in ['AE', 'BF', 'CG', 'DH']:
        for c in range(1,13):
            for r in pair:
                w96_2blocks.append(str(p) + '_' + r + str(c))

# For triple fragments and hetero-trimer layouts -- three vertical blocks of 32.
w96_3blocks = []
for p in range(1, 10):
    for group in ['1 5 9', '2 6 10', '3 7 11', '4 8 12']:
        for r in 'ABCDEFGH':
            for c in group.split():
                w96_3blocks.append(str(p) + '_' + r + c)

# For quadruple fragments and hetero-tetramer layouts -- four horizontal blocks of 24.
w96_4blocks = []
for p in range(1, 10):
    for pair in ['AE', 'BF', 'CG', 'DH']:
        for c in range(1,13):
            for r in pair:
                w96_4blocks.append(str(p) + '_' + r + str(c))

# Find fragments and heterooligomers.
design_names = df['design_name'].values
df['order_name'] = len(df) * [args.order_name]
df['is_frag'] = df['readin_order'].apply(lambda x: True if ('x' in x) or ('y' in x) else False )
df['is_hto'] = df['design_name'].apply(lambda x: True if '_isheterooligomer' in x else False)

# Split eBlocks into categories for individual reformatting.
single_eblocks = df[(df['is_frag']==False) & (df['is_hto']==False)]
frag_eblocks = df[(df['is_frag']==True) & (df['is_hto']==False)]
hto_eblocks = df[(df['is_frag']==False) & (df['is_hto']==True)]
htd_frag_eblocks = df[(df['is_frag']==True) & (df['is_hto']==True)]

n_designs = int(args.design_id)
n_plates = int(0)

 # Reformat single eBlocks.
if len(single_eblocks) > 0:
    design_ids = [args.design_prefix + str(i).zfill(4) for i in range(n_designs, n_designs + len(single_eblocks))]
    single_eblocks['design_id'] = design_ids
    single_eblocks['position'] = w96_zigzag[:len(single_eblocks)]
    single_eblocks['plate_id'] = single_eblocks['position'].apply(lambda x: int(x.split('_')[0]) + n_plates)
    single_eblocks['Well Position'] = single_eblocks['position'].apply(lambda x: x.split('_')[1])

    n_designs += int(len(single_eblocks))
    n_plates = int(single_eblocks['plate_id'].max())

# Re-format fragment eBlocks.
if len(frag_eblocks) > 0:
    unique_idx = np.unique([v[:-1] for v in frag_eblocks['readin_order'].values])
    oldidx2newidx = {old:n_designs + i for i, old in enumerate(unique_idx)}
    frag_eblocks['design_id'] = frag_eblocks['readin_order'].apply(lambda x: args.design_prefix + str(oldidx2newidx[x[:-1]]).zfill(4) + x[-1])
    frag_eblocks['position'] = w96_2blocks[:len(frag_eblocks)]
    frag_eblocks['plate_id'] = frag_eblocks['position'].apply(lambda x: int(x.split('_')[0]) + n_plates)
    frag_eblocks['Well Position'] = frag_eblocks['position'].apply(lambda x: x.split('_')[1])

    n_designs += int(len(frag_eblocks) / 2)
    n_plates = int(frag_eblocks['plate_id'].max())

# Reformat heterooligomer eBlocks.
if len(hto_eblocks) > 0:

    # Check oligomeric stoichiometries.
    chains = [dn.split('_')[-2] for dn in hto_eblocks['design_name'].values] # get all chains
    n_each_chain = [chains.count(ch) for ch in np.unique(chains)]
    if len(np.unique(n_each_chain)) > 1:
        print('No plate formatting currently implemented for mixed stoichiometry heterooligomers. Chains have different IDs and the layout is the same as for monomers.')
        design_ids = [args.design_prefix + str(i).zfill(4) for i in range(n_designs, n_designs + len(hto_eblocks))]
        hto_eblocks['design_id'] = design_ids
        hto_eblocks['position'] = w96_zigzag[:len(hto_eblocks)]
        hto_eblocks['plate_id'] = hto_eblocks['position'].apply(lambda x: int(x.split('_')[0]) + n_plates)
        hto_eblocks['Well Position'] = hto_eblocks['position'].apply(lambda x: x.split('_')[1])

    else:
        n_unique_ch = len(np.unique(chains))
        ch_idx = []
        for ch in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'[:n_unique_ch]:
            ch_idx += [args.design_prefix + str(i).zfill(4) + f'_{ch}' for i in range(n_designs, n_designs+int(len(hto_eblocks)/n_unique_ch))]

        hto_eblocks['design_id'] = np.hstack(np.reshape(ch_idx, (n_unique_ch, -1)).T)

        if n_unique_ch > 4:
            print('No plate formatting for heterooligomers beyond tetramers currently implemented. The layout is the same as for monomers.')
            hto_eblocks['position'] = w96_zigzag[:len(hto_eblocks)]

        else:
            if n_unique_ch == 2:
                hto_eblocks['position'] = w96_2blocks[:len(hto_eblocks)]
            elif n_unique_ch == 3:
                hto_eblocks['position'] = w96_3blocks[:len(hto_eblocks)]
            elif n_unique_ch == 4:
                hto_eblocks['position'] = w96_4blocks[:len(hto_eblocks)]

    hto_eblocks['plate_id'] = hto_eblocks['position'].apply(lambda x: int(x.split('_')[0]) + n_plates)
    hto_eblocks['Well Position'] = hto_eblocks['position'].apply(lambda x: x.split('_')[1])


    n_designs += int(len(hto_eblocks) / 2)
    n_plates = int(hto_eblocks['plate_id'].max())

# Reformat heterooligomer eBlocks that have been fragmented.
if len(htd_frag_eblocks) > 0:
    print('No plate formatting currently implemented for fragmented heterooligomers. Chains have different IDs and the layout is the same as for fragmented monomers.')
    unique_idx = np.unique([v[:-1] for v in htd_frag_eblocks['readin_order'].values])
    oldidx2newidx = {old:n_designs + i for i, old in enumerate(unique_idx)}
    htd_frag_eblocks['design_id'] = htd_frag_eblocks['readin_order'].apply(lambda x: args.design_prefix + str(oldidx2newidx[x[:-1]]).zfill(4) + x[-1])
    htd_frag_eblocks['position'] = w96_2blocks[:len(htd_frag_eblocks)]
    htd_frag_eblocks['plate_id'] = htd_frag_eblocks['position'].apply(lambda x: int(x.split('_')[0]) + n_plates)
    htd_frag_eblocks['Well Position'] = htd_frag_eblocks['position'].apply(lambda x: x.split('_')[1])

    n_designs += int(len(htd_frag_eblocks) / 2)
    n_plates = int(htd_frag_eblocks['plate_id'].max())


reformated_df = pd.concat([single_eblocks, frag_eblocks, hto_eblocks, htd_frag_eblocks])

if args.no_layout:
    if args.echo:

        def gen_design_id(r):
            if ('x' in r.readin_order) or ('y' in r.readin_order):
                return args.design_prefix + str(int(args.design_id) + int(r.readin_order[:-1])).zfill(4) + r.readin_order[-1]

            else:
                return args.design_prefix + str(int(args.design_id) + int(r.readin_order)).zfill(4)

        df['design_id'] = df.apply(gen_design_id, axis=1)
        df['position'] = w384[:len(df)]
        df['plate_id'] = df['position'].apply(lambda x: int(x.split('_')[0]))
        df['Well Position'] = df['position'].apply(lambda x: x.split('_')[1])

        reformated_df = df.copy()

        # Generate ECHO transfer file -- currently only to 1frag/1vector type of cloning.
        df['Source Plate Name'] = ''
        df['Source Well'] = df['Well Position']
        for p in df['plate_id'].unique():
            df.loc[df['plate_id']==p, 'Source Plate Name'] = f'{date}_{args.order_name}.{int(p)}_{args.species}_{enzyme}'

        # Assumes: 348w (eBlocks) --> 96w (cultures)
        split_idx = [df.index.values[x:x+96] for x in range(0, len(df), 96)]
        df['Destination Plate Name'] = ''
        df['Destination Well'] = ''
        for i, si in enumerate(split_idx):
            df.loc[si, 'Destination Plate Name'] = f'PCR96 - {i+1}'
            df.loc[si, 'Destination Well'] = w96_zigzag_single[:len(si)]

        # GG mastermix: hard-coded defintions.
        tot_rxn_vol = 1 # uL
        eblock_conc = 4 # ng/uL -- as delivered by IDT
        vector_conc = 100 # ng/uL -- assumes vector stocks have been normalised
        mol_vector = 4 # fmol, target amount
        insert_to_vector = 2 # target insert:vector molar ratio

        eblock_bp = df['length_eblock'].median()
        vector_bp = len([str(r.seq) for r in SeqIO.parse(glob.glob(f'./entry_vectors/{args.gg_vector[0]}*.fa')[0], 'fasta')][0]) # only considers first vector if multiple were specified
        vector_vol = mol_vector / (((vector_conc * 0.000000001) / ((vector_bp * 617.96) + 36.04)) * 1000000000000000)
        eblock_vol = 0.025 * np.ceil(((mol_vector * insert_to_vector) / (((eblock_conc * 0.000000001) / ((eblock_bp * 617.96) + 36.04)) * 1000000000000000)) / 0.025)
        mm_vol = tot_rxn_vol - eblock_vol

        echo_df = df[['Source Plate Name', 'Source Well', 'Destination Plate Name', 'Destination Well']].copy()
        echo_df['Transfer Volume'] = eblock_vol * 1e3 # transfer volumes are specified in nL
        echo_df2 = echo_df.copy() # for specifying GGMM destinations
        echo_df2['Source Well'] = 'X' # replace manually later
        echo_df2.loc[:, 'Transfer Volume'] = mm_vol * 1e3

        echo_df = pd.concat([echo_df, echo_df2])
        echo_df.to_csv(f'{filename}_ECHO.csv', index=False)


    else:
        design_ids = [args.design_prefix + str(i).zfill(4) for i in range(int(args.design_id), int(args.design_id) + len(df))]
        df['design_id'] = design_ids
        df['position'] = w96[:len(df)]
        df['plate_id'] = df['position'].apply(lambda x: int(x.split('_')[0]))
        df['Well Position'] = df['position'].apply(lambda x: x.split('_')[1])
        reformated_df = df.copy()

reformated_df = reformated_df.astype({'plate_id':'int32'})

# Change some characters in design names otherwise IDT complains.
change_char = {
'_isheterooligomer':'',
'+':'plus'
}
def name(row):

    clean_name = row['design_name']
    for k, v in change_char.items():
        clean_name = clean_name.replace(k, v)

    return row['design_id'] + '__' + row['Well Position'] + '__' + row['order_name'] + '.' + str(int(row['plate_id'])) + '__' + clean_name

reformated_df['Name'] = reformated_df.apply(name, axis=1)

# ============================================
# OUTPUTS
# ============================================
# Save dataframe as log.
reformated_df.to_csv(f'{filename}_df.csv', index=False)

# Save .xlsx for submitting the order to IDT (one sheet per plate).
with pd.ExcelWriter(f'{filename}.xlsx') as writer:

    for p in reformated_df['plate_id'].unique():

        sheet_df = reformated_df[reformated_df['plate_id']==p][['Well Position', 'Name', 'Sequence']]

        if args.echo: # need all wells up to P24 for IDT to understand it as a 384w plate.
            sheet_df = sheet_df.merge(w384_df, left_on='Well Position', right_on='Well Position', how='right')

        sheet_df.to_excel(writer, sheet_name=f'{date}_{args.order_name}.{int(p)}_{args.species}_{enzyme}', index=False)

# Save FASTA file(s) of expression products for protparam.
for gg_v in args.gg_vector:
    with open(f'{filename}_into_{gg_v}.fa', 'w') as f:
        for i, r in reformated_df.iterrows():
            if 'y' in r['design_id']:
                pass

            else:
                if 'x' in r['design_id']:
                    f.write('>' + r['design_id'].replace('x', '') + '__' + '__'.join(r['Name'].split('__')[1:]) + '\n')
                    f.write(r[f'exp_aa_seq_from_{gg_v}'] + '\n')

                else:
                    f.write('>' + r['Name'] +  '\n')
                    f.write(r[f'exp_aa_seq_from_{gg_v}'] + '\n')


# Save text file with plate layouts and cloning instructions.
txt = f'CLONING INSTRUCTIONS FOR {date}_{args.order_name}_{args.species}_{enzyme}:\n'
single_frag = reformated_df[reformated_df['is_frag']==False]
dual_frag = reformated_df[reformated_df['is_frag']==True]

if len(single_frag) > 0:
    txt += f' * Number of 1-fragment reactions: {len(single_frag)}\n'

if len(dual_frag) > 0:
    txt += f' * Number of 2-fragments reactions: {int(len(dual_frag)/2)}\n'

txt += '\n'

# Layouts.
w96 = np.reshape([r + str(c) for r in 'ABCDEFGH' for c in range(1,13)], (8, 12))
idx2row = {i:r for i, r in enumerate('ABCDEFGH')}
for p in reformated_df['plate_id'].unique():
    plate_df = reformated_df[reformated_df['plate_id']==p]
    txt += f'{date}_{args.order_name}.{int(p)}_{args.species}_{enzyme}\n'
    txt += '  ' + ''.join([str(c).center(11) for c in np.arange(1,13)]) + '\n'
    for i, row in enumerate(w96):
        txt += idx2row[i] + ' |'

        for j, well in enumerate(row):

            if well in plate_df['Well Position'].values:
                well_content = plate_df[plate_df['Well Position']==well]['design_id'].values[0].center(10)

            else:
                well_content = ''.center(10)

            txt += well_content + '|'
        txt += '\n'
    txt += '\n\n'

with open(f'{filename}_cloning_instructions.txt', 'w') as f:
    f.write(txt)
    
    
# Unless supressed, generate cloned plasmid maps.
if not args.no_plasmids:
    os.makedirs(output_folder + 'cloned_plasmids/', exist_ok=True) # make a subfolder to store the plasmid maps

    for gg_v in args.gg_vector:
        
        for i, r in reformated_df.iterrows():
            if 'y' in r['design_id']:
                pass
        
            else:
                plasmid_name = r['design_id'].replace('x', '') + '_into_' + gg_v
                
                with open(output_folder + 'cloned_plasmids/' + plasmid_name + '.fa' , 'w') as f:
                    f.write('>' + plasmid_name + '\n')
                    f.write(r[f'{gg_v}_cloned_plasmid_seq'] + '\n')


# Make symlink of PDBs with design_id + well + plate appended to the filename.
if args.order_pdbs != None:
    print('Making PDB symlinks...')

    for i, r in reformated_df.iterrows():

        if 'y' in r['design_id']: # don't make symlinks for y-fragments
            pass

        else:
            info_label = '__'.join(r['Name'].split('__')[:3]).replace('x', '')

            if '_isheterooligomer' in r['design_name']:
                pdb_name = '_'.join(r['design_name'].split('_')[:-2]) + '.pdb'

            else:
                pdb_name = r['design_name'] + '.pdb'

            if args.order_pdbs + pdb_name in glob.glob(f'{args.order_pdbs}*.pdb'):
                os.system(f'ln -s {pdb_name} {args.order_pdbs}{info_label}__{pdb_name}') # create symlink of PDB with protein ID, well ID and plate ID

if len(skipped) > 0:
    print(f'The following {len(skipped)} designs could not be converted to eBlocks and were skipped:')
    with open(f'{filename}_FAILED.fa', 'w') as f:
        for k, v in skipped.items():
            f.write(f'>{k}\n{v}\n')
            print(f'>{k}\n{v}')

else:
    print('All designs were successfully converted to eBlocks.')

print(f"Order contains {len(reformated_df)} eBlocks across {len(reformated_df['plate_id'].unique())} plate(s).")
print(f"Place order at https://www.idtdna.com/site/order/plate/eblocks\n")

if args.verbose:
    print(
    "John Bercow says:\n"
    "////////////////////////////////////////////+///////////////////////////////////////////////////////\n"
    "////////////////////////////////////+shdmNNNNmmhs+://///////////////////////////////////////////////\n"
    "/////////////////////////////////+ymNNMMMMMNNhssso+so+//////////////////////////////////////////////\n"
    "///////////////////////////////ohmNNNNMMMMMNhs++::++sso+////////////////////////////////////////////\n"
    "////////////////////////////+ymNNMNNNMMMMMNmhs+///:::+os+///////////////////////////////////////////\n"
    "//////////////////////////+ymNMMMNMMMMMMMNmhyo//::/++yysy///////////////////////////////////////////\n"
    "/////////////////////////shmMMMMMMMMMMMMMNdhyyoosshhhddhy+://///////////////////////////////////////\n"
    "/////////////////////////+hmNNNNNMMMNmdmddhyyyydNdNNNhs+/--/////////////////////////////////////////\n"
    "///////////////////////////oyhddhdddhhysyyyhhdNMMMMNmho:-../////////////////////////////////////////\n"
    "//////////////////////////:+yssyyyddmNNNmNmmmmNMMMMmhhys+..:////////////////////////////////////////\n"
    "///////////////////////////shhhhddmmNNNNNNNNNNmmNNNmhyys+:-/////////////////////////////////////////\n"
    "//////////////////////////+hdNddmmNNNNNNNNNNmmhysso+osys+:-/////////////////////////////////////////\n"
    "/////////////////////////:`./hmNNNNNNNNNNNmmds/:-----:/o-:://///////////////////////////////////////\n"
    "//////////////////////////.``.os//oymNNmNmmmdy:.``-:-/o+`.://///////////////////////////////////////\n"
    "//////////////////////////-.-`+h+` `-+oyyydmdy:..++:+/yy..://///////////////////////////////////////\n"
    "//////////////////////////:::/NNds+:-.:+odmdhysoooyshdd:`-//////////////////////////////////////////\n"
    "//////////////////////////:/sdMNmhyhhhyhhdhhhdhyooydmd/: `--:///////////////////////////////////////\n"
    "///////////////////////////-+dmmmdosyyysssyhhhy+://::..y+    `.--://////////////////////////////////\n"
    "/////////////////////////::..----:yyyyssooooo+/---:.``/NN`        ``..--:///////////////////////////\n"
    "//////////////////////:-.`   .:/sdmhhhys+/::-........-dMM:              ``.:////////////////////////\n"
    "//////////////////::-``       ``.-++syso+/-.....--..-hMMM-                 ``.-::///////////////////\n"
    "///////////////::.`           .`  ./o+oo+/---:::..-+dMMMN`                      `.-:////////////////\n"
    "/////////////:-`              :/-:+/+////::---.`-+dNMMMMs                          `:///////////////\n"
    "///////////:.`                `-:-``.-::-..`   `-+osydmm-                           `-//////////////\n"
    "/////////:-                        `...````   .---``+yyo                              -/////////////\n"
    "////////:`                               ..`  .----sNMMo                               -////////////\n"
    "///////:                                `+o/:.+y`:hMMMM-                                :///////////\n"
    "///////.                                 -mms/`/+yMMMMm                                 `://////////\n"
    "//////:                                   +Mm:./h/dMMM+                                  `://///////\n"
    "/////:`                                    om.h:://hNN.                                   `:////////\n"
    "/////.                                      .:--:ho+ho                                     `:///////\n"
    "////:                                        ..//s+/s`                                      .///////\n"
    "////`                                         `o+:///                                        -//////\n"
    "///:                                           -::sh`                                        `://///\n"
    "///.                                            :o+:                                          -/////\n"
    "//:`                                            -s/`                                           -////\n"
    "//-                                              ``                                            `:///\n"
    "//-                                                                                             -///\n"
    "//.                                                                                             `://\n"
    "/:`                                             ss/.`                                            `:/\n"
    "/-                                             `NMMMmy+:`                                         ./\n"
    "/.                          ``                  oMMMMMMMNdy/.          ``                          :\n"
    ":`                           +:                  hMMMMMMMMMMMmy+-`  `+dNMms/.                      .\n"
    "-                            .d:-/+oys           .NMMMMMMMMMMMMMMMdshNMMMMMMmh+.                   `\n"
    "`                        `:osshdNNNh/.            /MMMMMMMMMMMMMMMMMMMMMNNMNmNNm+`                  \n"
    "                        -sNNNmNmNNmho+-            sMMMMMMMMMMMMMMMMMMMMmhmmNNNNN+                  \n"
    "                       .++hhyo+::::/.``            `dMMMMMMMMMMMMMMMMMMmshhddyNNd-                  \n"
    "``````.....````````....-:+hhhys+..+y-`...........```:dddddmmdmmddddhyhds-:/soydy/.......``````......\n"
    "://+ydmNNNmdyo///:/mmmmmmmNNNNmmmho//:dmmmmmmmmmmmdyo/////mNNNNNNNNNmNNmy/+mmNNmmmmmmmmmdy+///mmmmmy\n"
    "/+dNMMMMMMMMMMmo///MMMMMMMMMMMMMMMMy//NMMMMMMMMMMMMMMms///MMMMMMMMMMMMMMh/oMMMMMMMMMMMMMMMNs//MMMMMh\n"
    "+NMMMMMNmmNMMMMNs//MMMMMdsssssmMMMMM//NMMMMmsssyhmMMMMMs//MMMMMdsssssssso/oMMMMMhsssssNMMMMm//NMMMMh\n"
    "dMMMMNs///omMMMMN//MMMMMNdddddNMMMMm//NMMMMd//////dMMMMN//MMMMMNddddddddy/oMMMMMmdddddNMMMMd//NMMMMy\n"
    "NMMMMd/////yMMMMM+/MMMMMMMMMMMMMMNd+//NMMMMd//////yMMMMM+/MMMMMMMMMMMMMMh/oMMMMMMMMMMMMMMNh///dMMMMo\n"
    "hMMMMMy+/+sNMMMMm//MMMMMmyymMMMMNo////NMMMMd/////oNMMMMm//MMMMMdsssssssss/oMMMMMdyhNMMMMN+////sNNNM/\n"
    "/mMMMMMMNMMMMMMNo//MMMMMy///mMMMMN+///NMMMMNhhhhmMMMMMN+//MMMMMdsssssssss/oMMMMMs//+mMMMMm+/////////\n"
    "//yNMMMMMMMMMMh+///MMMMMy////dMMMMNo//NMMMMMMMMMMMMMNy////MMMMMMMMMMMMMMh/oMMMMMs///+mMMMMN+//MMMMMh\n"
    "////oydmmmdhs//////ddddds/////hddddd+/hddddddddddhyo//////ddddddddddddddy/+dddddo/////hddddh//ddddds"
    )
