# John Bercow for DNA ordering

JohnBercow script for ordering gene fragments for SAPP.

Simply pass a folder of PDBs, or FASTA file (or both) and the script will generate a `.xlsx` file of reverse-translated sequences with the correct overhangs for Golden Gate Assembly into SAPP vectors ready to be uploaded to commerical DNA suppliers.

Explanation for the different scripts:
- `JohnBercow.py`: the main script allowing you to specify the different options.
- `domesticator.py`: the module for DNA optimization.
- `idt.py`: the module with the IDT API for querying synthesizability on the fly. Use `--skip_idt_query` if you don't want to use that option (e.g. using a supplier different than IDT).
- `entry_vectors`: contains the `.fa` files of the SAPP vectors in which cloning can be automated.

```
usage: JohnBercow.py [-h] [--order_pdbs ORDER_PDBS] [--order_fasta ORDER_FASTA] --order_name ORDER_NAME --gg_vector GG_VECTOR [GG_VECTOR ...] --species SPECIES --design_prefix DESIGN_PREFIX --design_id DESIGN_ID [--skip_idt_query] [--idt_score IDT_SCORE]
                     [--starting_kmers_weight STARTING_KMERS_WEIGHT] [--n_domesticator_steps N_DOMESTICATOR_STEPS] [--max_attempts MAX_ATTEMPTS] [--max_length MAX_LENGTH] [--print_hto] [--no_layout] [--no_plasmids] [--verbose] [--echo]

 * Generates an IDT-ready .xlsx file for ordering eBlocks from a folder of PDBs and/or a concatenated FASTA file.
 * Appropriate overhangs for Golden Gate cloning into entry vector(s) of interest are added automatically.
 * RECOMMENDED: check your GG assemblies at https://goldengate.neb.com/#!/
 * Wondering why the script is called John Bercow? https://www.youtube.com/watch?v=VYycQTm2HrM&ab_channel=TheSun

 * EXAMPLE AVAILABLE ENTRY VECTORS:
 *** see ./entry_vectors/ for the FULL list ***
 - LM0627 (BsaI): C-term SNAC-His | MSG[...]GSGSHHWGSTHHHHHH
 - LM0670 (BsaI): C-term His | MSG[...]GSHHHHHH
 - LM0671 (BsaI): N-term sfGFP and C-term His | MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFARYPDHMKQHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYKGSSG[...]GSHHHHHH
 - LM0673 (BsaI): N-term sfGFP-SNAC and C-term His | MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFARYPDHMKQHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYKGGSHHWSSG[...]GSHHHHHH
 - LM1371 (BsaI): N-term His | MSHHHHHHSG[...]GS
 - LM1425 (BsaI): C-term AviTag-His | MSG[...]GSGLNDIFEAQKIEWHESHHHHHH
 - BW1001 (BsaI): C-term smBiT-His | MSG[...]GSSGSGGSGGGGSGGSSSGGVTGYRLFEEILGSHHHHHH
 - BW1002 (BsaI): C-term smBiT-GB1-His | MSG[...]GSSGSGGSGGGGSGGSSSGGVTGYRLFEEILGSTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTEGSHHHHHH
 - BW1003 (BsaI): C-term lgBiT | MSG[...]GSSGSGGSGGGGSGGSSSGGVFTLEDFVGDWEQTAAYNLDQVLEQGGVSSLLQNLAVSVTPIQRIVRSGENALKIDIHVIIPYEGLSADQMAQIEEVFKVVYPVDDHHFKVILPYGTLVIDGVTPNMLNYFGRPYEGIAVFDGKKITVTGTLWNGNKIIDERLITPDGSMLFRVTINSGSHHHHHH
 - BW1004 (BsaI): N-term His-smBiT | MSGHHHHHHGSVTGYRLFEEILGGSGSGGSGGGGSGGSSSGG[...]GS
 - BW1005 (BsaI): N-term His-GB1-smBiT | MSGHHHHHHGSTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTEGSVTGYRLFEEILGGSGSGGSGGGGSGGSSSGG[...]GS
 - BW1006 (BsaI): N-term His-lgBiT | MSGHHHHHHGSVFTLEDFVGDWEQTAAYNLDQVLEQGGVSSLLQNLAVSVTPIQRIVRSGENALKIDIHVIIPYEGLSADQMAQIEEVFKVVYPVDDHHFKVILPYGTLVIDGVTPNMLNYFGRPYEGIAVFDGKKITVTGTLWNGNKIIDERLITPDGSMLFRVTINSGGSGSGGSGGGGSGGSSSGG[...]GS

options:
  -h, --help            show this help message and exit
  --order_pdbs ORDER_PDBS
                        path to a folder containing the PDBs you want to order. Either this option or --order_fasta needs to be specified. Both can specified.
  --order_fasta ORDER_FASTA
                        path to a FASTA file containing the AA sequences you want to order. Either this option or --order_pdbs needs to be specified. Both can be specified.
  --order_name ORDER_NAME
                        name of the order (appended to output files). NB date, plate IDs, organism, and enzyme are added automatically.
  --gg_vector GG_VECTOR [GG_VECTOR ...]
                        name of target vector(s) for Golden Gate cloning (determines the DNA adapters). Also determines the AA tags appended to the design in the FASTA output. Multiple vectors can be specified as space-separated values, but they need to be compatible (i.e. use the same enzyme and same overhangs).
  --species SPECIES     codon optimisation will be performed for this species (e.g. e_coli, s_cerevisiae, h_sapiens, etc...)
  --design_prefix DESIGN_PREFIX
                        designs get IDs with this prefix (e.g. LM0001, LM0002, etc...)
  --design_id DESIGN_ID
                        increment design indices from this number.
  --skip_idt_query      skip IDT website query that checks synthesisability of eBlocks.
  --idt_score IDT_SCORE
                        IDT complexity score threshold for accepting a reverse translated sequence (Default: 7). 0-7 is green, 7-15 is yellow, >15 is red and IDT will not make it.
  --starting_kmers_weight STARTING_KMERS_WEIGHT
                        starting value for the kmers_weight setting of Domesticator (Default: 10). This parameter is linearly ramped (up to 100) over --n_domesticator_steps.
  --n_domesticator_steps N_DOMESTICATOR_STEPS
                        maximum number of Domesticator steps attempted (Default: 10). The kmers_weight parameter (which increases synthesiability of repetitive sequences) is linearly ramped up to 100 over this number of steps.
  --max_attempts MAX_ATTEMPTS
                        maximum number of reverse translation attempts at each Domesticator step (Default: 20). Since Domesticator is stochastic, re-running the optimisation problem with the same parameters can lead to different solutions.
  --max_length MAX_LENGTH
                        maximum length of eBlocks. IDT's maximum is 1500 bp, but less can be specified if sequence complexity is a issue for synthesis and you want to force the generation of smaller fragments.
  --print_hto           print heterooligomers to stdout.
  --no_layout           do not apply automated layout formatting.
  --no_plasmids         do not generate the cloned plasmid maps.
  --verbose             increase the verbosity of th e output (recommended).
  --echo                generates outputs formated as 384w plates (for ordering into ECHO-qualified plates). Currently only available in conjuction with the no_layout option.
```

---
# Note on IDT credentials

Proper script execution requires IDT credentials (for the IDT API queries). These should be located at `~/idt_credentials`. Details on how to obtain these credentials can be found at https://eu.idtdna.com/pages/tools/apidoc

If you don't have IDT credentials and still want to use the script, use option `--skip_idt_query`



