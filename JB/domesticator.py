#!/software/containers/john_bercow.sif

# Ryan's domesticator with some minor modifications.

# ============================================
# LIBRARIES
# ============================================
import copy
import numpy as np
from collections import Counter
from Bio import SeqFeature

# Domesticator uses dnachisel
import dnachisel
from dnachisel import DnaOptimizationProblem, NoSolutionError
from dnachisel import DEFAULT_SPECIFICATIONS_DICT
from dnachisel import Location
from dnachisel import Specification, SpecEvaluation

# Special k-mer minimasation objective from Ryan's domesticator.
class MinimizeNumKmers(Specification):
    """Minimizes a kmer score."""

    best_possible_score = 0

    def __init__(self, k=8, location=None, boost=1.0):
        self.location = location
        self.k = k
        self.boost = boost

    def initialize_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)

    def evaluate(self, problem):
        """Return a customized kmer score for the problem's sequence"""
        sequence = self.location.extract_sequence(problem.sequence)
        all_kmers = [sequence[i : i + self.k] for i in range(len(sequence) - self.k)]
        number_of_non_unique_kmers = sum(
            [count for kmer, count in Counter(all_kmers).items() if count > 1]
        )
        score = -(float(self.k) * number_of_non_unique_kmers) / len(sequence)
        return SpecEvaluation(
            self,
            problem,
            score=score,
            locations=[self.location],
            message="Score: %.02f (%d non-unique %d-mers)"
            % (score, number_of_non_unique_kmers, self.k),
        )

    def label_parameters(self):
        return [("k", str(self.k))]

    def short_label(self):
        return f"Avoid {self.k}mers {self.boost}"

    def __str__(self):
        """String representation."""
        return "MinimizeNum%dmers" % self.k

DEFAULT_SPECIFICATIONS_DICT["MinimizeNumKmers"] = MinimizeNumKmers

# ============================================
# FUNCTIONS
# ============================================
def reverse_translate(
        amino_acid_sequence,
        kmers_weight=1.0,
        cai_weight=1.0,
        hairpins_weight=1.0,
        max_tries=10,
        species='e_coli',
        avoid=['GGTCTC', 'GAGACC']
    ):
    '''
    Ryan's domesticator.
    '''

    # Generate a naively reverse-translated DNA sequence first.
    naive_dna_sequence = dnachisel.reverse_translate(amino_acid_sequence)

    # Sequence optimisation will happen across the whole sequence.
    location = Location.from_biopython_location(SeqFeature.FeatureLocation(0, len(amino_acid_sequence) * 3))

    # Add optimisation objectives.
    objectives = []
    objectives.append(MinimizeNumKmers(k=8, boost=kmers_weight, location=location))
    objectives.append(dnachisel.builtin_specifications.AvoidHairpins(boost=hairpins_weight, location=location))
    objectives.append(dnachisel.builtin_specifications.MaximizeCAI(species=species, boost=cai_weight, location=location))

    # Add optimisation constraints.
    constraints = []
    constraints.append(dnachisel.builtin_specifications.EnforceTranslation(location=location))

    # A series of sequence patterns to remove
    constraints.append(dnachisel.builtin_specifications.AvoidPattern("AAAAA", location=location)) # terminator
    constraints.append(dnachisel.builtin_specifications.AvoidPattern("TTTTT", location=location)) # terminator
    constraints.append(dnachisel.builtin_specifications.AvoidPattern("CCCCCC", location=location)) # repetitions
    constraints.append(dnachisel.builtin_specifications.AvoidPattern("GGGGGG", location=location)) # repetitions
    constraints.append(dnachisel.builtin_specifications.AvoidPattern("ATCTGTT", location=location)) # T7/T3 RNA polymerase pausing
    constraints.append(dnachisel.builtin_specifications.AvoidPattern("GGRGGT", location=location)) # G-quadruplex?

    for seq in avoid: # GG enzyme recognition site.
        constraints.append(dnachisel.builtin_specifications.AvoidPattern(seq, location=location))

    # Organism-specific constraints
    if species == 'e_coli':
        constraints.append(dnachisel.builtin_specifications.AvoidPattern("GGAGG", location=location)) # E. coli Shine-Dalgarno
        constraints.append(dnachisel.builtin_specifications.AvoidPattern("TAAGGAG", location=location)) # strong RBS
        constraints.append(dnachisel.builtin_specifications.AvoidPattern("GCTGGTGG", location=location)) # Chi site in E. coli

        # Alternative start sites: G/A rich 6 nt upstream of ATG or GTG or TTG (cryptic start sites)
        # N = A or T or C or G
        # D = A or G or T
        # R = A or G
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC523762/?page=1
        # This constraint can be too harsh and lead to problems with some sequences.
        constraints_easier = copy.deepcopy(constraints)
        constraints.append(dnachisel.builtin_specifications.AvoidPattern("RRRRRNNNNNDTG", location=location)) # 5 nt spacing
        constraints.append(dnachisel.builtin_specifications.AvoidPattern("RRRRRNNNNNNDTG", location=location)) # 6 nt spacing
        constraints.append(dnachisel.builtin_specifications.AvoidPattern("RRRRRNNNNNNNDTG", location=location)) # 7 nt spacing

    elif species == 'h_sapiens':
        constraints.append(dnachisel.builtin_specifications.AvoidPattern("GCCRCCATGG", location=location)) # Kozak sequence

    # Global GC content (according to TWIST).
    constraints.append(dnachisel.builtin_specifications.EnforceGCContent(mini=0.25, maxi=0.65, location=location))

    # Local GC content (according to TWIST).
    constraints.append(dnachisel.builtin_specifications.EnforceGCContent(mini=0.35, maxi=0.65, window=50, location=location))

    # Start optimisation.
    solutions = []
    solution_found = False
    for i in range(max_tries):

        if solution_found:
            break

        try:
            if species == 'e_coli' and i >= max_tries/2:
                print('  [!] Preventing alternative start sites removed from the list of optimisation constraints.')
                initial_problem = DnaOptimizationProblem(naive_dna_sequence, constraints=constraints_easier, objectives=objectives, logger=None)

            else:
                initial_problem = DnaOptimizationProblem(naive_dna_sequence, constraints=constraints, objectives=objectives, logger=None)

            problem = copy.deepcopy(initial_problem)
            problem.resolve_constraints_by_random_mutations()
            problem.optimize()
            problem.resolve_constraints(final_check=True)
            solutions.append(problem)
            solution_found = True

        except NoSolutionError:
            initial_problem.max_random_iters += 1000
            solution_found = False

            continue

    if len(solutions) == 0:
        raise NoSolutionError(f"No solution found for {amino_acid_sequence}", initial_problem)

    # Return the best solution according to dnachisel if multiple attemptes were made.
    scores = [solution.objectives_evaluations().scores_sum() for solution in solutions]
    best_idx = np.argmin(scores)

    best_solution = solutions[best_idx]

    return best_solution.sequence
