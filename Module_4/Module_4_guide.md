"""
🧬 Module 4: Probabilistic Optimization and Stochastic Search
Concept: Moving beyond deterministic limits to handle biological noise and local optima.

This module addresses the "Zero-Probability" problem in DNA profiling and 
introduces randomized searching to ensure we find the globally optimal 
motifs across a genomic landscape.
"""

# 🏗️ Phase 4.1: Greedy Motif Search with Pseudocounts
# Concept: Implementing Laplace's Rule of Succession to stabilize the Profile Matrix.
#
# In standard counting, if a specific nucleotide is absent from a position, its 
# probability becomes 0. Since we multiply probabilities, one missing base 
# "kills" the entire score. Pseudocounts add a small buffer (+1) to every base 
# to ensure no biological signal is ignored due to a small sample size.

def ProfileWithPseudocounts(motifs):
    """
    🛠️ Implementation
    Calculates the profile matrix by applying a +1 pseudocount to all 
    nucleotide counts before normalization. This algorithm incorporates 
    a stabilized profile to maintain sensitivity in complex datasets.
    """
    t = len(motifs)
    k = len(motifs[0])
    # Initialize with 1 to prevent zero-probability multiplications
    profile = {nt: [1] * k for nt in "ACGT"}
    
    for motif in motifs:
        for i, nt in enumerate(motif):
            profile[nt][i] += 1
            
    for nt in "ACGT":
        # Normalize using (t + 4) to account for the four added pseudocounts
        profile[nt] = [x / (t + 4) for x in profile[nt]]
    return profile

def GreedyMotifSearchWithPseudocounts(DNA, k, t):
    """
    Iteratively refines motifs using pseudocount-stabilized profiles.
    """
    n = len(DNA[0])
    best_motifs = [DNA[i][:k] for i in range(t)]
    
    for i in range(n - k + 1):
        motifs = [DNA[0][i:i + k]]
        for j in range(1, t):
            # Use the stabilized profile from previously selected motifs
            profile = ProfileWithPseudocounts(motifs)
            motifs.append(ProfileMostProbableKmer(DNA[j], k, profile))
            
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs

    return best_motifs

"""
🏗️ Phase 4.2: Randomized Motif Search
Concept: Utilizing stochastic exploration to bypass local optima.

Deterministic algorithms often get trapped in "local minima"—they find a motif 
that looks acceptable within a small range but miss the global optimum. 
Randomized Motif Search starts with random k-mer coordinates and iteratively 
refines them. Running this random process multiple times allows for a 
comprehensive exploration of the genomic landscape.
"""

import random

def RandomMotifs(DNA, k):
    """
    🛠️ Implementation
    Randomized Selection: Unlike deterministic algorithms, this starts by 
    randomly selecting k-mers (substrings of length k) from each sequence.
    """
    motifs = []
    for sequence in DNA:
        # Select a random starting position for each DNA strand
        start = random.randint(0, len(sequence) - k)
        motifs.append(sequence[start:start + k])
    return motifs

def RandomizedMotifSearch(DNA, k, t):
    """
    Performs stochastic motif discovery through iterative refinement.
    
    BIOLOGICAL SIGNIFICANCE:
    ------------------------
    - This algorithm 'descends' toward the lowest mismatch score and stops 
      when it reaches a local peak. 
    - The randomness helps explore more possible motif combinations, 
      potentially avoiding local optima.
    """
    # Step 1: Initialize motifs randomly
    motifs = RandomMotifs(DNA, k)
    best_motifs = motifs

    while True:
        # Step 2: Build a profile from the current selection with pseudocounts
        profile = ProfileWithPseudocounts(motifs)
        
        # Step 3: Find the most probable k-mer in each sequence based on current profile
        motifs = [ProfileMostProbableKmer(seq, k, profile) for seq in DNA]
        
        # Step 4: Update best_motifs if the new set improves the score
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs
        else:
            # Step 5: Stop iterating when no further improvement is observed
            return best_motifs

# Example Execution Setup
dna_input = [
    "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
    "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
    "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
    "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
    "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
]

""" 
🏁 REPOSITORY MILESTONE: PROJECT COMPLETION

The algorithmic architecture for this genomic analysis toolkit is now complete. 
From initial sequence analysis to stochastic motif discovery, the 
'Concrete Slurry' has set, and the 'Final Touches' are applied.

🌟 PROJECT SUMMARY:
- Module 1: Foundational Genome Replication (Skew, Forward/Reverse Strands)
- Module 2: Pattern Identification (DnaA Boxes & Hidden Messages)
- Module 3: Comparative Alignment (Greedy Search & Consensus Discovery)
- Module 4: Probabilistic Refinement (Pseudocounts & Randomized Search)

🏆 BIOLOGICAL IMPACT:
These scripts provide the computational framework necessary to decode 
regulatory regions in circular bacterial genomes. By combining deterministic 
logic with probabilistic stability, we can now identify conserved signals 
that drive life at the molecular level.

THE END OF THE INITIAL JOURNEY. CODE COMMITTED. KNOWLEDGE ARCHIVED.
"""

def ProjectStatus():
    """ Returns the final status of the Bioinformatics journey. """
    modules_completed = 4
    repository_status = "STABLE"
    optimization_level = "MAXIMUM"
    return f"Journey complete across {modules_completed} Modules. Repo is {repository_status}."

# Final Execution
if __name__ == "__main__":
    print("Executing Randomized Motif Search...")
    # Discovered motifs will vary slightly due to stochastic nature
    print(f"Discovered Motifs: {RandomizedMotifSearch(dna_input, 8, 5)}")
    print(ProjectStatus())
