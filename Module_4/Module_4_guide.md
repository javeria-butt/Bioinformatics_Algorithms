"""
🧬 MODULE 4: PROBABILISTIC OPTIMIZATION AND STOCHASTIC SEARCH
------------------------------------------------------------
This module introduces robust statistical methods to handle biological 
variability and zero-probability errors in motif discovery.

PHASE 4.1: GREEDY MOTIF SEARCH WITH PSEUDOCOUNTS
------------------------------------------------
Concept: Implementing Laplace's Rule of Succession to stabilize the Profile Matrix.
In standard counting, a missing nucleotide (count=0) results in a probability 
of 0.0, which 'extinguishes' the score of an entire DNA strand. Pseudocounts 
add a small buffer (+1) to every base to ensure no signal is ignored.
"""

def ProfileWithPseudocounts(motifs):
    """
    Generates a profile matrix using a Laplace-style buffer.
    
    BIOLOGICAL SIGNIFICANCE:
    ------------------------
    - STABILIZATION: Initializing nucleotide counts at 1 instead of 0 ensures 
      that no genomic signal is completely ignored due to absence in small 
      sample sets.
    - NORMALIZATION: The total count is adjusted to (t + 4) to account 
      for the four virtual nucleotides (A, C, G, T) added to the matrix.
    """
    t = len(motifs)
    k = len(motifs[0])
    # Step 1: Initialize counts with 1 (the 'Pseudocount')
    profile = {nt: [1] * k for nt in "ACGT"}  
    
    # Step 2: Aggregate counts from the motif alignment
    for motif in motifs:
        for i, nt in enumerate(motif):
            profile[nt][i] += 1
            
    # Step 3: Normalize counts into probabilities
    for nt in "ACGT":
        profile[nt] = [x / (t + 4) for x in profile[nt]]
    return profile

def GreedyMotifSearchWithPseudocounts(DNA, k, t):
    """
    An iterative search that utilizes pseudocounts to find robust motifs.
    
    Args:
        DNA (list): Set of DNA sequences.
        k (int): Length of the target motif.
        t (int): Total number of sequences.
    """
    n = len(DNA[0])
    # Initialize with the first available k-mers
    best_motifs = [DNA[i][:k] for i in range(t)]

    for i in range(n - k + 1):
        # Begin with a candidate k-mer from the first sequence
        motifs = [DNA[0][i:i + k]]
        for j in range(1, t):
            # Generate the 'buffered' profile from current motif set
            profile = ProfileWithPseudocounts(motifs)
            # Find the most probable k-mer in the next sequence
            motifs.append(ProfileMostProbableKmer(DNA[j], k, profile))
            
        # Update if the new motif set yields a lower mismatch score
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs

    return best_motifs

"""
PHASE 4.2: RANDOMIZED MOTIF SEARCH
----------------------------------
Concept: Utilizing stochastic exploration to bypass local optima.
Deterministic algorithms can get 'trapped' in local minima. Randomized 
Search starts with random coordinates and iteratively refines them, 
allowing for a broader exploration of the genomic search space.
"""

import random

def RandomizedMotifSearch(DNA, k, t):
    """
    Performs stochastic motif discovery through iterative refinement.
    
    BIOLOGICAL SIGNIFICANCE:
    ------------------------
    - STOCHASTIC INITIALIZATION: Starting with random coordinates allows 
      the algorithm to explore motifs that might be missed by a linear scan.
    - ITERATIVE CONVERGENCE: The algorithm 'descends' toward the lowest 
      score (highest conservation) and stops at the local peak performance.
    """
    # Step 1: Initialize with random k-mers from the DNA strands
    motifs = RandomMotifs(DNA, k)
    best_motifs = motifs

    while True:
        # Step 2: Build a profile from the current selection
        profile = ProfileWithPseudocounts(motifs)
        
        # Step 3: Update the motif set based on the current profile probabilities
        motifs = [ProfileMostProbableKmer(seq, k, profile) for seq in DNA]
        
        # Step 4: Retain the set only if the mismatch score improves
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs
        else:
            # Return the best motifs found when no further improvement is possible
            return best_motifs

def RandomMotifs(DNA, k):
    """Selects one random k-mer from each sequence in the input list."""
    motifs = []
    for sequence in DNA:
        start = random.randint(0, len(sequence) - k)
        motifs.append(sequence[start:start + k])
    return motifs

# --- CORE UTILITY FUNCTIONS (Required for Module 4 Integration) ---

def Pr(Text, profile):
    """Calculates cumulative probability of a k-mer given a profile."""
    p = 1
    for i in range(len(Text)):
        p *= profile[Text[i]][i]
    return p

def ProfileMostProbableKmer(Text, k, profile):
    """Finds the k-mer with the highest probability in a sequence."""
    max_p, best_k = -1, Text[0:k]
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i+k]
        prob = Pr(kmer, profile)
        if prob > max_p:
            max_p, best_k = prob, kmer
    return best_k

def Score(motifs):
    """Calculates total mismatches between motifs and their consensus."""
    k = len(motifs[0])
    count = {nt: [0]*k for nt in "ACGT"}
    for m in motifs:
        for i, nt in enumerate(m): count[nt][i] += 1
    # Generate Consensus Sequence
    consensus = "".join([max("ACGT", key=lambda nt: count[nt][j]) for j in range(k)])
    # Sum mismatches
    return sum(1 for m in motifs for i in range(k) if m[i] != consensus[i])

# --- EXAMPLE EXECUTION ---

if __name__ == "__main__":
    dna_input = [
        "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
        "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
        "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
        "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
        "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
    ]
    k_val, t_val = 8, len(dna_input)
    
    print("Executing Randomized Motif Search...")
    results = RandomizedMotifSearch(dna_input, k_val, t_val)
    print(f"Optimal Motifs Found: {results}")
    print(f"Final Score: {Score(results)}")
    
