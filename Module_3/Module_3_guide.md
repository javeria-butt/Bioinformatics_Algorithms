# 🧬 Module 3 Guide: Structural Reinforcement
### *Pouring the Concrete Slurry into the Bioinformatics Skeleton*

In this module, we integrate our individual functions into "Umbrella" programs. We introduce **NumPy** for high-performance matrix manipulation and build toward the **Greedy Motif Search**—an algorithm that can actually "discover" hidden regulatory signals in unknown DNA sequences.

---

### 🏗️ Phase 3.0: NumPy Foundations
**Concept:** Using high-efficiency arrays to handle DNA Motif Matrices.

A Motif Matrix represents aligned DNA sequences where each row is a sequence and each column is a genomic position. NumPy allows us to "slice" these columns instantly without complex loops.

```python
import numpy as np

# 1. Define a motif matrix (aligned DNA sequences)
motif = [
    ['A', 'T', 'G', 'C'],
    ['A', 'G', 'T', 'G'],
    ['T', 'C', 'C', 'A'],
    ['G', 'T', 'A', 'G'],
    ['C', 'G', 'A', 'T']
]

# 2. Convert to NumPy Array for high-speed slicing
motif_array = np.array(motif)

# 3. Extracting the 3rd column (Index 2)
# The syntax [:, 2] means "All rows, but only the 2nd index column"
extracted_col = motif_array[:, 2]

print(f"Extracted Column: {extracted_col}")
print(f"Full Matrix:\n{motif_array}")
```
### 🏗️ Phase 3.1: Nucleotide Frequency Counting
**Concept:** Calculating the "conservation" of A, C, G, and T at every position.

Before we can find a motif, we must count how often each nucleotide appears in a specific alignment. This **Count** function is the fundamental step in generating **Position Weight Matrices (PWMs)** and helps us visualize which parts of a DNA sequence are biologically preserved.

#### 🧠 The Logic
1.  **Initialize:** Create a dictionary with keys `A`, `C`, `G`, and `T`.
2.  **Map:** Each key points to a list of zeros equal to the length of the sequences.
3.  **Iterate:** Scan each sequence in the motif alignment.
4.  **Update:** For every nucleotide found at position `i`, increment the corresponding counter in the dictionary.

---

#### 🛠️ Implementation
This code demonstrates how to build a frequency map from a set of aligned DNA motifs.

```python
def Count(motif):
    """
    Counts the frequency of each nucleotide (A, C, G, T) at each position in the motif.
    
    Args:
        motif (list): A list of aligned DNA sequences (all same length).
        
    Returns:
        dict: A dictionary where keys are nucleotides and values are lists of frequencies.
    """
    # Step 1: Initialize the count dictionary with zeros for each nucleotide
    # The length of each list matches the length of the sequences (motif[0])
    k = len(motif[0])
    count = {nt: [0] * k for nt in "ACGT"}
    
    # Step 2: Iterate through each sequence in the motif set
    for sequence in motif:
        # Step 3: Iterate through each position (i) and nucleotide in the sequence
        for i, nucleotide in enumerate(sequence):
            # Increment the count for the corresponding nucleotide and position
            if nucleotide in count:
                count[nucleotide][i] += 1
    
    return count

# Example Motif Input (Aligned Sequences)
motif_set = [
    'AACGTA',
    'CCCGTT',
    'CACCTT',
    'GGATTA',
    'TTCCGG'
]

# Calculate and display nucleotide frequency counts
nucleotide_counts = Count(motif_set)
print("Nucleotide Frequency Count Mapping:")
for nt, counts in nucleotide_counts.items():
    print(f"{nt}: {counts}")

# Expected Output:
# A: [1, 2, 1, 0, 0, 2]
# C: [2, 1, 4, 2, 0, 0]
# G: [1, 1, 0, 2, 1, 1]
# T: [1, 1, 0, 1, 4, 2]
```
#### 💡 Significance
The resulting dictionary shows us the distribution of nucleotides. For instance, if position 2 has a high count for C (e.g., 4), it suggests that Cytosine is highly conserved at that spot, making it a likely part of a functional biological signal.

### 🏗️ Phase 3.2: The Profile Matrix
**Concept:** Converting nucleotide counts into probabilistic weights.

While raw counts tell us how often a nucleotide appears, a **Profile Matrix** converts those counts into probabilities. This transformation is the "Concrete Slurry" that allows us to score any new DNA sequence against our known motif.

#### 🧠 The Logic
1.  **Count:** Use the `Count` function (from Phase 3.1) to get the raw frequency of each nucleotide.
2.  **Calculate Total ($t$):** Identify the total number of sequences in the motif set.
3.  **Divide:** For every nucleotide at every position, divide the raw count by $t$.
4.  **Result:** A dictionary where every value is a probability between $0.0$ and $1.0$.

---

#### 🛠️ Implementation
This function extends our previous counting logic to generate a probability-based profile.

```python
def Count(motif):
    # Standard count logic from Phase 3.1
    count = {nt: [0] * len(motif[0]) for nt in "ACGT"}
    for sequence in motif:
        for i, nucleotide in enumerate(sequence):
            count[nucleotide][i] += 1
    return count

def Profile(motif):
    """
    Calculates the profile matrix for a given motif by converting 
    nucleotide counts to probabilities.
    
    Args:
        motif (list): A list of DNA sequences of the same length.
        
    Returns:
        dict: A dictionary where keys are nucleotides and values 
              are lists of probabilities (0.0 to 1.0).
    """
    t = len(motif)      # Total number of sequences
    k = len(motif[0])   # Length of each sequence
    profile = {}        # Initialize the profile matrix dictionary
    
    # Step 1: Get raw nucleotide counts
    count_matrix = Count(motif)
    
    # Step 2: Convert counts to probabilities for each nucleotide
    for nt in "ACGT":
        # Divide each count in the list by the total number of sequences (t)
        profile[nt] = [count / t for count in count_matrix[nt]]
    
    return profile

# Example Motif Input
motif_set = [
    'AACGTA',
    'CCCGTT',
    'CACCTT',
    'GGATTA',
    'TTCCGG'
]

# Calculate and display the profile matrix
profile_matrix = Profile(motif_set)
print("Profile Matrix (Probabilities):")
for nt, probabilities in profile_matrix.items():
    print(f"{nt}: {probabilities}")
```
#### 💡 Significance
The `Profile Matrix` is the DNA "Signature." If a position has a probability of 0.8 for C, it means that any sequence we find in the future that has a C at that spot is very likely to be part of this specific motif family.

### 🏗️ Phase 3.3: Consensus Sequence Discovery
**Concept:** Identifying the "most likely" representative sequence for a motif.

The **Consensus Sequence** is the ultimate summary of a DNA alignment. It is constructed by picking the most frequent nucleotide at each position. This provides a single, "ideal" sequence that represents the biological signal found across multiple DNA strands.

#### 🧠 The Logic
1.  **Count:** Use the `Count` function to calculate the frequency of A, C, G, and T at every position.
2.  **Determine Max:** For each column in the alignment, identify which nucleotide has the highest count.
3.  **Construct:** Append that "winning" nucleotide to a growing string.
4.  **Result:** A single DNA sequence that reflects the most common bases.

---

#### 🛠️ Implementation
This function extracts the dominant biological signal from a set of aligned motifs.

```python
def Count(motif):
    """
    Counts the frequency of each nucleotide (A, C, G, T) 
    at each position in the motif.
    """
    count = {nt: [0] * len(motif[0]) for nt in "ACGT"}
    for sequence in motif:
        for i, nucleotide in enumerate(sequence):
            count[nucleotide][i] += 1
    return count

def Consensus(motif):
    """
    Finds the consensus sequence for a given motif.
    
    Args:
        motif (list): A list of DNA sequences of the same length.
        
    Returns:
        str: The consensus sequence representing the most likely bases.
    """
    k = len(motif[0])    # Length of sequences
    count = Count(motif) # Get raw nucleotide counts
    consensus = ""       # Initialize the consensus sequence string
    
    # Iterate through each position (j) in the alignment
    for j in range(k):
        max_count = -1
        frequent_symbol = ""
        
        # Compare A, C, G, and T to find the winner for this position
        for symbol in "ACGT":
            if count[symbol][j] > max_count:
                max_count = count[symbol][j]
                frequent_symbol = symbol
        
        # Append the most frequent nucleotide to our consensus
        consensus += frequent_symbol
        
    return consensus

# Example Motif Input
motif_set = [
   'AACGTA',
   'CCCGTT',
   'CACCTT',
   'GGATTA',
   'TTCCGG'
]

# Calculate and display the consensus sequence
consensus_result = Consensus(motif_set)
print(f"Consensus Sequence: {consensus_result}")
# Expected Output: CACCTA
```
#### 💡 Significance
The Consensus Sequence `(e.g., CACCTA)` acts as the "Standard" for our motif. It allows us to compare other sequences against a single benchmark. In later phases, we use this sequence to calculate how much "variation" or "mutation" exists within our biological building.

### 🏗️ Phase 3.4: Scoring DNA Motifs
**Concept:** Measuring how much each sequence in a motif deviates from the consensus.

Scoring is the "Stress Test" for our biological structure. It identifies the total number of mismatches between the individual DNA strands in an alignment and their ideal **Consensus Sequence**. 

---

#### 🛠️ Implementation
This function calculates the level of conservation within a motif by summing every mismatch across all sequences.

```python
def Count(motif):
    """Counts frequencies of nucleotides at each position."""
    count = {nt: [0] * len(motif[0]) for nt in "ACGT"}
    for sequence in motif:
        for i, nt in enumerate(sequence):
            count[nt][i] += 1
    return count

def Consensus(motif):
    """Finds the most frequent nucleotide for each position."""
    k = len(motif[0])
    count = Count(motif)
    consensus = ""
    for j in range(k):
        frequent_symbol = max("ACGT", key=lambda nt: count[nt][j])
        consensus += frequent_symbol
    return consensus

def Score(motif):
    """
    Calculates the total mismatch score of a motif alignment.
    
    Returns:
        int: Total number of mismatches relative to the consensus.
        
    BIOLOGICAL SIGNIFICANCE:
    ------------------------
    - Low Score: High Conservation. This indicates a strong biological signal 
      (e.g., a regulatory site that evolution has kept consistent).
    - High Score: Low Conservation. This suggests the sequences are 
      diverse or potentially random noise.
    """
    # Step 1: Find the ideal consensus sequence
    consensus = Consensus(motif)
    k = len(motif[0])
    total_score = 0
    
    # Step 2: Compare each sequence at every position to the consensus
    for j in range(k):
        for sequence in motif:
            if sequence[j] != consensus[j]:
                total_score += 1  # Add 1 for every mismatch
    
    return total_score

# Example Execution
motif_data = [
    'AACGTA',
    'CCCGTT',
    'CACCTT',
    'GGATTA',
    'TTCCGG'
]

print(f"Consensus Sequence: {Consensus(motif_data)}")
print(f"Motif Score: {Score(motif_data)}")
# Output: Consensus: CACCTA, Score: 14
```
### 🏗️ Phase 3.5: Profile-Most Probable K-mer
**Concept:** Searching for the "Best Fit" genomic window using probabilistic weights.

The **Profile-Most Probable K-mer** function is the bridge between a known biological signal (the Profile) and an unknown DNA sequence (the Text). It identifies which $k$-length substring is most likely to have been generated by our specific profile matrix.

---

#### 🛠️ Implementation
This function scans a sequence and uses probability multiplication to find the $k$-mer that most closely matches the expected motif signature.

```python
def Pr(Text, Profile):
    """
    Calculates the probability of a specific k-mer being 
    generated by a given profile matrix.
    
    Args:
        Text (str): A k-mer substring.
        Profile (dict): A probability matrix for each position.
        
    Returns:
        float: The cumulative probability (0.0 to 1.0).
    """
    p = 1  # Initialize probability at 1 (Neutral for multiplication)
    for i in range(len(Text)):
        # Multiply probabilities for each nucleotide at its specific position
        p *= Profile[Text[i]][i]
    return p

def ProfileMostProbableKmer(Text, k, Profile):
    """
    Finds the single most likely k-mer in Text based on a profile matrix.
    
    Returns:
        str: The k-mer with the highest probability.
        
    BIOLOGICAL SIGNIFICANCE:
    ------------------------
    - In motif discovery, we use this to 'hunt' for a regulatory site 
      within a long strand of DNA. 
    - It allows us to pick the sequence that most closely 'resembles' 
      the profile we've built from other known motifs.
    """
    max_probability = -1  # Initialize with a value lower than possible
    most_probable_kmer = Text[0:k]  # Default to the first k-mer

    # Step 1: Slide through the Text to consider all possible k-mers
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i + k]  
        
        # Step 2: Calculate the probability for this specific window
        probability = Pr(kmer, Profile)

        # Step 3: Update the winner if a higher probability is found
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer

    return most_probable_kmer

# Example Execution
dna_text = "ACGTACGTGACG"
k_length = 3
profile_map = {
    'A': [0.2, 0.4, 0.3],
    'C': [0.4, 0.3, 0.1],
    'G': [0.3, 0.1, 0.5],
    'T': [0.1, 0.2, 0.1]
}

result_kmer = ProfileMostProbableKmer(dna_text, k_length, profile_map)
print(f"Most Probable {k_length}-mer: {result_kmer}")
# Expected Output: ACG
```
### 🏗️ Phase 3.6: Greedy Motif Search
**Concept:** The complete "Umbrella" algorithm for biological motif discovery.

This is the culmination of **Module 3**. The **Greedy Motif Search** is an iterative algorithm that identifies hidden regulatory signals across multiple DNA sequences. It "greedily" builds a profile from the best k-mers it finds, refining its search sequence-by-sequence to minimize the total score.

---

#### 🛠️ Implementation
This algorithm acts as a master controller, integrating `Profile`, `Score`, and `ProfileMostProbableKmer` into a single, high-level search engine.

```python
def GreedyMotifSearch(DNA, k, t):
    """
    Performs an iterative search for the best motif set across multiple sequences.
    
    Args:
        DNA (list): Multiple DNA strands to search.
        k (int): Desired motif length.
        t (int): Number of sequences in the DNA list.
        
    Returns:
        list: The set of k-mers that minimize the consensus score.
        
    BIOLOGICAL SIGNIFICANCE:
    ------------------------
    - This algorithm is a 'conglomeration' of our modular functions. 
    - It allows us to discover conserved sequences (like TATA boxes or 
      transcription factor binding sites) even when we don't know the 
      consensus sequence beforehand.
    - It works by 'guessing' a k-mer in the first strand and seeing how well 
      subsequent strands align with it.
    """
    n = len(DNA[0])
    # Step 1: Initialize best_motifs with the first k-mers from each sequence
    best_motifs = [DNA[i][0:k] for i in range(t)]

    # Step 2: Iterate through every possible k-mer in the FIRST sequence
    for i in range(n - k + 1):
        # Start a candidate motif set with one k-mer from the first sequence
        motifs = [DNA[0][i:i + k]]
        
        # Step 3: Build the motif set for the remaining (t-1) sequences
        for j in range(1, t):
            # Create a profile based ONLY on the motifs selected so far
            current_profile = Profile(motifs[:j])
            
            # Find the most probable k-mer in the NEXT sequence using that profile
            next_motif = ProfileMostProbableKmer(DNA[j], k, current_profile)
            motifs.append(next_motif)
            
        # Step 4: Check if this new set of motifs is better than our previous best
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs

    return best_motifs

# --- Helper Functions Required for Integration ---
def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return p

def ProfileMostProbableKmer(Text, k, Profile):
    max_p = -1
    best_k = Text[0:k]
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i+k]
        prob = Pr(kmer, Profile)
        if prob > max_p:
            max_p, best_k = prob, kmer
    return best_k

def Profile(motifs):
    t, k = len(motifs), len(motifs[0])
    profile = {nt: [0]*k for nt in "ACGT"}
    for m in motifs:
        for i, nt in enumerate(m):
            profile[nt][i] += 1
    for nt in "ACGT":
        profile[nt] = [x / t for x in profile[nt]]
    return profile

def Score(motifs):
    k = len(motifs[0])
    count = {nt: [0]*k for nt in "ACGT"}
    for m in motifs:
        for i, nt in enumerate(m):
            count[nt][i] += 1
    consensus = ""
    for j in range(k):
        consensus += max("ACGT", key=lambda nt: count[nt][j])
    
    total_score = 0
    for m in motifs:
        for i in range(k):
            if m[i] != consensus[i]:
                total_score += 1
    return total_score

# Example Execution
dna_list = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
result = GreedyMotifSearch(dna_list, 3, 5)
print(f"Best Motifs Discovered: {result}")
# Expected Output: ['CAG', 'CAG', 'CAA', 'CAA', 'CAA']
```
