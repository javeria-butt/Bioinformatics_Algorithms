# 🧬 Module 2 Guide: Structural Framework
### *Adding Steel Rods to the Bioinformatics Skeleton*

In this module, we transition from basic sequence manipulation to complex algorithmic structures. We are no longer just looking at strings; we are analyzing **genomic distribution**, **replication origins**, and **evolutionary distance**. These are the "steel rods" that allow our bioinformatics building to stand tall against large-scale genomic data.

---

### 🏗️ Phase 2.1: The Symbol Array
**Concept:** Analyzing nucleotide enrichment in circular genomes.

The Symbol Array identifies how frequently a specific nucleotide (e.g., "A" or "C") appears within a sliding window. By extending the genome string, we simulate a **circular genome**, ensuring that we don't miss patterns that "wrap around" the end of the sequence.

```python
def count_pattern(pattern, text):
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i+len(pattern)] == pattern:
            count += 1
    return count

def get_symbol_array(genome, symbol):
    array = {}
    n = len(genome)
    # Circularized extension: appending the first half to the end
    extended_genome = genome + genome[0:n//2]
    
    for i in range(n):
        # Analyzing a window of half the genome length
        window = extended_genome[i : i + (n//2)]
        array[i] = count_pattern(symbol, window)
    return array

# Example Execution
genome_sample = "ATGATAGTCCGAAA"
print(f"Symbol Array for 'A': {get_symbol_array(genome_sample, 'A')}")
```
### 🏗️ Phase 2.2: The Skew Array
**Concept:** Calculating the cumulative G-C imbalance to track DNA replication.

The **Skew** function measures the balance between Guanine (G) and Cytosine (C) across a genome. In Bioinformatics, this is a critical tool for identifying the "direction" of DNA strands.
* **G** increases the skew value by **+1**.
* **C** decreases the skew value by **-1**.
* **A** and **T** are neutral (0).

#### 🧠 The Logic of Skewness
* **Skew = 0**: The distribution is balanced.
* **Skew > 0**: There is a higher concentration of Guanine (G).
* **Skew < 0**: There is a higher concentration of Cytosine (C).

---

#### 🛠️ Implementation - Method 1: List-Based Approach
This method is ideal for sequential data processing and creating visualizations.

```python
def SkewArray(Genome):
    # Initialize skew array starting at 0 balance
    skew = [0]
    
    for i in range(len(Genome)):
        if Genome[i] == 'G':
            skew.append(skew[-1] + 1)  # Increment for G
        elif Genome[i] == 'C':
            skew.append(skew[-1] - 1)  # Decrement for C
        else:
            skew.append(skew[-1])      # No change for A or T
            
    return skew

# Example Execution
print(SkewArray("CATGGGCATCGGCCATACGCC"))
```
#### 🛠️ Implementation - Method 2: Dictionary-Based Approach
This method allows for quick coordinate-based lookups (Key = Position, Value = Skew).
```
def SkewArrayDict(Genome):
    # Initialize dictionary with position 0 mapped to skew 0
    skew = {0: 0}
    
    for i in range(len(Genome)):
        if Genome[i] == 'G':
            skew[i + 1] = skew[i] + 1
        elif Genome[i] == 'C':
            skew[i + 1] = skew[i] - 1
        else:
            skew[i + 1] = skew[i]
            
    return skew

# Example Execution
print(SkewArrayDict("GATACACTTCCCGAGTAGGTACTG"))
```
#### 💡 Key Takeaway
The Skew Array provides a cumulative trail of nucleotide distribution. The point where this array reaches its minimum value is a major biological landmark, often indicating the potential Origin of Replication (oriC), which we will explore in the next Phase.

### 🏗️ Phase 2.3: Finding the Minimum Skew
**Concept:** Pinpointing the potential "Origin of Replication" (*oriC*).

The positions of **Minimum Skew** represent the regions where the concentration of Cytosine is at its highest relative to Guanine. In many bacterial genomes, this specific point marks the biological landmark where DNA replication begins.

#### 🧠 The Logic
1.  **Calculate:** Generate the complete Skew Array (using the logic from Phase 2.2).
2.  **Identify:** Find the global minimum value within that array.
3.  **Collect:** Record every genomic position that hits that minimum value.

---

#### 🛠️ Implementation
This function integrates the `SkewArray` logic to extract the specific coordinates of the minimum skew.

```python
def SkewArray(Genome):
    # Helper: Calculates the cumulative skew
    skew = {0: 0}
    for i in range(len(Genome)):
        if Genome[i] == 'G':
            skew[i + 1] = skew[i] + 1
        elif Genome[i] == 'C':
            skew[i + 1] = skew[i] - 1
        else:
            skew[i + 1] = skew[i]
    return skew

def MinimumSkew(Genome):
    """
    Identifies the exact positions in the genome where the skew is at its minimum.
    
    Returns:
        list: A collection of indices where the replication origin is likely located.
    """
    positions = []
    # Step 1: Generate the full skew map
    skew_data = SkewArray(Genome)
    
    # Step 2: Find the lowest value in the dictionary
    min_val = min(skew_data.values())
    
    # Step 3: Collect all keys (positions) that share this minimum value
    for position, value in skew_data.items():
        if value == min_val:
            positions.append(position)
            
    return positions

# Example Execution
sample_genome = "GATACACTTCCCGAGTAGGTACTG"
print(f"Minimum Skew Position(s): {MinimumSkew(sample_genome)}")
# Expected Output: [12]
```

#### 💡 Biological Significance
> The output `[12]` signifies that at position 12, the cumulative **G-C imbalance** is at its lowest. These coordinates are the primary candidates for researchers looking to identify the **Origin of Replication** (*oriC*) in a circular bacterial chromosome.

```python
# The result [12] indicates that at genomic coordinate 12, 
# the skew reaches its global minimum. 

# Biologically, this point is often the location where 
# DNA replication initiates (the Origin of Replication).
```
### 🏗️ Phase 2.4: Hamming Distance
**Concept:** Quantifying evolutionary "distance" through single nucleotide mutations.

The **Hamming Distance** between two DNA sequences of equal length is the total count of positions where the corresponding nucleotides differ. In Bioinformatics, this metric is the fundamental way we measure mutations and genetic divergence between species.

#### 🧠 The Logic
1.  **Iterate:** Traverse both sequences simultaneously from the first position to the last.
2.  **Compare:** Check if the base at index `i` in the first sequence matches the base at index `i` in the second.
3.  **Count:** Increment a counter for every mismatch found.

---

#### 🛠️ Implementation
This function provides a direct calculation of the genomic distance between two sequences.

```python
def HammingDistance(seq1, seq2):
    """
    Calculates the Hamming Distance (number of mismatches) between two equal-length DNA strands.
    
    Args:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence
        
    Returns:
        int: Total count of differing positions (Mutations)
    """
    # Step 1: Initialize a counter for mismatches
    mismatch_counter = 0

    # Step 2: Compare bases at each position
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            mismatch_counter += 1  # Increment for each mutation found

    return mismatch_counter

# Example Execution
sequence_a = "CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG"
sequence_b = "ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT"

distance = HammingDistance(sequence_a, sequence_b)
print(f"Total Mutations (Hamming Distance): {distance}")
# Expected Output: 43
```
### 🏗️ Phase 2.5: Approximate Pattern Matching
**Concept:** Identifying motifs while allowing for evolutionary mutations.

In nature, DNA sequences rarely match a target pattern perfectly due to mutations. **Approximate Pattern Matching** allows us to find locations where a pattern almost matches a specific region, provided the number of differences is within a defined limit (the Hamming distance threshold $d$).

#### 🧠 The Logic
1.  **Define Inputs:** The target **Text**, the **Pattern** to search for, and the maximum allowed mismatches (**d**).
2.  **Iterate:** Extract every possible substring of the same length as the pattern using a sliding window.
3.  **Compare:** Calculate the **Hamming Distance** between the current window and the target pattern.
4.  **Record:** If the Hamming Distance is $\le d$, store the starting coordinate of that window.

---

#### 🛠️ Implementation
This function integrates the `HammingDistance` logic to perform a search that accounts for genetic variations.

```python
def HammingDistance(seq1, seq2):
    # Counts occurrences of mismatches between two equal-length sequences
    counter = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            counter += 1
    return counter

def ApproximatePatternMatching(Text, Pattern, d):
    """
    Finds all positions where a pattern matches the text with up to 'd' mismatches.
    
    Args:
        Text (str): The DNA sequence to search through.
        Pattern (str): The motif or sequence to search for.
        d (int): Maximum allowed mismatches (Hamming distance).
        
    Returns:
        list: Positions where the pattern approximately matches the text.
    """
    positions = []  # Store matching positions
    n = len(Text)
    m = len(Pattern)

    # Check all substrings of Text with length equal to the pattern
    for i in range(n - m + 1):
        # Extract the substring and calculate Hamming distance
        if HammingDistance(Text[i:i + m], Pattern) <= d:
            positions.append(i)  # Add position if within mismatch threshold

    return positions

# Example Execution
Text = "ATGCATATGACTACTAGATACTGATACTGATACATA"
Pattern = "ATAG"
d = 1
print(f"Approximate Matches: {ApproximatePatternMatching(Text, Pattern, d)}")
# Expected Output: [4, 13, 17, 23, 29]
```
#### 💡 Biological Significance
The output [4, 13, 17, 23, 29] tells us that the pattern "ATAG" appears approximately at these positions in the sequence, where each match has exactly 1 mismatch ($d=1$). This is vital for Motif Discovery, as it allows us to find functional signals even if they have slightly evolved over time.
