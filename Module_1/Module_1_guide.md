# 🧬 Module 1: Substratum & Excavation
### *Digging into the Foundations of Sequence Analysis*

Welcome to the first layer of the Bio-Architect project. Before we can build complex search algorithms, we must understand how to manually handle the raw "dirt" of bioinformatics: DNA strings. 

In this module, we bypass the "cheat codes" of BioPython to write the manual logic for pattern recognition and genetic mapping.

---

### 🏗️ Phase 1.1: The Pattern Counter
**Concept:** Identifying the core frequency of a specific motif within a parent sequence.

To find how many times a specific pattern (e.g., "ATA") exists in a genome, we slide a window across the text and compare segments.

```python
# Dataset
genome_text = 'ATGCATATGACTACTAGATACTGATACTGATACATA'
target_pattern = 'ATA'

def count_pattern(text, pattern):
    count = 0
    # Sliding window approach
    for i in range(len(text) - len(pattern) + 1):
        if text[i : i + len(pattern)] == pattern:
            count += 1 
    return count
```
### 🏗️ Phase 1.2: Genomic Frequency Mapping
**Concept:** Constructing a dictionary-based map of every k-mer (sequence of length k).

Frequency maps serve as the "blueprints" for discovering hidden motifs. By mapping every possible k-mer to its frequency count, we can visualize the distribution of genetic data across a sequence.

```python
def get_frequency_map(text, k):
    frequencies = {}
    n = len(text)
    for i in range(n - k + 1):
        kmer = text[i : i + k]
        # Using .get() allows us to handle new k-mers and 
        # increment existing ones in a single, clean line.
        frequencies[kmer] = frequencies.get(kmer, 0) + 1
    return frequencies

# Example Usage:
sample_dna = 'ACGTTGCATGTCGCATGAGCATGAGAGCT'
k_length = 3
print(get_frequency_map(sample_dna, k_length))
```
### 🏗️ Phase 1.3: Identifying Most Frequent K-mers
**Concept:** Extracting the "Peaks" from our frequency map.

By combining the mapping logic from Phase 1.2 with a filter for maximum values, we can pinpoint the most significant motifs (the "Highest Points") in any genomic sequence. This allows us to identify potential biological signals or binding sites.

```python
def find_most_frequent_kmers(text, k):
    frequent_words = []
    # Step 1: Generate the blueprint (Frequency Map)
    freq_map = get_frequency_map(text, k)
    
    # Step 2: Identify the maximum frequency count
    max_count = max(freq_map.values())
    
    # Step 3: Extract all k-mers that hit that peak
    for kmer, count in freq_map.items():
        if count == max_count:
            frequent_words.append(kmer)
    
    return frequent_words

# Example Usage:
dna_input = 'ACGTTGCATGTCGCATGAGCATGAGAGCT'
k_size = 3
print(f"Top Motifs Found: {find_most_frequent_kmers(dna_input, k_size)}")
# Expected Output: ['GCA', 'CAT', 'ATG']
```
### 🏗️ Phase 1.4: Manual DNA Complementation
**Concept:** Replicating the logic of `Bio.Seq.complement()`.

While BioPython allows this in a single line, writing the logic manually helps us understand the Python-to-Biology interface—specifically the A ↔ T and C ↔ G relationship. By building the string one base at a time, we mimic the base-pairing rules of DNA replication.

```python
def generate_complement(pattern):
    # Mapping dictionary for high-speed base-pair lookup
    pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # We build the complement string by iterating through each base
    complement_sequence = ''
    for base in pattern:
        complement_sequence += pairs.get(base, 'N') # 'N' for unknown bases
        
    return complement_sequence

# Example Usage:
original_dna = 'ACGTTGCATGTCGCATGAGCATGAGAGCT'
print(f"Sequence:   {original_dna}")
print(f"Complement: {generate_complement(original_dna)}")
# Expected Output: TGCAACGTACAGCGTACTCGTACTCTCGA
```
### 🏗️ Phase 1.5: Pattern Localization (Matching)
**Concept:** Identifying the exact starting coordinates of a pattern within a genome.

Knowing the *frequency* of a motif is vital, but knowing the *location* is crucial for identifying gene regions, promoters, or binding sites within a large genome. By searching the sequence and recording the starting index of every match, we can pinpoint precisely where these biological signals are located.

```python
def locate_pattern(pattern, genome):
    start_positions = []
    
    # We iterate through the genome, sliding a window of the same length as our pattern
    for i in range(len(genome) - len(pattern) + 1):
        # We compare the current window segment to our target pattern
        if genome[i : i + len(pattern)] == pattern:
            start_positions.append(i) # Store the index of the match
            
    return start_positions

# Example Usage:
search_text = 'CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC'
search_motif = 'CGCG'

print(f"Coordinates found: {locate_pattern(search_motif, search_text)}")
# Expected Output: [0, 16, 28, 38]
```
