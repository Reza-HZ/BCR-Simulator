```mermaid
graph TD
    A[Start] --> B[Initialize BCRSimulator]
    B --> C[Generate root BCR (naive)]
    C --> D{Set root affinity and abundance}
    D --> E[Initialize all_bcrs with root]
    E --> F[Set current_generation = [root]]
    F --> G[generation = 1]
    G --> H{generation <= max_generations?}
    H -->|Yes| I[Initialize next_generation = []]
    I --> J{For each parent in current_generation}
    J --> K[Set parent_abd = round(parent['abundance'])]
    K --> L{i = 0; i < parent_abd}
    L --> M[Mutate parent sequence<br>(mutated_seq, num_mutations)]
    M --> N[Translate mutated_seq to a_sequence<br>Check for stop codons (min_stop)]
    N --> O[Calculate affinity based on align_type]
    O --> P[Set temp_aff_thresh =<br>max((affinity_threshold + 3*parent_aff) / 4, affinity_threshold)]
    P --> Q{Determine selection probability (selection_p)<br>based on affinity and min_stop}
    Q --> R{random.random() > selection_p?}
    R -->|Yes| S[i += 1]
    S --> L
    R -->|No| T[Calculate expected_abundance =<br>round(max_abundancy ** (affinity ** 3))]
    T --> U[Set abundance = max(np.random.poisson(expected_abundance), 1)]
    U --> V{num_mutations == 0?}
    V -->|Yes| W[parent_abd += abundance]
    V -->|No| X[Create child BCR with:<br>id, sequence, a_sequence, generation,<br>parent, affinity, abundance, mutations]
    X --> Y[Append child to all_bcrs and next_generation]
    Y --> Z[i += 1]
    Z --> L
    L -->|i >= parent_abd| AA[Update parent['abundance'] = parent_abd]
    AA --> J
    J -->|No more parents| BB{next_generation empty?}
    BB -->|Yes| CC[Merge sequences (merge_sequences)]
    BB -->|No| DD[current_generation = next_generation]
    DD --> EE[generation += 1]
    EE --> H
    H -->|No| CC
    CC --> FF[Return all_bcrs, merged_all_bcrs]
    FF --> GG[End]
