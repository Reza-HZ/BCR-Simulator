import os
import csv
import re
import numpy as np
import random
import matplotlib.pyplot as plt
from collections import defaultdict
from Bio import Phylo
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
from Bio import BiopythonWarning
import warnings
from ete3 import Tree
from matplotlib import cm
import mplcursors
from itertools import product
from concurrent.futures import ThreadPoolExecutor
warnings.filterwarnings("ignore", message="Attempting to set identical low and high xlims*")
warnings.simplefilter('ignore', BiopythonWarning)


class Pathogen:
    def __init__(self, name, n_epitopes, a_epitopes):
        self.name = name
        self.n_epitopes = n_epitopes
        self.a_epitopes = a_epitopes
        self.n_aligner = self._init_n_aligner()
        self.a_aligner = self._init_a_aligner()

    def _init_n_aligner(self):
        aligner = PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -1
        return aligner

    def _init_a_aligner(self):
        aligner = PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -1
        aligner.mode = 'local'
        return aligner

    def calculate_affinity_n(self, bcr_nsequence):
        max_affinity = 0
        for epitope in self.n_epitopes:
            alignment_score = self._local_n_alignment_score(bcr_nsequence, epitope)
            max_possible_score = len(epitope) * 2 
            normalized_score = alignment_score / max_possible_score if max_possible_score > 0 else 0
            max_affinity = max(max_affinity, normalized_score)
        return max_affinity
    
    def calculate_affinity_a(self, bcr_asequence):
        max_affinity = 0
        for epitope in self.a_epitopes:
            alignment_score = self._local_a_alignment_score_normalized(bcr_asequence, epitope)
            max_affinity = max(max_affinity, alignment_score)
        return max_affinity
    
    # def calculate_affinity(self, bcr_sequence):
    #     best_energy = float('-inf')
    #     for epitope in self.epitopes:
    #         alignment = self.aligner.align(bcr_sequence, epitope)[0]
    #         # Simple energy model: score - gap penalties
    #         energy = alignment.score - (alignment.gap_openings * 10 + alignment.gaps * 1)
    #         best_energy = max(best_energy, energy)
    #     return best_energy

    def _local_n_alignment_score(self, n_seq1, n_seq2):
        return self.n_aligner.score(n_seq1, n_seq2)

    def _local_a_alignment_score_normalized(self, a_seq1, a_seq2):
        raw_score = self.a_aligner.score(a_seq1, a_seq2)
        if raw_score == 0:
            return 0.0
        # shorter_seq = a_seq1 if len(a_seq1) <= len(a_seq2) else a_seq2
        # max_self_score = self.a_aligner.score(shorter_seq, shorter_seq)
        max_self_score = self.a_aligner.score(a_seq2, a_seq2)
        if max_self_score == 0:
            max_self_score = 1.0
        normalized_score = raw_score / max_self_score
        return max(0.0, min(1.0, normalized_score))


class BCRSimulator:
    def __init__(self):
        self.vdj_segments = self._initialize_vdj_segments()
        self.mutation_patterns = self._initialize_mutation_patterns()
        self.bcr_counter = 1

    def _initialize_vdj_segments(self):
        segments = {
            "V": read_fasta(r"K:\BCellPhylo\ClusteringBCtrees\Bcr_simulator\IGH genes\V.fasta"),
            "D": read_fasta(r"K:\BCellPhylo\ClusteringBCtrees\Bcr_simulator\IGH genes\D.fasta"),
            "J": read_fasta(r"K:\BCellPhylo\ClusteringBCtrees\Bcr_simulator\IGH genes\J.fasta")
        }
        return segments

    def _initialize_mutation_patterns(self):
        hotspots = {
            "WRC": 4,
            "GYW": 4,
            "RGYW": 5,
            "WRCY": 5
        }
        return hotspots

    def _generate_bcr_id(self):
        bcr_id = f"seq{self.bcr_counter}"
        self.bcr_counter += 1
        return bcr_id

    def _generate_initial_bcr(self, preproc_genes_index):
        if preproc_genes_index:
            weights = [item[-1] for item in preproc_genes_index]
            v_index, d_index, j_index, _ = random.choices(preproc_genes_index, weights=weights, k=1)[0]
            v_segment = self.vdj_segments["V"][v_index]
            d_segment = self.vdj_segments["D"][d_index]
            j_segment = self.vdj_segments["J"][j_index]
        else:
            v_segment = random.choice(self.vdj_segments["V"])
            d_segment = random.choice(self.vdj_segments["D"])
            j_segment = random.choice(self.vdj_segments["J"])

        while True:
            n1 = ''.join(random.choices("acgt", k=random.randint(0, 10)))
            n2 = ''.join(random.choices("acgt", k=random.randint(0, 10)))
            bcr_sequence = v_segment + n1 + d_segment + n2 + j_segment
            non_stop_regions = [(0, len(v_segment)), (len(v_segment)+len(n1), len(v_segment)+len(n1)+len(d_segment)), (len(bcr_sequence)-len(j_segment),len(bcr_sequence))]
            bcr_a_sequence, min_stop = translate_nucleotide_to_protein_min_stops(bcr_sequence, non_stop_regions)
            if min_stop == 0:
                break
        return {
            "sequence": bcr_sequence,
            "a_sequence": bcr_a_sequence,
            "v_gene": self.vdj_segments["V"].index(v_segment),
            "d_gene": self.vdj_segments["D"].index(d_segment),
            "j_gene": self.vdj_segments["J"].index(j_segment),
            "v_len": len(v_segment),
            "d_len": len(d_segment),
            "j_len": len(j_segment),
            "n1_len": len(n1),
            "n2_len": len(n2),
            "generation": 0,
            "mutations": 0,
            "parent": None,
            "affinity": 0.0,
            "id": "naive"
        }

    
    def _mutate_sequence(self, sequence, mutation_rate):
        nucleotides = list(sequence)
        mutations = 0

        HOTSPOT_TARGETS = {
        "WRC": 2,
        "GYW": 0,
        "RGYW": 1,
        "WRCY": 2
        }

        coldspots = {
            "SYC": 0.6,
            "GRS": 0.5
        }
        seq_len = len(nucleotides)
        for i in range(seq_len):
            raw_score = 0.0
            target_idx = i
            for motif, weight in self.mutation_patterns.items():
                sub_seq = sequence[i:i+len(motif)]
                if len(sub_seq) == len(motif) and self._match_hotspot(sub_seq, motif):
                    target_idx = i + HOTSPOT_TARGETS[motif]
                    raw_score += weight

            mutation_prob = raw_score * mutation_rate if raw_score > 0 else mutation_rate

            for cold, penalty in coldspots.items():
                sub_seq = sequence[i:i+len(cold)]
                if len(sub_seq) == len(cold) and self._match_hotspot(sub_seq, cold):
                    mutation_prob *= penalty
            
            for cold_motif, penalty in coldspots.items():
                motif_len = len(cold_motif)
                for offset in range(-motif_len + 1, 1):
                    start = target_idx + offset
                    end = start + motif_len
                    if start >= 0 and end <= len(sequence):
                        window = sequence[start:end]
                        if self._match_hotspot(window, cold_motif):
                            mutation_prob *= penalty
                            break

            if random.random() < mutation_prob and target_idx < seq_len:
                original = nucleotides[target_idx]
                if original in ['a','c','g','t']:
                    mutated = biased_mutation(original)
                    nucleotides[target_idx] = mutated
                    mutations += 1
        mutated_sequence = ''.join(nucleotides)
        if mutated_sequence == sequence:
            mutations = 0
        return mutated_sequence, mutations

    
    def _match_hotspot(self, sequence, hotspot_pattern):
        if len(sequence) != len(hotspot_pattern):
            return False

        for base, code in zip(sequence, hotspot_pattern):
            valid_bases = IUPAC_CODES.get(code, {code})
            if base not in valid_bases:
                return False
        return True


    def generate_repertoire_dynamically(self, pathogen, affinity_threshold=0.2, 
                                        max_generations=8, mutation_rate=0.001, preproc_genes_index = []):
        root = self._generate_initial_bcr(preproc_genes_index)
        root["affinity"] = pathogen.calculate_affinity_a(root["a_sequence"])
        root["abundance"] = 1
        all_bcrs = [root]
        current_generation = [root]
        generation = 1

        print(root["affinity"])

        while generation <= max_generations:
            next_generation = []
            for parent in current_generation:
                parent_aff = parent["affinity"]
                parent_abd = int(round(parent["abundance"]))
                i = 0
                while i < parent_abd:
                    for _ in range(2):
                        mutated_seq, num_mutations = self._mutate_sequence(parent["sequence"], mutation_rate)
                        v_l = root['v_len']
                        d_l = root['d_len']
                        j_l = root['j_len']
                        n1_l = root['n1_len']
                        if num_mutations > 0:
                            non_stop_regions = [(0, v_l), (v_l + n1_l, v_l + n1_l + d_l), (len(mutated_seq) - j_l, len(mutated_seq))]
                            mutated_a_seq, min_stop = translate_nucleotide_to_protein_min_stops(mutated_seq, non_stop_regions)
                            affinity = pathogen.calculate_affinity_a(mutated_a_seq)
                        else:
                            mutated_a_seq = parent['a_sequence']
                            min_stop = mutated_a_seq.count('*')
                            affinity = parent['affinity']
                        #temp_aff_thresh = max((affinity_threshold + 3*parent_aff) / 4, affinity_threshold)
                        min_aff, max_aff = sorted([parent_aff, affinity_threshold])
                        temp_aff_thresh = (((max_generations - generation) * min_aff) + generation * max_aff) / max_generations
                        abundance = 1
                        
                        if affinity >= temp_aff_thresh and min_stop == 0:
                            selection_p = 1
                        elif affinity < temp_aff_thresh and min_stop == 0:
                            selection_p = 0.4
                        elif affinity >= temp_aff_thresh and not min_stop == 0:
                            selection_p = 0.8 / min_stop
                        else:
                            selection_p = 0.1 / min_stop
                        
                        if random.random() > selection_p:
                            i += 1
                            continue

                        if num_mutations == 0:
                            parent_abd += abundance
                        else:
                            child = {
                                "id": self._generate_bcr_id(),
                                "sequence": mutated_seq,
                                "a_sequence":mutated_a_seq,
                                "generation": generation,
                                "parent": parent["id"],
                                "affinity": affinity,
                                "abundance": abundance,
                                "mutations": num_mutations
                            }
                            
                            all_bcrs.append(child)
                            next_generation.append(child)
                    i += 1
                parent["abundance"] = parent_abd

            if not next_generation:
                break

            current_generation = next_generation
            generation += 1

        merged_all_bcrs = merge_sequences(all_bcrs)
        return all_bcrs, merged_all_bcrs

def calculate_stats(bcrs):
    affinities = [b['affinity'] for b in bcrs]
    generations = [b['generation'] for b in bcrs]
    abundances = [b['abundance'] for b in bcrs]
    unique_seqs = set(b['sequence'] for b in bcrs)
    return {
        "naive_affinity": bcrs[0]['affinity'],
        "mean_affinity": np.mean(affinities),
        "median_affinity": np.median(affinities),
        "max_affinity": max(affinities),
        "min_affinity": min(affinities),
        "mean_abundance": np.mean(abundances),
        "max_abundance": max(abundances),
        "min_abundance": min(abundances),
        "total_abundance": sum(abundances),
        "max_generation": max(generations),
        "num_bcrs": len(bcrs),
        "num_unique_sequences": len(unique_seqs)
    }

def visualize_newick_tree(newick_file, output_file="lineage_tree_from_newick.png"):
    tree = Phylo.read(newick_file, "newick")
    fig = plt.figure(figsize=(12, 8))
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=axes)
    plt.title("BCR Lineage Tree from Newick")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


def plot_newick_bcellTree(newick, node_weights, output_file, title_name = 'Hi' ,title_num = 0):
    def rectangular_layout(tree):
        positions = {}
        x_offset = 0
        y_offset = 0
        level_spacing = 50  # Vertical spacing between levels
        sibling_spacing = 100  # Horizontal spacing between siblings

        def assign_positions(node, x, y):
            nonlocal x_offset
            if node.is_leaf():
                positions[node.name] = (x_offset, y)
                x_offset += sibling_spacing
            else:
                child_positions = []
                for child in node.children:
                    assign_positions(child, x, y - level_spacing)
                    child_positions.append(positions[child.name][0])
                positions[node.name] = (sum(child_positions) / len(child_positions), y)

        assign_positions(tree, x_offset, y_offset)
        return positions
    node_info = {}
    for key, val in node_weights.items():
        node_info[key] = f'Node_name: {key}\n\nAbundancy: {int(val)}'
    t = Tree(newick, format=1)
    tree_nodes = {node.name for node in t.traverse() if node.name}
    node_info.update({key: 'Node_name: ''\n\nAbundancy: 1' for key in tree_nodes if key not in node_info})
    node_weights.update({key: 1 for key in tree_nodes if key not in node_weights})
    max_weight = max(node_weights.values())
    node_sizes = {node: (weight / max_weight) * 1000 for node, weight in node_weights.items()}
    weights = np.array(list(node_weights.values()))
    normalized_weights = (weights - weights.min()) / (weights.max() - weights.min())  # Scale between 0 and 1
    colors = cm.viridis(normalized_weights)
    
    # Create a mapping of nodes to colors
    node_colors = {node: colors[idx] for idx, node in enumerate(node_weights.keys())}
    pos = rectangular_layout(t)
    
    # Draw the tree with rectangular edges
    fig = plt.figure(figsize=(8, 5))
    fig.canvas.manager.set_window_title(f'{title_name} T{title_num+1}')

    x_coords = [x for x, y in pos.values()]
    y_coords = [y for x, y in pos.values()]
    x_range = max(x_coords) - min(x_coords)
    y_range = max(y_coords) - min(y_coords)
    x_margin = x_range * 0.3 if x_range > 0 else 1.0
    y_margin = y_range * 0.3 if y_range > 0 else 1.0
    plt.gca().set_xlim(min(x_coords) - x_margin, max(x_coords) + x_margin)
    plt.gca().set_ylim(min(y_coords) - y_margin, max(y_coords) + y_margin)
    
    # Draw nodes
    scatter = plt.scatter([], [], s=[], alpha=0.9, color=[], edgecolor="black", zorder=2)
    scatter.set_offsets([pos[node] for node in tree_nodes])
    scatter.set_sizes([node_sizes[node] for node in tree_nodes])
    scatter.set_color([node_colors[node] for node in tree_nodes])

    for i, node in enumerate(tree_nodes):
        if node == "naive":
            plt.scatter(
                pos[node][0], pos[node][1],
                s=20,
                alpha=0.9,
                color="black",
                edgecolor="black",
                zorder=2, 
                marker="^" 
            )
    
    for node in t.traverse("postorder"):
        if not node.is_root():
            parent = node.up
            x_start, y_start = pos[parent.name]
            x_end, y_end = pos[node.name]
            plt.plot([x_start, x_end], [y_start, y_start], color="black", lw=1, zorder=1)
            plt.plot([x_end, x_end], [y_start, y_end], color="black", lw=1, zorder=1)


    plt.axis("off")
    
    cursor = mplcursors.cursor(scatter, hover=True)

    @cursor.connect("add")
    def on_add(sel):
        node_name = list(tree_nodes)[sel.index]
        info_text = node_info.get(node_name, "No information available")
        # width = 20  # Adjust the width as per your preference
        # centered_text = "\n".join(line.center(width) for line in info_text.splitlines())
        # sel.annotation.set_text(centered_text)
        sel.annotation.set_text(info_text)
        sel.annotation.set_multialignment('left') 
        bbox_properties = dict(
        boxstyle="round,pad=0.7",
        edgecolor="black",
        facecolor="yellow",
        linewidth=1,
        alpha=0.7
        )
        sel.annotation.set_bbox(bbox_properties)
        sel.annotation.arrowprops = None
        sel.annotation.get_bbox_patch()
        
        @sel.annotation.axes.figure.canvas.mpl_disconnect
        def remove_annotation(event):
            sel.annotation.set_visible(False)
            sel.annotation.axes.figure.canvas.draw_idle()
    plt.savefig(output_file, format='png')
    plt.close()

def export_to_fasta(bcrs, filename):
    with open(filename, 'w') as f:
        for b in bcrs:
            f.write(f">{b['id']}@{b['abundance']}\n{b['sequence']}\n")

def export_to_newick(bcrs, filename):
    from collections import defaultdict
    tree = defaultdict(list)
    seq_map = {}
    parent_map = {}
    
    for b in bcrs:
        b_id = b['id']
        parent = b['parent']
        if parent:
            tree[parent].append(b_id)
        seq_map[b_id] = b['sequence']
        parent_map[b_id] = parent

    roots = [b['id'] for b in bcrs if b['parent'] is None]
    if not roots:
        return

    def mutation_distance(seq1, seq2):
        return sum(a != b for a, b in zip(seq1, seq2))

    def recurse(node_id):
        children = tree[node_id]
        label = f"{node_id}"
        if parent_map[node_id] is None:
            label = "naive"
            dist = 0
        else:
            parent_seq = seq_map[parent_map[node_id]]
            child_seq = seq_map[node_id]
            dist = mutation_distance(parent_seq, child_seq)
        if not children:
            return f"{label}:{dist}"
        else:
            subtree = ",".join(recurse(c) for c in children)
            return f"({subtree}){label}:{dist}"

    newick = '('+recurse(roots[0])+')'+';'
    with open(filename, 'w') as f:
        f.write(newick)
    return newick

def translate_nucleotide_to_protein(nuc_seq):
    seq_obj = Seq(nuc_seq.upper().replace('U', 'T'))
    longest_protein = ""
    for frame in range(3):
        protein_seq = seq_obj[frame:].translate(to_stop=True)
        protein_str = str(protein_seq)
        if len(protein_str) > len(longest_protein):
            longest_protein = protein_str
    return longest_protein

def translate_nucleotide_to_protein_min_stops(nuc_seq, regions):
    seq_obj = Seq(nuc_seq.upper().replace('U', 'T'))
    best_protein = None
    min_stops = None
    for frame in range(3):
        protein_seq = seq_obj[frame:].translate(to_stop=False)
        protein_str = str(protein_seq)
        total_stops = 0
        if regions:
            for start_nt, end_nt in regions:
                aa_start = max(0, (start_nt - frame + 2) // 3)
                aa_end = max(0, (end_nt - frame) // 3)
                if aa_start <= aa_end and aa_end <= len(protein_str):
                    subregion = protein_str[aa_start:aa_end + 1]
                else:
                    subregion = ""
                total_stops += subregion.count('*')
        else:
            total_stops = protein_seq.count('*')
        if min_stops is None or total_stops < min_stops:
            min_stops = total_stops
            best_protein = protein_str
    return best_protein, min_stops


def group_by_sequence(seq_list):

    def extract_seq_num(seq_id):
        match = re.search(r'(\d+)', seq_id)
        return int(match.group(1)) if match else float('inf')
    
    seq_to_ids = defaultdict(list)
    for entry in seq_list:
        seq_to_ids[entry['sequence']].append(entry['id'])

    result = {}
    for seq, ids in seq_to_ids.items():
        sorted_ids = sorted(ids, key=extract_seq_num)
        key_id = sorted_ids[0]
        result[key_id] = sorted_ids
    return result


from collections import defaultdict

def merge_sequences(seq_list):
    grouped = defaultdict(list)
    for entry in seq_list:
        key = (entry['sequence'], entry['parent'])
        grouped[key].append(entry)

    merged_map = {} 
    merged_list = []
    for entries in grouped.values():
        if len(entries) == 1:
            merged_list.append(entries[0])
        else:
            entries.sort(key=lambda x: int(x['id'][3:]))
            merged_id = entries[0]['id']
            total_abundance = sum(e['abundance'] for e in entries)
            merged_entry = entries[0].copy()
            merged_entry['id'] = merged_id
            merged_entry['abundance'] = total_abundance
            merged_list.append(merged_entry)
            for e in entries:
                merged_map[e['id']] = merged_id

    for entry in merged_list:
        parent = entry['parent']
        if parent in merged_map:
            entry['parent'] = merged_map[parent]
    for entry in merged_list:
        if entry['id'] not in merged_map:
            merged_map[entry['id']] = entry['id'] 

    for entry in merged_list:
        if entry['parent'] in merged_map:
            entry['parent'] = merged_map[entry['parent']]

    return merged_list

IUPAC_CODES = {
    'A': {'a'}, 'C': {'c'}, 'G': {'g'}, 'T': {'t'},
    'R': {'a', 'g'}, 'Y': {'c', 't'}, 'W': {'a', 't'},
    'S': {'g', 'c'}, 'K': {'g', 't'}, 'M': {'a', 'c'},
    'B': {'c', 'g', 't'}, 'D': {'a', 'g', 't'}, 'H': {'a', 'c', 't'},
    'V': {'a', 'c', 'g'}, 'N': {'a', 'c', 'g', 't'}
}

def biased_mutation(base):
    transitions = {'a': 'g', 'g': 'a', 'c': 't', 't': 'c'}
    transversions = {
        'a': ['c', 't'],
        'g': ['c', 't'],
        'c': ['a', 'g'],
        't': ['a', 'g']
    }
    base = base.lower()
    if random.random() < 0.7:
        return transitions.get(base)
    else:
        return random.choice(transversions.get(base))

def read_fasta(filepath):
            return [str(record.seq) for record in SeqIO.parse(filepath, "fasta")]

def preproc_initial_bcr(vdj_folder_path):
    segments = {
            "V": read_fasta(os.path.join(vdj_folder_path, 'V.fasta')),
            "D": read_fasta(os.path.join(vdj_folder_path, 'D.fasta')),
            "J": read_fasta(os.path.join(vdj_folder_path, 'J.fasta'))
            }
    v_segments = segments["V"]
    d_segments = segments["D"]
    j_segments = segments["J"]

    v_translations = [translate_nucleotide_to_protein_min_stops(v, [])[0] for v in v_segments]
    d_translations = [translate_nucleotide_to_protein_min_stops(d, [])[0] for d in d_segments]
    j_translations = [translate_nucleotide_to_protein_min_stops(j, [])[0] for j in j_segments]

    combos = [
        (v_idx, d_idx, j_idx, v_a + d_a + j_a)
        for (v_idx, v_a), (d_idx, d_a), (j_idx, j_a) in product(
            enumerate(v_translations),
            enumerate(d_translations),
            enumerate(j_translations)
        )
    ]

    def compute_score(combo):
        v_idx, d_idx, j_idx, combined_seq = combo
        score = pathogen.calculate_affinity_a(combined_seq)
        return (v_idx, d_idx, j_idx, score)

    with ThreadPoolExecutor() as executor:
        naive_bcrs = list(executor.map(compute_score, combos))

    return naive_bcrs

if __name__ == "__main__":
    sim = BCRSimulator()
    
    with open(r'K:\BCellPhylo\ClusteringBCtrees\Bcr_simulator\IGH genes\epitope.txt', 'r') as file:
        epitope_list = [line.strip() for line in file.readlines()]
    

    # random.seed(50)
    # n_epitopes = random.randint(20, 40)
    # a_epitopes = random.sample(epitope_list, n_epitopes)
    # pathogen = Pathogen("patho", n_epitopes=[], a_epitopes=a_epitopes)
    # random.seed(None)

    pathogen = Pathogen("patho", n_epitopes=[], a_epitopes=epitope_list)

    affinity_threshold = 0.8
    max_generations = 12
    mutation_rate = 0.001
    n_runs = 1
    min_uniq_seq = 15
    prec = 0# 1 or 0       1 if preproccessing is True, 0 if not
    base_output_dir = os.getcwd()
    vdj_files_folder = 'IGH genes'

    if prec == 1:
        vdj_files_folder_path = os.path.join(base_output_dir, vdj_files_folder)
        preproc_genes_index = preproc_initial_bcr(vdj_files_folder_path)
    else:
        preproc_genes_index=[]
    current_run = 1
    while(current_run <= n_runs):
        print(f"\n--- RUN {current_run} ---")
        run_folder = os.path.join(base_output_dir, f"run_{current_run}")
        os.makedirs(run_folder, exist_ok=True)

        repertoire, merged_repertoire = sim.generate_repertoire_dynamically(
            pathogen,
            affinity_threshold,
            max_generations,
            mutation_rate, 
            preproc_genes_index
        )

        #stats = calculate_stats(repertoire)
        stats = calculate_stats(merged_repertoire)
        if stats['num_unique_sequences'] >= min_uniq_seq and np.abs(stats['max_affinity'] - affinity_threshold) <= 0.1:
            print('Repertoire was successfully produced.')
            with open(os.path.join(run_folder, "Statistics.txt"), 'w') as file:
                for key, value in stats.items():
                    file.write(f"{key}: {value}\n")
            export_to_fasta(merged_repertoire, os.path.join(run_folder, "repertoire.fasta"))
            newick = export_to_newick(merged_repertoire, os.path.join(run_folder, "repertoire.nk"))
            #visualize_newick_tree(os.path.join(run_folder, "repertoire.newick"), os.path.join(run_folder, "lineage_tree.png"))
            node_weights = {item['id']: item['abundance'] for item in merged_repertoire}
            plot_newick_bcellTree(newick, node_weights, os.path.join(run_folder, "lineage_tree.png"), title_name = f"RUN {current_run}" ,title_num = current_run)
            keys_to_include = ['id', 'generation', 'parent', 'affinity', 'abundance', 'mutations']
            with open(os.path.join(run_folder, "repertoire_info.csv"), 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=keys_to_include)
                writer.writeheader()
                for entry in merged_repertoire:
                    filtered_entry = {k: entry[k] for k in keys_to_include if k in entry}
                    writer.writerow(filtered_entry)

            if stats['num_unique_sequences'] != stats['num_bcrs']:
                unique_seq_list = group_by_sequence(merged_repertoire)
                with open(os.path.join(run_folder, "duplicated_seqs.txt"), 'w') as file:
                    for key, value in unique_seq_list.items():
                        if len(value) > 1:
                            file.write(f"{key}: {value}\n")
            print(f"\nSaved outputs in: {run_folder}")
            current_run += 1
        