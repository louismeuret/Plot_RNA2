
def get_contact_map_and_signature(
    native_structure_file_name,
    atoms_selection='nucleicbackbone',
    min_residues_apart=3,
    distance_scale=7.5,
    atom_index_file='atom_indices.txt',
    heatmap_file='contact_map_heatmap.png'
):
    """
    Computes the contact map and signature of the given native structure file, writes atom correspondence, logs debugging information, and generates a heatmap.
    
    Parameters
    ----------
    native_structure_file_name: str
        A pdb/gro file containing the native structure.
    atoms_selection: str
        A string specifying the atoms to keep in computation. Default is 'nucleicbackbone'.
    min_residues_apart: int
        The minimum number of residues that must separate atom pairs in the output signature. Default is 3.
    distance_scale: float
        The characteristic distance scale (Angstroms). Default is 7.5.
    atom_index_file: str
        Name of the file to save atom-index correspondence.
    heatmap_file: str
        Name of the file to save the contact map heatmap.
    
    Returns
    -------
    contact_map: np.ndarray
        The contact map of the native structure.
    signature: np.ndarray
        The signature of atom-atom contacts in the native structure.
    """

    # Step 1: Load structure with MDAnalysis
    u = mda.Universe(native_structure_file_name)

    # Step 2: Select backbone atoms in MDAnalysis
    backbone_selection = u.select_atoms(atoms_selection)
    
    # Step 3: Get the indices of the selected backbone atoms
    backbone_indices = backbone_selection.indices
    print(f"Number of backbone atoms: {len(backbone_indices)}")
    
    # Step 4: Load the same structure into MDTraj
    native_structure = md.load(native_structure_file_name)
    
    # Step 5: Write the correspondence of atom indices
    with open(atom_index_file, 'w') as f:
        for idx in backbone_indices:
            atom = native_structure.topology.atom(idx)
            f.write(f"{idx}: {atom}\n")
    
    # Step 6: Map residue indices
    residue_indices = [
        atom.residue.index for atom in native_structure.topology.atoms
        if atom.index in backbone_indices
    ]
    
    # Step 7: Generate signature based on minimum residue separation
    signature = []
    for i in range(len(backbone_indices)):
        for j in range(i + 1, len(backbone_indices)):
            if residue_indices[j] - residue_indices[i] >= min_residues_apart:
                signature.append([backbone_indices[i], backbone_indices[j]])
    signature = np.array(signature)
    
    # Step 8: Compute distances
    distances = md.compute_distances(native_structure, signature)[0]
    distances = distances / (distance_scale / 10.0)  # Convert Angstroms to nm
    
    # Step 9: Compute the contact map
    contact_map = np.zeros(distances.shape)
    where_zero = distances == 0
    where_flexus = np.abs(distances - 1) < 1e-3
    contact_map[where_zero] = 1.0
    contact_map[where_flexus] = 0.6
    mask = ~ (where_zero | where_flexus)
    contact_map[mask] = (1 - distances[mask] ** 6) / (1 - distances[mask] ** 10)
    
    # Step 10: Generate heatmap
    heatmap_matrix = np.zeros((len(backbone_indices), len(backbone_indices)))
    for idx, (i, j) in enumerate(signature):
        heatmap_matrix[i, j] = contact_map[idx]
        heatmap_matrix[j, i] = contact_map[idx]
    
    plt.figure(figsize=(10, 8))
    plt.title("Contact Map Heatmap")
    plt.imshow(heatmap_matrix, cmap='viridis', interpolation='nearest')
    plt.colorbar(label="Contact Strength")
    plt.xlabel("Atom Index")
    plt.ylabel("Atom Index")
    plt.savefig(heatmap_file)
    print(f"Heatmap saved to {heatmap_file}")
    print(contact_map.shape)
    
    return contact_map, signature
