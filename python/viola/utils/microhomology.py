def get_microhomology_from_positions(
    chrom1: str,
    pos1: int,
    chrom2: str,
    pos2: int,
    strand1: str,
    strand2: str,
    fa,
    max_homlen = 200
):
    """
    get_microhomology_from_positions(chrom1, pos1, chrom2, pos2, strand1, strand2, fa)
    Infer microhomology length and sequence from SV breakends. This function depends on Biopython.
    
    Parameters
    ------------
    chrom1: str
        Chromosome of the left breakend.
    pos1: int
        1-origin position of the left breakend.
    chrom2: str
        Chromosome of the right breakend.
    pos2: int
        1-origin position of the right breakend.
    strand1: str
        Strand of left breakend. '+' or '-'.
    strand2: str
        Strand of right breakend.
    fa: Fasta
        Fasta object.
    max_homlen: int, default 200
        Maximum length of microhomology to be considered.
    
    Note
    -------
    In the arguments "pos1" and "pos2", input 1-origin position.
    This is quite important to identify accurate microhomology.
    This function requires Biopython installation in your environment.
    """
    
    # make position 0 origin
    pos1 -= 1
    pos2 -= 1
    
    if (strand1 == '+') & (strand2 == '-'):
        left_connected = fa[chrom1][pos1-max_homlen: pos1+1]
        left_extended = fa[chrom1][pos1+1: pos1+max_homlen+1]
        right_connected = fa[chrom2][pos2: pos2+max_homlen]
        right_extended = fa[chrom2][pos2-max_homlen: pos2]
    elif (strand1 == '+') & (strand2 == '+'):
        left_connected = fa[chrom1][pos1-max_homlen: pos1+1]
        left_extended = fa[chrom1][pos1+1: pos1+max_homlen+1]
        right_connected = fa[chrom2][pos2-max_homlen: pos2+1].reverse_complement()
        right_extended = fa[chrom2][pos2+1: pos2+max_homlen+1].reverse_complement()
    elif (strand1 == '-') & (strand2 == '-'):
        right_connected = fa[chrom1][pos1: pos1+max_homlen]
        right_extended = fa[chrom1][pos1-max_homlen: pos1]
        left_connected = fa[chrom2][pos2: pos2+max_homlen].reverse_complement()
        left_extended = fa[chrom2][pos2-max_homlen: pos2].reverse_complement()
    else:
        right_connected = fa[chrom1][pos1: pos1+max_homlen]
        right_extended = fa[chrom1][pos1-max_homlen: pos1]
        left_connected = fa[chrom2][pos2-max_homlen: pos2+1]
        left_extended = fa[chrom2][pos2+1: pos2+max_homlen+1]
    
    left_count = 0
    for a, b in zip(left_connected[::-1], right_extended[::-1]):
        if a.lower() == b.lower():
            left_count += 1
        else:
            break
    
    left_seq = left_connected[::-1][:left_count][::-1].upper()
            
    right_count = 0
    for a, b in zip(right_connected, left_extended):
        if a.lower() == b.lower():
            right_count += 1
        else:
            break
    
    right_seq = right_connected[:right_count].upper()
    
    homlen = left_count + right_count
    homseq = right_seq  + left_seq
    #print(left_connected[::-1])
    #print(right_extended[::-1])
    #print(right_connected)
    #print(left_extended)

    return (homlen, str(homseq))