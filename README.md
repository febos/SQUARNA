# SQUARNA, version 0.9 [01.09.2023]

[D.R. Bohdan, G.I. Nikolaev, J.M. Bujnicki, E.F. Baulin (2023) SQUARNA - an RNA secondary structure prediction method based on a greedy stem formation model. bioRxiv. DOI: 10.1101/2023.08.28.555103](https://doi.org/10.1101/2023.08.28.555103)

# Usage

pathto/python3 pathto/SQUARNA.py i=inputfile [OPTIONS]

# Installation & Dependencies

1) Clone the GitHub repository by typing:

	git clone https://github.com/febos/SQUARNA
	
2) For C++ version compile the code by typing:

	TODO
	
SQUARNA requires Python3 of at least 3.8 version along with
a (hopefully) arbitrary version of NumPy library. 

# Usage examples

    1) python3 SQUARNA.py i=examples/seq_input.fas c=alt.conf
    
    Demonstration example.
    
    2) python3 SQUARNA.py i=datasets/SRtrain150.fas if=qf 
    
    An example reproducing the benchmarks.
    
    3) python3 SQUARNA.py i=examples/ali_input.afa if=q
    
    An example of running single-sequence predictions for a set
    of aligned sequences. "if=q" tells SQUARNA to ignore all the
    default input lines and read only the sequences.
    
    4) python3 SQUARNA.py i=examples/ali_input.afa a
    
    An example of running alignment-based predictions.
    
    5) python3 SQUARNA.py i=examples/ali_input.afa if=q a v
    
    An example of running alignment-based predictions in the
    verbose mode.

# Input format

    For input SQUARNA uses a fasta-like format with the "name" lines
    starting with ">" symbol and the following lines treated as the data
    lines. The order of lines in which SQUARNA will read the data 
    is defined by the inputformat (if) parameter, see below. By default,
    the order is "qtrf", meaning the first line will be read as the
    seQuence, the second line as the reacTivities line, the third line
    as the Restraints line, the fourth line as the reFerence structure 
    line, and all the following lines will be ignored until the new 
    "name" line. 
    
    The starting lines in the input file faced before the first ">"
    symbol will be treated as default reactivities/restraints/reference
    lines according to the inputformat value. The default lines will be 
    used for the individual sequences of appropriate length if the 
    matching individual line is empty or absent. See ali_input.afa file
    in the examples sub-folder for an example of default lines. 
    
    Sequence is the only mandatory field. Symbols "AaCcGgUu" will be 
    treated as the four types of RNA bases in the case-insensitive
    manner. Symbols "Tt" will be replaced with "U". Symbols ";&" will be
    treated as the separators of different RNA chains. Symbols ".-~"
    will be treated as gaps and ignored accordingly (the sequence along 
    with the other data lines will be unaligned for the prediction and 
    the predicted structures will be then realigned back to the initial 
    sequence). All the other symbols will be treated as bases that cannot
    form any base pairs unless they are present in the bpweights parameter
    as defined within the used config file.
    
    Reactivities can be given either as a space/tab separated sequence 
    of float values from 0.0 to 1.0 (with the number of values equal 
    to the sequence length, including the separator positions whose 
    values will be ignored), or as an encoded line of sequence length, 
    see the description of reactformat (rf) parameter below (the mix 
    of the encoded values is allowed), or be an empty line. 
    
    Restraints line should be either a sequence length line or an empty
    line. All pairs of brackets ((),[],{},<>) and pairs of latin letters 
    (Aa,Bb,...Yy,Zz) will be treated as base pairs. The underscore sign
    "_" will be treated as an unpaired base. The slash sign "/" will be
    treated as the base that cannot form any base pairs "to the left" 
    (i.e. with the counterpart being closer to the 5'-end). In contrast,
    the backslash sign "\" will be treated as the base that cannot form
    any base pairs "to the right" (i.e. with the counterpart being closer
    to the 3'-end). All the other symbols will be treated as unrestrained
    positions.
    
    Reference line should be either a sequence length line or an empty 
    line. In the reference line all pairs of brackets ((),[],{},<>) and 
    pairs of latin letters (Aa,Bb,...Yy,Zz) will be treated as base pairs 
    and all the other characters will be ignored. 
    
    For examples of an appropriate default input file format see 
    examples/seq_input.fas. 
    
    Alternatively, SQUARNA can read standard Fasta, Stockholm 
    and Clustal (.aln) formats. The input format is recognized 
    automatically. In the case of Stockholm format the "SS_cons"
    structure will be treated as default reference line. In the case
    of Fasta or Clustal formats only the sequences will be processed. 

# Output format (single-sequence mode)

    The output format is a fasta-like format with the "name" lines 
    starting with ">" sign and followed by a number of data sections:
    (a) the input sequence; (b) the input data lines with the 
    appropriate marks (reactivities/restraints/reference), the scores
    for the reference structure are printed if the reference is 
    specified (if non-canonical base pairs are present in the reference
    structure, they are considered with 0.0 weight); (c) a break line 
    of underscores ("_"-line); (d) the predicted consensus structure with 
    the mark top-X_consensus, where X is defined with conslim parameter, 
    see below. If a reference was specified the metrics values will be 
    printed in the same line (TP - number of correctly predicted base pairs; 
    FP - number of wrongly predicted base pairs; FN - number of missed 
    base pairs; FS - F-score; PR - precision; RC=recall); (e) a break 
    line of equality signs ("="-line); (f) N lines with the predicted 
    structures, where N is defined with outplim parameter, see below. 
    The structures are followed by a tab-separated list of values: 
    the rank of the structure (starting with "#" sign), total_score, 
    structure_score, reactivity_score, name of the generative parameter 
    set, and (if a reference was specified) the metrics values will be 
    printed for the best of the top-K structures (the format is the same 
    as for the consensus structure with the only addition of RK (rank) 
    value), where K is defined with toplim parameter, see below. The chain 
    separators are introduced into all the lines as they appear 
    in the sequence.
    
# Output format (alignment-based mode)

    The output format consists of three main sections: (1) intermediate
    information (in verbose mode only); (2) processed default input 
    lines; (3) the three (steps 1-3) predicted secondary structures in 
    the dot-bracket format. Between the sections 1 and 2 there is 
    the first separator line ("="-line), and between the sections 
    2 and 3 there is the second separator line ("_"-line).
    In the verbose mode the intermediate information includes 
    the following: (1) ">Step 1, Iteration 1", the conserved base pairs,
    first one by one in the dot-bracket format along with their scores,
    then assembled into a number of secondary structures; (2) ">Step 1,
    Iteration 2", the conserved base pairs from the second iteration, 
    the format is the same as in the first iteration; (3) output of
    restrained single-sequence predictions, see section "Output format
    (single-sequence mode)" for more details; (4) ">Step 2, Populated
    base pairs", the base pairs from the single-sequence predictions
    listed one by one in dot-bracket format along with the numbers of 
    their source sequences; (5) ">Step 2, Consensus", the step-2 
    consensus structures built from the populated base pairs along with
    the used frequency thresholds.

# Options 

    i=FILENAME / input=FILENAME [REQUIRED OPTION]
    
        Path to an input file in fasta-like format, see "Input format" 
        section for details.

    c=FILENAME / config=FILENAME [DEFAULT: c=./def.conf]
    
        Path to a config file, see files "def.conf" and "alt.conf" 
        for the format details. In the alignment-based mode,
        the default config file is c=./ali.conf.
        def.conf is the default config.
        alt.conf is recommended for short pseudoknotted RNAs.
        500.conf is recommended for RNAs longer 500 nts.
        1000.conf is recommended for RNAs longer 1000 nts.
        sk.conf is recommended with SHAPE data input.

    a / ali / alignment [DEFAULT: FALSE]
    
        Run SQUARNA in the alignment-based mode. If specified,
        ali.conf will be used as the config file by default, 
        unless another config file is explicitly specified 
        by the user. The bpweights, minlen, and minbpscore 
        parameters for step-1 will be derived from the first 
        parameter set in the config file.

    if={qtrf} / inputformat={qtrf} [DEFAULT: if=qtrf]
    
        The order of the lines in the input file. By default, SQUARNA 
        reads the first line of the entry (among the lines after 
        the ">" line) as the seQuence (q), the second line as the 
        reacTivities (t), the third line as the Restraints (r), 
        the fourth line as the reFerence (f), and all the further lines 
        are ignored. inputformat should be a subset of qtrf letters 
        in any order, with q being mandatory.
  
    rb={rsd} / rankby={rsd} [DEFAULT: rb=s]
    
        How to rank the predicted structures. rankby should be a subset of
        letters r, s, and d in any order (r / s / rs / rd / sd / rsd).
        If both r and s are present, the structures will be ranked according
        to the total_score = structure_score * reactivity_score. If only 
        r is present, the structures will be ranked by the reactivity_score,
        and if only s is present, the structures will be ranked by the 
        structure_score. Independently, if d is present, the mutually 
        divergent structures will be put first.  
        
    fl=INT / freqlim=INT [DEFAULT: fl=0.35]
    
        Ignored in the single-sequence mode.
        The percentage of sequences required to contain a base pair,
        in order for it to be added to the predicted consensus structure
        at step-2. The consensus will include all the base pairs present
        in at least "fl" share of the sequences given that the base pair
        is not in conflict (does not share a position) with a more 
        populated base pair.

    ll=INT / levlim=INT [DEFAULT: ll=(3 - len(seq)>500)]
    
        Ignored in the single-sequence mode.
        The allowed number of pseudoknot levels. All the base pairs
        of the higher levels will be removed from the structure predicted
        at step-1 and from the structure predicted at step-2. By default, 
        ll=3 for short alignments of no more than 500 columns, 
        and ll=2 for longer alignments.
                
    tl=INT / toplim=INT [DEFAULT: tl=5]
    
        How many top-N structures will be subject to comparison with the
        reference.
        
    ol=INT / outplim=INT [DEFAULT: ol=tl]
    
        How many top-N structures will be printed into the stdout.
        
    cl=INT / conslim=INT [DEFAULT: cl=1]
    
        How many top-N structures will be used to derive the predicted
        structure consensus.
        
    pl=INT / poollim=INT [DEFAULT: pl=1000]
    
        Maximum number of structures allowed to populate the current 
        structure pool (if exceeded, no bifurcation will occur anymore).

    s3={i,u,1,2} / step3={i,u,1,2} [DEFAULT: s3=u]
    
        Ignored in the single-sequence mode.
        Defines the structure that will be printed at step-3. If s3=1,
        the structure from step-1 will be printed, and the step-2 will
        be skipped completely, meaning the prediction will be super fast.
        If s3=2, the structure from step-2 will be printed. If s3=u, the
        union of base pairs of the two structures will be printed. 
        If s3=i, the intersection of base pairs of the two structures 
        will be printed.
        
    msn=INT / maxstemnum=INT [DEFAULT: None]
    
        Maximum number of stems to predict in each structure. By default,
        maxstemnum is defined in a config file for each parameter set. 
        If specified in the command line it will overwrite the maxstemnum 
        values for all the parameter sets. 

    rf={3,10,26} / reactformat={3,10,26} [DEFAULT: rf=3]
    
        Encoding used to output the reactivities line.
        
        rf=3:   0.0  <= "_" <   1/3; 
                1/3  <= "+" <   2/3;
                2/3  <= "#" <=  1.0;
                
        rf=10:  0.0  <= "0" <   0.1;
                ....................
                0.5  <= "5" <   0.6;
                ....................
                0.9  <= "9" <=  1.0;
                
        rf=26:  0.00 <= "a" <  0.02;
                0.02 <= "b" <  0.06;
                ....................
                0.50 <= "n" <  0.54;
                ....................
                0.94 <= "y" <  0.98;
                0.98 <= "z" <= 1.00.

    eo / evalonly [DEFAULT: FALSE]
    
        Ignored in the alignment mode.
        If specified, no predictions are made and just the reference structure
        scores are returned provided the reference is specified. 
        If non-canonical base pairs are present in the reference structure, 
        they will be considered with 0.0 weight).

    hr / hardrest [DEFAULT: FALSE]
    
        If specified, all the base pairs from the restraints line will be
        forced to be present in the predicted structures. However, it will
        not affect the structure scores, as the forced base pairs won't
        contribute to the structure score unless they were predicted without
        forcing as well.
        
    ico / interchainonly [DEFAULT: FALSE]
    
        Allow only inter-chain base pairs to be predicted.  
        
    t=INT / threads=INT [DEFAULT: t=cpu_count]
    
        Number of CPUs to use. 
        
    v / verbose [DEFAULT: FALSE]
    
        Ignored in the single-sequence mode.
        Run SQUARNA in the verbose mode. 


# Contacts

Eugene F. Baulin, *e-mail: ebaulin@iimcb.gov.pl, efbaulin@gmail.com* 
