# SQUARNA

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

# Input format

    For input SQUARNA uses a fasta-like format with the "name" lines
    starting with ">" sign and the following lines treated as the data
    lines. The order of lines in which SQUARNA will read the data 
    is defined by the inputformat (if) parameter, see below. By default,
    the order is "qtrf", meaning the first line will be read as the
    seQuence, the second line as the reacTivities line, the third line
    as the Restraints line, the fourth line as the reFerence structure 
    line, and all the following lines will be ignored until the new 
    "name" line. 
    
    Sequence is the only mandatory field. Symbols "AaCcGgUu" will be 
    treated as the four types of RNA bases in the case-insensitive
    manner. Symbols "Tt" will be replaced with U. Symbols ";&" will be
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

# Output format

    The output format is a fasta-like format with the "name" lines 
    starting with ">" sign and followed by a number of data sections:
    (a) the input sequence; (b) the input data lines with the 
    appropriate marks (reactivities/restraints/reference); (c) a 
    break line of underscores ("_"-line); (d) the predicted consensus
    structure with the mark top-X_consensus, where X is defined with 
    conslim parameter, see below. If a reference was specified the 
    metrics values will be printed in the same line (TP - number of 
    correctly predicted base pairs; FP - number of wrongly predicted
    base pairs; FN - number of missed base pairs; FS - F-score; PR - 
    precision; RC=recall); (e) a break line of equality signs ("="-line);
    (f) N lines with the predicted structures, where N is defined with
    outplim parameter, see below. The structures are followed by a 
    tab-separated list of values: the rank of the structure (starting
    with "#" sign), total_score, structure_score, reactivity_score,
    name of the generative parameter set, and (if a reference was 
    specified) the metrics values will be printed for the best of
    the top-K structures (the format is the same as for the consensus 
    structure with the only addition of RK (rank) value), where K is 
    defined with toplim parameter, see below. The chain separators are 
    introduced into all the lines as they appear in the sequence.

# Options 

    i=FILENAME / input=FILENAME [REQUIRED OPTION]
    
        Path to an input file in fasta-like format, see "Input format" 
        section for details.

    c=FILENAME / config=FILENAME [DEFAULT: c=./def.conf]
    
        Path to a config file, see files "def.conf" and "alt.conf" 
        for the format details.

    if={qtrf} / inputformat={qtrf} [DEFAULT: if=qtrf]
    
        The order of the lines in the input file. By default, SQUARNA 
        reads the first line of the entry (among the lines after 
        the ">" line) as the seQuence (q), the second line as the 
        reacTivities (t), the third line as the Restraints (r), 
        the fourth line as the reFerence (f), and all the further lines 
        are ignored. inputformat should be a subset of qtrf letters 
        in any order, with q being mandatory.
  
    rb={rsd} / rankby={rsd} [DEFAULT: rb=rs]
    
        How to rank the predicted structures. rankby should be a subset of
        letters r, s, and d in any order (r / s / rs / rd / sd / rsd).
        If both r and s are present, the structures will be ranked according
        to the total_score = structure_score * reactivity_score. If only 
        r is present, the structures will be ranked by the reactivity_score,
        and if only s is present, the structures will be ranked by the 
        structure_score. Independently, if d is present, the mutually 
        divergent structures will be put first.  
        
    tl=INT / toplim=INT [DEFAULT: tl=5]
    
        How many top-N structures will be subject to comparison with the
        reference.
        
    ol=INT / outplim=INT [DEFAULT: ol=tl]
    
        How many top-N structures will be printed into the stdout.
        
    cl=INT / conslim=INT [DEFAULT: cl=1]
    
        How many top-N structures will be used to derive the predicted
        structure consensus.

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


# Contacts

Eugene F. Baulin, *e-mail: ebaulin@iimcb.gov.pl, efbaulin@gmail.com* 
