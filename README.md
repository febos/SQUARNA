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

# Output format

# Options 

    i=FILENAME / input=FILENAME [REQUIRED OPTION]
    
        Path to an input file in fasta-like format, 
        see "Input format" section for details.

	c=FILENAME / config=FILENAME [DEFAULT: c=./def.conf]
    
        Path to a config file, see files "def.conf" 
        and "alt.conf" for the format details.



    inputformat = "qtrf"               # Input line order, q=seQuence,t=reacTivities,r=Restraints,f=reFerence

    maxstemnumset = False              # do we overwrite the maxstemnum from configfile
    maxstemnum    = 10**6              # maximum number of stems for each structure

    threads       = os.cpu_count()     # Number of cpus to use

    rankbydiff     = False             # Output diverse structures first
    rankby         = "rs"              # Rank by, r / s / rs / dr / ds / drs, r=reactscore,s=structscore,d=rankbydiff
    hardrest       = False             # Force bp-restraints into predicted structures 
    interchainonly = False             # Forbid intra-chain base pairs

    toplim        = 5                  # Top-N to print
    outplim       = toplim             # Top-N structs used for metrics calculations if reference
    conslim       = 1                  # Top-N structs used for consensus

    reactformat   = 3                  # 3 / 10 / 26


    q=FILENAME [REQUIRED OPTION]
        Path to a query structure, the one that ARTEM superimposes to 
        the reference, in PDB/mmCIF format.

    matchrange=FLOAT [DEFAULT: matchrange=3.0]
    	The threshold used for searching the mutually closest residues. Only 
    	those pairs of residues that have their centers of mass at a distance 
    	under the specified matchrange value can be added to the subset 
    	of the mutually closest residues. The higher matchrange value 
    	will produce more "noisy" matchings but won't miss anything. The lower 
    	matchrange value will produce more "clean" matchings but 
    	can miss something.

    rformat=KEYWORD, qformat=KEYWORD [DEFAULT: rformat=PDB,qformat=PDB] 
        The specification of the input coordinate file formats 
        (case-insensitive). By default, ARTEM tries to infer the format 
        from the extensions of the input filenames. ".pdb", ".cif", 
        and ".mmcif" formats can be recognized (case-insensitive). In the case 
        of any other extension ARTEM will treat the file as the PDB-format 
        file by default. If the "rformat" ("qformat") parameter is specified 
        and it's either "PDB", "CIF", or "MMCIF" (case-insensitive), 
        ARTEM will treat the reference (query) coordinate file
        as the specified format.

    rmsdmin=FLOAT, rmsdmax=FLOAT [DEFAULT: rmsdmin=0,rmsdmax=1e10]
    	The specification of minimum and maximum RMSD thresholds. 
    	ARTEM will print and save only the superpositions that satisfy 
    	the specified thresholds. 
    	
    rmsdsizemin=FLOAT, rmsdsizemax=FLOAT [DEFAULT: rmsdsizemin=0,rmsdsizemax=1e10]
        The specification of minimum and maximum RMSD/SIZE ratio thresholds. 
        ARTEM will print and save only the superpositions that satisfy 
        the specified thresholds. 

    resrmsdmin=FLOAT, resrmsdmax=FLOAT [DEFAULT: resrmsdmin=0, resrmsdmax=1e10]
    	The specification of minimum and maximum per-residue RMSD threshold.
    	ARTEM will print and save only the superpositions that satisfy 
    	the specified thresholds.

    rres=STRING, qres=STRING [DEFAULT: rres="#1",qres="#1"]
        The specification of the input reference (rres) and query (qres) 
        structures. Only the specified residues will considered as part 
        of the structure and all the other residues will be ignored. 
        See the format description at the end of the OPTIONS section.

    rresneg=STRING, qresneg=STRING [DEFAULT: None]
    	The specification of the input reference (rresneg) and query (qresneg) 
    	structures. The specified residues will be ignored and all the other 
    	residues considered as part of the structure. If both "rres" 
    	and "rresneg" (or "qres" and "qresneg") are specified simultaneusly, 
    	ARTEM will ignore "rres" ("qres") and consider only "rresneg" 
    	("qresneg").
        See the format description at the end of the OPTIONS section.

    rseed=STRING, qseed=STRING [DEFAULT: rseed=rres,qseed=qres]
        The specification of the reference and query residues that ARTEM can use
        for single-residue matching seeds.
        See the format description at the end of the OPTIONS section.

    saveformat=KEYWORD [DEFAULT: saveformat=qformat]
        The specification of the format of the output coordinate files. 
        By default, ARTEM will save the coordinate files in the same format 
        as the query input file. If the "saveformat" parameter is specified 
        and it's either "PDB", "CIF", or "MMCIF" (case-insensitive), ARTEM 
        will save the output coordinate files in the specified format.

    saveres=STRING [DEFAULT: saveres=qres]
    	The specification of the query structure residues that will be saved 
    	in the output coordinate files.
        See the format description at the end of the OPTIONS section.

    saveto=FOLDER [DEFAULT: None]
        Path to the output folder to save the coordinate files 
        of the superimposed query structures along with the mutually 
        closest residue subsets. If the specified folder does not exist, 
        ARTEM will create it. If the folder is not specified, 
        nothing will be saved.

    sizemin=FLOAT, sizemax=FLOAT [DEFAULT: sizemin=1,sizemax=1e10]
        The specification of minimum and maximum SIZE thresholds. 
        ARTEM will print and save only the superpositions that satisfy 
        the specified thresholds. If sizemin is set to zero, ARTEM will 
        output the dummy 0-size superpositions matching the reference 
        and query residues lacking the 5-atom representation specified. 
        The user can specify custom atom representations of any residues 
        via editing the lib/nar.py text file.

    trim=BOOL [DEFAULT: trim=None]
    	When specified, for each subset of mutually closest residues ARTEM will 
    	iteratively remove residues with the worst per-residue RMSD from the 
    	subset one by one with the following re-superpositioning based on the 
    	remaining residues of the subset until the specified thresholds for rmsdmax,
    	rmsdsizemax, resrmsdmax are reached or the subset size is less than sizemin.

    t=INT / threads=INT [DEFAULT: t=cpu_count]
    
        Number of CPUs to use. 


# Contacts

Eugene F. Baulin, *e-mail: ebaulin@iimcb.gov.pl, efbaulin@gmail.com* 
