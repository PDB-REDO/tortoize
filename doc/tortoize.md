# NAME

tortoize - Calculate ramachandran z-scores

# SYNOPSIS

tortoize \[OPTION\] input \[output\]

# DESCRIPTION

Tortoize validates protein structure models by checking the Ramachandran
plot and side-chain rotamer distributions. Quality Z-scores are given at
the residue level and at the model level (ramachandran-z and torsion-z).
Higher scores are better. To compare models or to describe the
reliability of the model Z-scores jackknife- based standard deviations
are also reported (ramachandran-jackknife-sd and torsion-jackknife-sd).

# OPTIONS

The input file can be either mmCIF or PDB format and the file may be
gzip compressed.

The output is a json file, if no file name is specified the output is
written to *stdout*.

**-h**, **\--help**

:   Display a help message and exit.

**\--version**

:   Print the version and exit.

**-v**, **\--verbose**

:   Increase the verbosity of the diagnostic output.

**\--dict**=\<file\>

:   Specify a dictionary file containing restraints for residues
    specific to this file. This option may be specified more than once.

**\--log**=\<file\>

:   Write a log with diagnostic information to this file. An output file
    must also be specified when this option is used.

# REFERENCES

References:

**Sobolev et al.**

:   A Global Ramachandran Score Identifies Protein Structures with
    Unlikely Stereochemistry, Structure (2020), DOI:
    *https://doi.org/10.1016/j.str.2020.08.005*

**Van Beusekom et al.**

:   Homology-based loop modeling yields more complete crystallographic
    protein structures, IUCrJ (2018), DOI:
    *https://doi.org/10.1107/S2052252518010552*

**Hooft et al.**

:   Objectively judging the quality of a protein structure from a
    Ramachandran plot, CABIOS (1993), DOI:
    *https://doi.org/10.1093/bioinformatics/13.4.425*

# AUTHOR

Written by Maarten L. Hekkelman \<maarten@hekkelman.com\>

# REPORTING BUGS

Report bugs at https://github.com/PDB-REDO/tortoize/issues
