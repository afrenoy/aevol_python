aevol_python
============

Python code to analyze data from aevol

## Curves.py

Produce curves with statistical data for each experiment (experiment = set of replicates with different seeds).

Usage:

    ./curves.py [-a] [-s start_gen] [-e end_gen] [-g genetic_unit] exp1 exp2 .. expn

    -a: write one curve per replicate instead of the average curve with standard error
    -s: start generation
    -e: end generation
    -g: which genetic units to treat (0: all (default), 1: chromosome, 2: plasmid)

## Movies.py

Produces a movie for any property written in a dump file.

Usage:

    ./movies.py replicatefolder

Where repliactefolder contains one replicate, with dumps each generation.
Everything (options, property to record) is hardcoded inside the file for now.

