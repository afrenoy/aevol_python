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

    movies.py -i inputname [-s firstgen] [-e lastgen] [-p patchsize] [-m minvalue] [-M maxvalue] [-c color] [-k] simulation

Where 'simulation' is the path of a simulation folder, with dumps for each generation, and 'inputname' is the name of the property we want to use in the movie, for example 'fitness_metabolic' (check the avaibility of files stats/dump/intputname_nnnn.out).

Example:

    ./movies.py -i fitness_metabolic -p 40 -c '255 255 255' -m 0 -M 1 -s 30 -e 100 ~/myaevolsimulation

See source code for documentation of other parameters.
 
