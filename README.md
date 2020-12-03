# mutation_report
Genetic analysis tool written for C-Path.

## Purpose
Critical Path Institute is involved in a world-wide effort to identify drugs which are most effective against
Tuberculosis. C-Path gathers TB samples from partner clinics, along with drug toxicity reports. The samples
are sequenced and then run through a pipeline that identifies specific mutations in the sample's DNA and matches
the sample with its drug reports.

These tools were written to take several of the files produced by the pipeline and extract relevent data for
presentation in a more readable format.

## Structure
There are two main programs in this package. One extracts the data of interest and creates a JSON file of all
the pertinent information. The second reads the JSON file and can create an output file as either a PDF report
or a CSV file, depending on command line options. The program that creates the JSON file has several components
to extract specific types of information. These were created as seperate modules so that they might be used
in other projects in the future.
