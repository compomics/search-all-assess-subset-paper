The Mzid to Tsv Converter converts an .mzid file created by MS-GF+ to a 
tab-delimited text file.  Although MS-GF+ has this option (see MzidToTsv.html 
at http://htmlpreview.github.io/?https://github.com/sangtaekim/msgfplus/blob/master/doc/MzidToTsv.html )
MzidToTsvConverter.exe can convert the .mzid file faster, using less memory.

Usage: MzidToTsvConverter -mzid:"mzid path" [-tsv:"tsv output path"] [-unroll|-u] [-showDecoy|-sd]

Required parameters:
    -mzid:path 
        - Path to the .mzid or .mzid.gz file.  If the path has spaces, it must be in quotes.

Optional parameters:
    -tsv:path
        - Path to the tsv file to be written. If not specified, will be created in the same location as the .mzid file.
    -unroll or -u
        - Signifies that results should be unrolled: one line per unique peptide/protein combination in each spectrum identification
    -showDecoy or -sd
        - Signifies that decoy results should be included in the output .tsv file.
        - Decoy results have protein names that start with XXX_
