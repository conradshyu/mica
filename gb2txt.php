#!/usr/bin/php -q
<?php
/*
 * GB2TXT.PHP
 *
 * PHP script that parse GenBank files into the plain text format
 *
 * Written by Conrad Shyu, June 5, 2004
 * shyu4751@uidaho.edu
 *
 * Initiative for Bioinformatics and Evolutionary Stuides (IBEST)
 * Department of Bioinformatics and Computational Biology (BCB)
 * Department of Biological Sciences
 * University of Idaho, Moscow, ID 83844
 *
 * Copyrights 2004. All rights reserved.
 *
 * revised on March 13, 2006
 * revised on June 17, 2020
*/

function strReverse($str) {
    return(strrev($str));
}   // find the reverse of the string

function strComplement($str) {
    return(strtr($str, "ATUGCYRSWKMBDHVN", "TAACGRYSWMKVHDBN"));
}   // find the complement of the string (DNA)

function querySafe($query) {
    $query = ereg_replace("%",  "", $query);
    $query = ereg_replace("_",  "", $query);
    $query = ereg_replace("'",  "", $query);
    $query = ereg_replace("&", " ", $query);
    $query = ereg_replace(";",  "", $query);
    $query = ereg_replace("\"", "", $query);
    $query = ereg_replace("\n", "", $query);

    return($query);
}   // to make sure the query string is compatible with Oracle SQL server

// a function to make opening a file easier
function openFile($filename, $open_mode) {
    $file_handle = fopen($filename, $open_mode);

    if (!$file_handle) {
        echo "file $filename cannot be opened\n";
        exit(1);
    }

    return($file_handle);
}   // end of openFile()

/*
 * the following functions extract information from GenBank files
 * currently only locus, accession, organism names, and sequences are extracted
*/
function trimString($str, $len) {
    if (!(strlen($str) > $len)) {
        return($str);
    }

    $str = substr($str, 0, ($len - 3));
    $str = "$str...";
    return($str);
}   // end of trimString()


/*
 * retrieve the locus information from the GenBank file
 *
 * the locus section typicall looks like this:
 *
 * ...
 * LOCUS       S000002564   1500 bp   RNA      25-MAR-2002
 *             ~~~~~~~~~~
 *             ^ information we are looking for ...
*/
function getLocus($record) {
    // loop until we find the keyword LOCUS
    for ($i = 0; $i < sizeof($record); $i++) {
        // found the LOCUS record
        if (!(strncmp($record[ $i ], 'LOCUS', 5))) {
            $locus = strtok($record[ $i ], " \n");
            $locus = trim(strtok(" \n"));

            return($locus);
        }
    }

    // if nothing can be found ...
    return('none');
}

/*
 * retrieve the organism information from the GenBank file
 *
 * the organism section typically looks like this:
 *
 *
 * FEATURES  Location/Qualifiers
 *           source          1..200
 *           /organism="Aeromonas punctata"
 *                      ~~~~~~~~~~~~~~~~~~
 *                      ^ the information we are looking for ...
*/
function getOrganism($record) {
    // loop until we find the keyword ORGANISM
    for ($i = 0; $i < sizeof($record); $i++) {
        // found the ORGANISM record
        if (!(strncmp($record[ $i ], 'DEFINITION', 10))) {
            $organism = strtok($record[ $i ], " \n");
            $organism = trim(strtok("\n"));
            $organism = querySafe($organism);
            return(trimString($organism, 60));
        }

        // try different place
        if (($organism = strstr($record[ $i ], '/organism='))) {
            querySafe($organism);
            return(trimString($organism, 60));
        }
    }

    // if nothing can be found ...
    return('none');
}

/*
 * retrieve the organism information from the GenBank file
 *
 * the organism section typically looks like this:
 *
 * FEATURES  Location/Qualifiers
 *   source  1..200
 *           /organism="Aeromonas punctata"
 *           /strain="JCM 1060"
 *
 * NOTE: the defintion can span more than one line
*/
function getStrain($record) {
    // loop through the entire record
    for ($i = 0; $i < sizeof($record); $i++) {
        // found the organism record
        if (($strain = strstr($record[ $i ], '/strain='))) {
            return(substr($strain, 8, 68));
        }
    }

    // if nothing can be found ...
    return('none');
}

/*
 * retrieve the Genbank accession number from the GenBank file
 *
 * the accession number section typically looks like this:
 *
 * ...
 * COMMENT
 *             Corresponding GenBank entry: AF365529
 *                                          ~~~~~~~~
 *                                          ^ the information we are looking for ...
 *
 * Note: SQL doesn't like ampersand (&), it should be replaced by a space
*/
function getAccession($record) {
    for ($i = 0; $i < sizeof($record); $i++) {
        // look for accession from two different locations
        if (($accession = strstr($record[ $i ], 'GenBank entry:'))) {
            $accession = substr($accession, 15);
            $accession = strtok($accession, " |\t\n\r");
            return($accession);
        }

        if (!(strncmp($record[ $i ], 'ACCESSION', 9))) {
            $accession = strtok($record[ $i ], " \n");
            $accession = trim(strtok(" \n"));

            return($accession);
        }
    }

    return('none');
}

/*
 * retrieve the sequence from the GenBank file
 *
 * the sequence section typically looks like this:
 *
 * ...
 * ORIGIN
 *         1 CGAACGCTGG CGGCGTGCCT AATACATGCA AGTCGAGCGA AGTTTTTCTG GTGCTTGCAC
 *        61 TAGAAAAACT TAGCGGCGAA CGGGTGAGTA ACACGTAAAG AACCTGCCTC ATAGACTGGG
 *       121 ACAACTATTG GAAACGATAG CTAATACCGG ATAACAGCAT TAACTGCATG GTTGATGTTT
 *       181 GAAAGTTGGT TTTGCTAACA CTATGAGATG GCTTTGCGGT GCATTAGCTA GTTGGTGGGG
 *       241 TAAAGGCCTA CCAAGGCGAC GATGCATAGC CGACCTGAGA GGGTGATCGG CCACACTGGG
 *       301 ACTGAGACAC GGCCCAGACT CCTACGGGAG GCAGCAGTAG GGAATCTTCC GCAATGGGCG
 *       361 AAAGCCTGAC GGAGCAACGC CGCGTGAGTG AAGAAGGATT TCGGTTCGTA AAGCTCTGTT
 *       421 GTTAGGGAAG AATGATTATG TAGTAACTAT ACATAGTAGA GACGGTACCT AACCAGAAAG
 *       481 CCACGGCTAA CTACGTGCCA GCAGCCGCGG TAAT
 *
 * Note: need to check for invalid characters and replace with 'X' or 'N'
*/
function getOrigin($record) {
    $i = 0; $record_size = sizeof($record);

    for (; $i < $record_size; $i++) {
        if (!(strncmp($record[ $i ], 'ORIGIN', 6))) {
            break;
        }
    }

    $i++;           // advance the pointer to the next entry
    $origin = "";

    // now, assemble the sequences
    for (; $i < $record_size; $i++) {
        // extract each segment of sequences from the file
        $pieces = explode(" ", trim($record[ $i ]));

        for ($j = 1; $j < sizeof($pieces); $j++) {
            $origin .= $pieces[ $j ];
        }
    }

    return(strtr(strtoupper($origin), "UX", "TN"));
}

/*
 * main procedure starts here ...
*/

$cmd_param = $_SERVER[ 'argv' ];        // get the command line parameter

// we need at least  parameters
if (sizeof($cmd_param) < 3) {
    echo $cmd_param[ 0 ] . ": genbank_file output_file\n";
    echo "\nThis PHP script will extract data from the Genbank file and\n";
    echo "enter them into the database.\n";
    echo "\nThe parameters are:\n";
    echo "genbank_file - the GenBank file to be processed.\n";
    echo "               note: there are no error checking on the format.\n";
    echo " output_file - the name of the output file.\n";
    echo "For more information about this script, please consult a magician.\n";
    exit();
}

// open the Genbank and SQL file
$gb_file    = openFile(getenv("PWD") . '/' . $cmd_param[ 1 ], "r");
$out_file   = openFile(getenv("PWD") . '/' . $cmd_param[ 2 ], "w");

$add_data   = false;
$index      = 0;
$total_rec  = 0;
$length     = 0;
$ambiguity  = 0;
$max_length = 0;
$min_length = 2000;

// keep reading the GenBank file until EOF is found
while (!feof($gb_file)) {
    $buffer = fgets($gb_file, 512);       // retrieve one line at a time

    if (strlen($buffer) == 0) {
        continue;
    }

    // if we see a double forward slashes, then we should complete one record
    if (!(strncmp($buffer, '//', 2))) {
        $total_rec++;
        $add_data   = false;
        $locus      = querySafe(getLocus($data));
        $organism   = querySafe(getOrganism($data));
        $accession  = querySafe(getAccession($data));
        $origin     = querySafe(getOrigin($data));

        $l = strlen($origin);

        $max_length = ($max_length < $l) ? $l : $max_length;
        $min_length = ($min_length > $l) ? $l : $min_length;

        $length += $l; $ambiguity += substr_count($origin, "N");

/*
        // attach the forward and reverse primers, if needed
        $forward = "AGAGTTTGATCCTGGCTCAG";
        $reverse = strReverse(strComplement("GYTACCTTGTTACGACTT"));
        $origin = $forward . $origin . $reverse;
*/

        // output the data into the plain text format
        $query_cmd = "$organism|$accession|$locus|$origin\n";
        fwrite($out_file, $query_cmd);
    }

    // it is assumed that LOCUS is the first entry in the file
    if (!(strncmp($buffer, 'LOCUS', 5))) {
        $data = array(); $add_data = true; $index = 0;
    }

    if ($add_data) {
        $data[ $index++ ] = $buffer;
    }
}

echo "   sequence data statistics\n";
echo "---------------------------\n";
echo "    total number of records: $total_rec\n";
echo "maximum length of sequences: $max_length\n";
echo "minimum length of sequences: $min_length\n";
echo "average length of sequences: " . ($length / $total_rec) . "\n";
echo "average number of ambiguity: " . ($ambiguity / $total_rec) . "\n";

fclose($out_file);
fclose($gb_file);     // properly close the file
?>
