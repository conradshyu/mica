#!/usr/bin/php -q

<?php
/*
 * fas2txt.php
 *
 * PHP script that parse FASTA files into the plain text format
 *
 * Copyright (C) 2015   Conrad Shyu (shyu4751@yahoo.com)
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * revised on July 26, 2006
 * revised on May 1, 2015
 * revised on June 17, 2020
 *
*/

define("DELIMITER", ";");

/*
 * note: there is not standard format for the annotation on fasta file
 * it is necessary to manually identify the annotation section and extract
 * required information
*/
class fasta
{
    private $key; private $dna;

    public function __construct($k, $s) {
        $this->dna = trim(str_replace(array("%", "_", "'", "&", ";", "\"", "\n"), "", $s));
        $this->key = trim($k);
    }   // default constructor

    public function __destruct() {
    }   // default destructor

    public function origin() {
        return($this->dna);
    }   // return the sequence

    public function length() {
        return(strlen($this->dna));
    }   // length of the string

    public function accession() {
        $t = explode(DELIMITER, $this->key);
        $a = explode(" ", $t[0]);
        return(trim($a[0]));
    }   // get the accession number

    public function locus() {
        $t = explode(DELIMITER, $this->key);
        return(trim($t[count($t) - 2]));
    }   // get the locus number

    public function strain() {
        $t = explode(DELIMITER, $this->key);
        return(trim($t[count($t) - 1]));
    }   // get the strain name

    public function organism() {
        return($this->strain());
    }   // get the strain name

    public function ambiguity() {
        return(substr_count($this->dna, "N"));
    }   // number of degenerate bases

    public function degenerate() {
        return($this->ambiguity());
    }   // number of degenerate bases

    public function reverse() {
        return(strrev($this->dna));
    }   // find the reverse of a string

    public function antisense() {
        return(strtr($this->dna, "ATUGCYRSWKMBDHVN", "TAACGRYSWMKVHDBN"));
    }   // find the complement of a string

    public function complement() {
        return($this->antisense());
    }   // find the complmenet of a string
}   // smart container implementation

/*
 * this file will load all fasta entries into memory
*/
function load($f) {
    $h = fopen($f, "r"); $list = array(); $fas = ""; $key = "";

    while (!feof($h)) {
        $b = trim(fgets($h, 8192), " \n\t\r\0\x0B");

        if (strlen($b) < 2) {
            continue;
        }   // skip comments and empty line

        if (!($b[0] == '>')) {
            $fas .= strtoupper($b); continue;
        }   // accumulate the sequence

        if (strlen($fas) > 0) {
            $list[$key] = $fas; $fas = "";
        }   // sequence is available

        $key = trim($b, ">"); $list[$key] = 0;
     }   // run through the entire file

    fclose($h); $list[$key] = $fas; return($list);
}   // end of load()

/*
 * main procedure starts here
*/
if ($argc < 3) {
    echo $argv[0] . ": fasta_file output_file\n";
    echo "\nThis PHP script will extract data from the FASTA file and\n";
    echo "convert them into the MiCA database format.\n";
    echo "\nThe parameters are:\n";
    echo " fasta_file - the FASTA file to be processed.\n";
    echo "              note: there are no error checking on the format.\n";
    echo "output_file - the name of the output file.\n";
    echo "For more information about this script, please consult a magician.\n";
    exit();
}   // print out the help

echo "loading file ... ";
$list = load($argv[1]); $db = array();
echo "completed\n";

echo "processing sequence(s) ... ";
foreach (array_keys($list) as $k) {
    $db[] = new fasta($k, $list[$k]);
}   // load the content into memory
echo "completed\n";

$h = fopen($argv[2], "w");

echo "writing database ... ";
foreach ($db as $i) {
    fprintf($h, "%s|%s|%s|%s\n", $i->organism(), $i->accession(), $i->locus(), $i->origin());
}   // write the database file
echo "completed\ntotal number of records: " . count($db) . "\n";

fclose($h);
?>
