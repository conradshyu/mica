<?php
  require "shared.help.inc";    
  DrawHeader( "MiCA: How it Works" );
?>

<p><b>Binary Encoding Scheme</b><br>Unlike traditional character-based approaches, the binary-encoding
scheme utilizes four parallel bit streams to represent the nucleic acid characters. This approach 
exploits the fact that there are only four possible bases (adenine, cytosine, guanine, and thymine). 
Each bit stream reflects the existence of such a nucleotide in the sequence. The bit pattern of '1' 
signals the existence of such a nucleotide at a given position while '0' signals its absence. Since 
nucleotides are individually represented, any possible ambiguities in the sequence can be fully encoded
and resolved. Consequently, instead of comparing each nucleotide explicitly, the entire process is 
reduced into a series of binary operations on each bit stream. Once the bit patterns of integers are 
compared using the binary AND operator, the resulting matrices are added with a binary OR operator. At 
this stage, the names of the nucleotides are no longer relevant because only the number of matches in a 
given subsequence is of interest. In addition, the binary OR operator reduces the resulting matrix into 
a single binary representation of the nucleic matches found at all positions in the subsequence. In 
addition, the binary encoding scheme also makes the calculation of complementary bases trivial. The 
calculation only requires swapping the order of the representing bit streams for the nucleic acids. The 
following figure shows a hypothetical forward primer and short section of template sequences to 
illustrate how the transformation and search are performed.</p>

<p><?php DrawFigure( "Hypothetical template and primer sequences", "images/sequence.png" ); ?></p>

<p><b>Figure 1: </b>The sequence on the top is a hypothetical DNA template and bottom is a forward primer
(5'-KGAGTTTGATCMTGGCTCAG-3'). The subsequence in brackets is the region that the sarch algorithm is 
currently examining. The template and primer sequences will be transformed into binary-coded streams 
before the search proceeds.</p>

<p>The transformation procedures first split up each DNA sequence into four parallel subsequences. Each 
subsequence then represents a particular nucleic acid and encodes the locations where such a nucleotide 
can be found in the sequence. Since each nucleotide is individually represented, any ambiguities can be 
fully resolved and represented. The subsequences are then further simplified to a series of bit streams, 
and each bit pattern encodes the relative locations and occurrences of a nucleotide in the sequence.</p>

<p><?php DrawFigure( "The transformation procedures", "images/split.png" ); ?></p> 

<p><b>Figure 2: </b>The transformation procedures first split up the sequence into four subsequences 
based on the nucleic acides. Each subsequence encodes the locations where such a nucleotide can be found. 
This approach enables ambiguities to be fully resolved and represented. For example, the ambiguitiy code, 
K (at the very beginning of the primer), is a redundant position that encodes both thymine (T) and guanine 
(G), and M (at the 12th position in the primer) is a second redundant position that encodes both adenine 
(A) and cytosine (C). The subsequences are then further simplified to a series of bit streams.</p>

<p>The search begins with a binary AND operation on the bit streams that represent the same nucleotide; 
that is, the bit streams that represent adenine, for example, for the template sequence will be compared 
with the bit streams that represent the same nucleotide for the primer sequence. The binary encoding 
scheme simplifies the search by reducing the complex character-to-character comparisons into a series of 
single step binary operations. The resulting bit streams indicate where a given nucleotide can be found 
on both the template and primer sequences. The following shows the binary AND operations that match 
corresponding bit patterns.</p>

<p><?php DrawFigure( "The binary AND operation on the bit streams", "images/binary_and.png" ); ?></p>

<p><b>Figure 3: </b>The binary AND operator identifies positions where a given nucleotide can be found 
on both the template (top row) and primer (second row) sequences. The binary encoding scheme reduces the 
conventional complex character-to-character comparisons into a series of simple binary operations.</p>

<p>After the bit patterns that represent individual nucleotides are compared, a binary OR operator then 
summarizes all resulting bit streams. At this stage, the names of the nucleotides are no longer relevant
because only the number of matches in a given subsequence is of interest. In addition, the binary OR 
operator reduces the resulting matrix into a single binary representation of the nucleic matches found
at all positions in the subsequence.</p>

<p><?php DrawFigure( "The binary OR operation on the bit streams", "images/binary_or.png" ); ?></p>

<p><b>Figure 4: </b>The binary OR operator summarizes the bit streams that encode the locations where 
given nucleotides can be found on the template and primer sequences. The names of nucleotides are no 
longer relevant because only the number of matches is of interest.</p>

<p>To advance the search along the template sequence, the entire matrix is shifted to the left by one 
bit. This approach is very efficient because all previous bit patterns are still preserved by the 
representation, and it is not necessary to repeat the transformation process for the remainder of 
sequence. An analogous process is used to identify restriction enzyme recognition sites. Once the primer
and restriction enzyme sites have been located, the length of the intervening DNA (in base pairs) can be 
calculated and reported.</p>

<p><?php DrawFigure( "The next nculeotide is being added to the bit streams", 
"images/binary_shift.png" ); ?></p>

<p><b>Figure 5: </b>The next nucleotide is added to the rightmost positions of the bit streams, and the
search advances to the next nucleotide. This example shows that an adenine is being added to the first 
bit streams, with the left-most base (in this example an A) being removed from the stream. All previous
bit patterns are still preserved by the representation.</p>

<p><b>Restriction Sites Search</b><br>
A restriction enzyme cleaves double-stranded DNA at all sites in the molecule where the base sequence 
matches a particular short nucleotide sequence called the restriction site of the enzyme. In naming 
restriction enzymes, the first three letters are italicized because they abbreviate the organism in 
which the enzyme was first discovered, for example, the <em>Eco</em> in <em>Eco</em>RI stands for 
<em>Escherichia coli</em>. These special enzymes recognize specific sequences in DNA molecules, and 
themselves display palindromic sequences. Restriction enzymes are typically of four, five, or six 
nucleotides in length. They function as microbial immune systems that protect bacteria from infection
by viruses. The presence of additional enzymes generally shows a multiplicative effect; a cell with 
multiple independent restriction enzymes could be virtually impregnable. The strict nucleotide sequence
specificity dictates the utilization of restriction enzymes. In evolutionary studies, for example,
restriction sites are commonly used to determine if DNA sequences vary between isolates. The cut made
by restriction enzymes is usually staggered such that the two strands of the double helix are cut a 
few bases apart. This creates single-stranded overhangs called sticky ends at either end of the digested
DNA molecule. The overhangs are said to be sticky because they are able to bind to a complementary 
single-stranded region. The figure shows that the restriction enzyme, <em>Eco</em>RI (5'-GAATTC-3'), 
leaves sticky ends at both ends of a sequence.</p>

<p><?php DrawFigure( "Restriction enzyme EcoRI leaves a sticky end", "images/enzyme_ecori.png" ); ?></p>

<p><b>Figure 6: </b>The restriction enzyme, <em>Eco</em>RI, causes breaks in the backbone of the DNA 
sequence and leaves a staggered cut because it results in fragments with single-stranded ends. The 
single-stranded ends are said to be sticky because they are able to bind to a complementary 
single-stranded region.</p>

<p>Some restriction enzymes are also capable of producing fragments without single-stranded overhangs.
This commonly refers to as a blunt cut. As opposed to sticky cutters, no complementarities are 
required. The following figure shows an example of a blunt cut. The figure illustrates a blunt cut 
made by the restriction enzyme <em>Eco</em>RV (5'-GATATC-3').</p>

<p><?php DrawFigure( "Restriction enzyme EcoRV leaves a blunt end", "images/enzyme_ecorv.png" ); ?></p>

<p><b>Figure 7: </b>The restriction enzyme, <em>Eco</em>RV (5'-GATATC-3') produces fragments without
single-stranded overhangs. This commonly refers to as a blunt cut. As opposed to sticky cutters, no 
complementarities are required.</p>

<?php DrawFooter(); ?>
