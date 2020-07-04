<?php

$handle = fopen( "tips.txt", "r" );
$tips = array();

for ( $i = 0; !feof( $handle ); $i++ )
{
    $tips[ $i ] = trim( fgets( $handle, 8192 ) );
    echo "index $i, size: " . strlen( $tips[ $i ] ) . "\n";
}

fclose( $handle );
echo "number of elements: $i\n";
echo "count function: " . count( $tips ) . "\n";

?>
