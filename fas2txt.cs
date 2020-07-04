/*
 * fas2txt.cs
 * This program converts sequence file in the fasta format into mica database.
 * Fasta format is generally not consistent, therefore it is necessary to update
 * the code in order to accommendate difference in the annotation section.
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
*/

/*
 * note: there is no standard format for the annotation on fasta file
 * it is necessary to manually identify the annotation section and extract
 * required information
*/
class fsa2txt
{
	/*
	 * this implementation loads the entire file into memory
	*/
    static System.Int32 fasta(
        System.Collections.Generic.Dictionary<System.String, System.String> _l,
        System.String _f )
    {
        System.String k = "";
        System.Text.StringBuilder s = new System.Text.StringBuilder();

        try
        {
            using ( System.IO.StreamReader r = new System.IO.StreamReader( _f ) )
            {
                System.String b = "";

                while ( ( b = r.ReadLine() ) != null )
                {
                    if ( b.Length < 2  )
                    {
                        continue;
                    }   // skip comments and empty line

                    if ( !( b[ 0 ] == '>' ) )
                    {
                        s.Append( b ); continue;
                    }   // accumulate the sequence

                    if ( ( s.ToString() ).Length > 0 )
                    {
                        _l[ k ] = ( s.ToString() ).Trim(); s.Length = 0;
                    }   // sequence is available

                    k = b.Trim( '>' ); _l.Add( k, "" );
                }   // set the taxonomic identification and rank
            }	// file is automatically closed
        }
        catch ( System.Exception e )
        {
            System.Console.WriteLine( "Exception: {0}", e.Message );
        }   // an exception has occurred

        _l[ k ] = ( s.ToString() ).Trim(); return( _l.Count );
    }   // extract sequences from fasta file

    public static int Main( System.String[] args )
    {
        if ( args.Length < 2 )
        {
            System.Console.WriteLine( "require: fasta output" ); return( 0 );
        }   // check the required parameters

        System.Collections.Generic.Dictionary<System.String, System.String> list =
            new System.Collections.Generic.Dictionary<System.String, System.String>();

        System.Console.Write( "processing: {0} ... ", args[ 0 ] );
        fasta( list, args[ 0 ] ); System.Console.WriteLine( "completed" );

        try
        {
            using ( System.IO.StreamWriter w = new System.IO.StreamWriter( args[ 1 ], false ) )
            {
                foreach( System.Collections.Generic.KeyValuePair<System.String, System.String> a in list )
                {
                    System.String[] id = ( a.Key ).Split( '|' );
                    w.WriteLine( "{0}|{1}|{2}|{3}",
					    ( ( ( id[ 4 ] ).Split( ',' ) )[ 0 ] ).Trim(),
                        ( id[ 1 ] ).Trim(), ( id[ 3 ] ).Trim(), ( a.Value ).Trim() );
                }   // process each entry
            }   // export the sequence into a file
        }
        catch ( System.Exception e )
        {
            System.Console.WriteLine( "Exception: {0}", e.Message );
        }   // an exception has occurred

        return( 1 ); 
    }
}
