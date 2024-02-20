package ru.biosoft.physicell.core;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;

public class CellCSVReader
{

    static void load_cells_csv_v1(BufferedReader br, Microenvironment m) throws Exception
    {
        //        try (BufferedReader br = new BufferedReader( new InputStreamReader( inputStream ) ))
        //        {
        String s = br.readLine();
        while( s != null )
        {
            String[] data = s.split( "," );
            if( data.length != 4 )
                System.out.println( "Error! Importing cells from a CSV file expects each row to be x,y,z,typeID." );

            double[] position = {Double.parseDouble( data[0] ), Double.parseDouble( data[1] ), Double.parseDouble( data[2] )};
            int my_type = (int)Math.round( Double.parseDouble( data[3] ) );
            CellDefinition pCD = CellDefinition.getCellDefinition( my_type );
            if( pCD != null )
            {
                //                System.out.println( "Creating " + pCD.name + " (type=" + pCD.type + ") at " + VectorUtil.print( position ) );
                Cell pCell = Cell.createCell( pCD, m, position );
            }
            else
            {
                System.out.println( "Warning! No cell definition found for index " + my_type + "!" + "\tIgnoring cell at position "
                        + VectorUtil.print( position ) );
            }
            s = br.readLine();
        }
        //        }
    }

    static Cell process_csv_v2_line(String line, String[] labels, Microenvironment m)
    {
        String[] tokens = line.split( "," );
        double[] position = new double[3];
        position[0] = Double.parseDouble( tokens[0] );
        position[1] = Double.parseDouble( tokens[1] );
        position[2] = Double.parseDouble( tokens[2] );

        // the cell type 
        String celltype = tokens[3];
        CellDefinition pCD = CellDefinition.getCellDefinition( celltype );
        if( pCD == null )
        {
            System.out.println( "Warning! CSV file requests creating cell type " + celltype + "\tat " + position
                    + "but I don't recognize that type. Skipping cell!" );
            return null;
        }

        // create the cell IF the definition was found 
//        System.out.println( "Creating " + pCD.name + " (type=" + pCD.type + ") at " + VectorUtil.print( position ) );

        Cell pCell = Cell.createCell( pCD, m, position );
        return pCell;
        // now write any extra data 
        //        for( int k=4 ; k < tokens.length; k++ )
        //        {
        //            double dval = Double.parseDouble(tokens[k] ); 
        //            boolean processed = false; 
        //            boolean skip = false; 


        //            // if the string is empty, skip 
        //            if( tokens[k].isEmpty() )
        //            { skip = true; }
        //            else
        //            {
        //                char c = tokens[k].c_str()[0]; 
        //                // use 's' or 'S' to skip the entry 
        //                if( c == 's' || c == 'S' )
        //                { skip = true; }
        //            }

        // special cases: 

        //                // volume 
        //            if( labels[k] == "volume" && skip == false )
        //            { 
        //                pCell.setTotalVolume( dval ); 
        //                processed = true; 
        //            }
        //
        //            // check behavior dictionary 
        //
        //            if( processed == false && skip == false )
        //            {
        //                // if the behavior is found in the dictionary, process it 
        //                if( find_behavior_index( labels[k] ) > -1 )
        //                {
        //                    set_single_behavior( pCell , labels[k] , dval ); 
        //                    processed = true; 
        //                }
        //            }
        //
        //            // warning message for any unprocessed variables 
        //            if( processed == false && skip == false )
        //            {
        //                std::cout << "\tWarning: I don't know how to process " << labels[k] 
        //                << " so I skipped it." << std::endl;
        //            }
        //            // give a notation for any intentinoally skipped variables 
        //            if( skip == true )
        //            {
        //                std::cout << "\tNote: Skipping " << labels[k] 
        //                << " for this cell." << std::endl;
        //            }


        //        }
        //
        //        return pCell;  
    }

    static void load_cells_csv_v2(BufferedReader br, Microenvironment m) throws Exception
    {
        //        File f = new File( filename );
        //        if( !f.exists() )
        //            throw new Exception( "Error: " + filename + " not found during cell loading. Quitting." );
        //        System.out.println( "Loading cells from simple (v2) CSV file " + filename + " ... " );
        //        try (BufferedReader br = new BufferedReader( new InputStreamReader( inputStream ) ))
        //        {
        String s = br.readLine();
        String[] labels = s.split( "," );
        s = br.readLine();
        while( s != null )
        {
            process_csv_v2_line( s, labels, m );
            s = br.readLine();
        }
        //        }
    }

    public static void load_cells_csv(InputStream inputStream, Microenvironment m) throws Exception
    {
        try (BufferedReader br = new BufferedReader( new InputStreamReader( inputStream ) ))
        {
            String s = br.readLine();
            char c = s.charAt( 0 );
            if( c == 'X' || c == 'x' )
            {
                // v2
                load_cells_csv_v2( br, m );
            }
            else
            {
                // v1
                load_cells_csv_v1( br, m );
            }
        }
    }
}