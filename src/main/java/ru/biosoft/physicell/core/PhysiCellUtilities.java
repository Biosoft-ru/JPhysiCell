package ru.biosoft.physicell.core;

import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.Random;

import ru.biosoft.physicell.biofvm.GeneralMesh;
import ru.biosoft.physicell.biofvm.Microenvironment;

public class PhysiCellUtilities
{
    private static double seed = 0;
    private static Random r = new Random();
    public static void setSeed(long val)
    {
        seed = val;
        r.setSeed( val );
    }

    public static double getSeed()
    {
        return seed;
    }

    public static double[] UniformOnUnitSphere()
    {
        double[] output = {0, 0, 0};
        double z = 2 * UniformRandom() - 1;// Choose z uniformly distributed in [-1,1].
        double theta = UniformRandom() * 2 * Math.PI;// Choose theta uniformly distributed on [0, 2*pi).
        double r = Math.sqrt( 1 - z * z ); // Let r = sqrt(1-z^2).
        output[0] = r * Math.cos( theta );
        output[1] = r * Math.sin( theta );
        output[2] = z; // (r*cos(theta) , r*sin(theta) , z )
        return output;
    }

    public static double[] UniformOnUnitCircle()
    {
        double[] output = {0, 0, 0};
        double theta = UniformRandom(); //  BioFVM::uniform_random();
        theta *= 2 * Math.PI;//two_pi; // Choose theta uniformly distributed on [0, 2*pi).
        output[0] = Math.cos( theta );
        output[1] = Math.sin( theta ); // (cos(t) , sin(t) , 0 )
        return output;
    }

    public static double UniformRandom(double min, double max)
    {
        return min + ( max - min ) * r.nextDouble();
    }

    public static double NormalRandom()
    {
        return r.nextGaussian();
    }

    public static double NormalRandom(double m, double var)
    {
        return m + var * r.nextGaussian();
    }

    public static double restrict(double val, double min, double max)
    {
        if( val < min )
            return min;
        if( val > max )
            return max;
        return val;
    }

    public static int restrict(int val, int min, int max)
    {
        if( val < min )
            return min;
        if( val > max )
            return max;
        return val;
    }

    public static String getCurrentTime()
    {
        Calendar cal = Calendar.getInstance();
        SimpleDateFormat sdf = new SimpleDateFormat( "[ HH:mm:ss ] " );
        return sdf.format( cal.getTime() );
    }

    public static double UniformRandom()
    {
        return r.nextDouble();
        //        thread_local std::uniform_real_distribution<double> distribution(0.0,1.0);
        //        if( local_pnrg_setup_done == false )
        //        {
        //            // get my thread number 
        //            int i = omp_get_thread_num(); 
        //            physicell_PRNG_generator.seed( physicell_random_seeds[i] ); 
        //            local_pnrg_setup_done = true; 
        //    /*
        //            #pragma omp critical 
        //            {
        //            std::cout << "thread: " << i 
        //            << " seed: " << physicell_random_seeds[i]  << std::endl; 
        //            std::cout << "\t first call: " << distribution(physicell_PRNG_generator) << std::endl; 
        //            }
        //    */
        //        }
        //        return distribution(physicell_PRNG_generator);
        //
        //        // helpful info: https://stackoverflow.com/a/29710970
        //    /*
        //
        //        static std::uniform_real_distribution<double> distribution(0.0,1.0); 
        //        double out;
        //        out = distribution(physicell_PRNG_generator);
        //        return out; 
        //    */  
    }

    public static double[] getCenter(GeneralMesh mesh)
    {
        double[] center = new double[3];
        center[0] = ( mesh.boundingBox[0] + mesh.boundingBox[3] ) / 2;
        center[1] = ( mesh.boundingBox[1] + mesh.boundingBox[4] ) / 2;
        center[2] = ( mesh.boundingBox[2] + mesh.boundingBox[5] ) / 2;
        return center;
    }

    public static double print(double v)
    {
        if( Double.isInfinite( v ) )
            return v;
        return Math.round( v * 2 ) / 2;
    }

    public static void placeInBox(double[] box, CellDefinition cd, int number, Microenvironment m)
    {
        for( int i = 0; i < number; i++ )
        {
            double[] position = {0, 0, 0};
            position[0] = PhysiCellUtilities.UniformRandom( box[0], box[3] );
            position[1] = PhysiCellUtilities.UniformRandom( box[1], box[4] );
            position[2] = PhysiCellUtilities.UniformRandom( box[2], box[5] );
            Cell.createCell( cd, m, position );
        }
    }

    public static void place(Microenvironment m, String type, int number)
    {
        double[] box = m.mesh.boundingBox.clone();
        if( m.options.simulate2D )
        {
            box[2] = 0.0;
            box[5] = 0.0;
        }
        placeInBox( box, type, number, m );
    }

    public static void placeInBox(double[] box, String type, int number, Microenvironment m)
    {
        CellDefinition cd = CellDefinition.getCellDefinition( type );
        for( int i = 0; i < number; i++ )
        {
            double[] position = {0, 0, 0};
            position[0] = PhysiCellUtilities.UniformRandom( box[0], box[3] );
            position[1] = PhysiCellUtilities.UniformRandom( box[1], box[4] );
            position[2] = PhysiCellUtilities.UniformRandom( box[2], box[5] );
            Cell.createCell( cd, m, position );
        }
    }

    public static List<Cell> getNeighbors(Cell pCell)
    {
        List<Cell> neighbors = new ArrayList<>();
        for( Cell neighbor : pCell.get_container().agentGrid.get( pCell.get_current_mechanics_voxel_index() ) )
        {
            neighbors.add( neighbor );
        }
        for( int ind : pCell.get_container().mesh.moore_connected_voxel_indices[pCell.get_current_mechanics_voxel_index()] )
        {
            for( Cell neighbor : pCell.get_container().agentGrid.get( ind ) )
            {
                neighbors.add( neighbor );
            }
        }
        return neighbors;
    }
}
