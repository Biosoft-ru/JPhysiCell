package ru.biosoft.physicell.core;

import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

import ru.biosoft.physicell.biofvm.GeneralMesh;
import ru.biosoft.physicell.biofvm.Microenvironment;

public class PhysiCellUtilities
{
    //    private static long seed = (long) -4.5786084958200968E16;//new Random().nextLong();//3.33256335198003E18
    //    private static Random r = new Random( seed );
    //
    //    public static void setSeed(long val)
    //    {
    //        seed = val;
    //        r.setSeed( val );
    //    }
    //
    //    public static double getSeed()
    //    {
    //        return seed;
    //    }

    public static double[] UniformOnUnitSphere(Model model)
    {
        double[] output = {0, 0, 0};
        double z = 2 * model.getRNG().UniformRandom() - 1;// Choose z uniformly distributed in [-1,1].
        double theta = model.getRNG().UniformRandom() * 2 * Math.PI;// Choose theta uniformly distributed on [0, 2*pi).
        double r = Math.sqrt( 1 - z * z ); // Let r = sqrt(1-z^2).
        output[0] = r * Math.cos( theta );
        output[1] = r * Math.sin( theta );
        output[2] = z; // (r*cos(theta) , r*sin(theta) , z )
        return output;
    }

    public static double[] UniformOnUnitCircle(Model model)
    {
        double[] output = {0, 0, 0};
        double theta = model.getRNG().UniformRandom(); //  BioFVM::uniform_random();
        theta *= 2 * Math.PI;//two_pi; // Choose theta uniformly distributed on [0, 2*pi).
        output[0] = Math.cos( theta );
        output[1] = Math.sin( theta ); // (cos(t) , sin(t) , 0 )
        return output;
    }

    //    public static double UniformRandom(double min, double max)
    //    {
    //        return min + ( max - min ) * r.nextDouble();
    //    }
    //
    //    public static double NormalRandom()
    //    {
    //        return r.nextGaussian();
    //    }
    //
    //    public static double NormalRandom(double m, double var)
    //    {
    //        return m + var * r.nextGaussian();
    //    }
    //
    //    public static double NormalRestricted(double m, double var, double min, double max)
    //    {
    //        return restrict( m + var * r.nextGaussian(), min, max );
    //    }

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

    //    public static boolean checkRandom(double probability)
    //    {
    //        if (probability <= 0)
    //            return false;
    //        else if (probability >= 1)
    //            return true;
    //        return r.nextDouble() < probability;
    //    }
    
    //    public static double UniformRandom()
    //    {
    //        return r.nextDouble();
    //    }

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

    public static void placeInBox(double[] box, CellDefinition cd, int number, Model model) throws Exception
    {
        RandomGenerator rng = model.getRNG();
        for( int i = 0; i < number; i++ )
        {
            double[] position = {0, 0, 0};
            position[0] = rng.UniformRandom( box[0], box[3] );
            position[1] = rng.UniformRandom( box[1], box[4] );
            position[2] = rng.UniformRandom( box[2], box[5] );
            Cell.createCell( cd, model, position );
        }
    }

    public static void place2D(Model model, String type, int number) throws Exception
    {
        double[] box = model.m.mesh.boundingBox.clone();
        box[2] = 0.0;
        box[5] = 0.0;
        placeInBox( box, type, number, model );
    }

    public static void place(Model model, CellDefinition cd, int number) throws Exception
    {
        Microenvironment m = model.getMicroenvironment();
        double[] box = m.mesh.boundingBox.clone();
        if( m.options.simulate2D )
        {
            box[2] = 0.0;
            box[5] = 0.0;
        }
        placeInBox( box, cd, number, model );
    }

    public static void place(Model model, String type, int number) throws Exception
    {
        double[] box = model.m.mesh.boundingBox.clone();
        if( model.m.options.simulate2D )
        {
            box[2] = 0.0;
            box[5] = 0.0;
        }
        placeInBox( box, type, number, model );
    }

    public static void placeInBox(double[] box, String type, int number, Model model) throws Exception
    {
        CellDefinition cd = model.getCellDefinition( type );
        RandomGenerator rng = model.getRNG();
        for( int i = 0; i < number; i++ )
        {
            double[] position = {0, 0, 0};
            position[0] = rng.UniformRandom( box[0], box[3] );
            position[1] = rng.UniformRandom( box[1], box[4] );
            position[2] = rng.UniformRandom( box[2], box[5] );
            Cell.createCell( cd, model, position );
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
