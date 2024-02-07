package ru.biosoft.physicell.core;

import java.util.Random;

import ru.biosoft.physicell.biofvm.GeneralMesh;

public class PhysiCellUtilities
{
    private static double seed = 0;
    private static Random r = new Random();
    public static void setSeed(long val)
    {
        seed = val;
        r.setSeed( val );
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

}
