package ru.biosoft.physicell.core;

import java.util.Random;

public class RandomGenerator
{
    private long seed = (long) -4.5786084958200968E16;//new Random().nextLong();//3.33256335198003E18
    private Random r = new Random( seed );

    public void setSeed(long val)
    {
        seed = val;
        r.setSeed( val );
    }

    public long getSeed()
    {
        return seed;
    }

    //    public double[] UniformOnUnitSphere()
    //    {
    //        double[] output = {0, 0, 0};
    //        double z = 2 * UniformRandom() - 1;// Choose z uniformly distributed in [-1,1].
    //        double theta = UniformRandom() * 2 * Math.PI;// Choose theta uniformly distributed on [0, 2*pi).
    //        double r = Math.sqrt( 1 - z * z ); // Let r = sqrt(1-z^2).
    //        output[0] = r * Math.cos( theta );
    //        output[1] = r * Math.sin( theta );
    //        output[2] = z; // (r*cos(theta) , r*sin(theta) , z )
    //        return output;
    //    }
    //
    //    public double[] UniformOnUnitCircle()
    //    {
    //        double[] output = {0, 0, 0};
    //        double theta = UniformRandom(); //  BioFVM::uniform_random();
    //        theta *= 2 * Math.PI;//two_pi; // Choose theta uniformly distributed on [0, 2*pi).
    //        output[0] = Math.cos( theta );
    //        output[1] = Math.sin( theta ); // (cos(t) , sin(t) , 0 )
    //        return output;
    //    }

    public double UniformRandom(double min, double max)
    {
        return min + ( max - min ) * UniformRandom();
    }
    
    public double LogUniformRandom(double min, double max)
    {
        min = Math.log( min);
        max = Math.log( max );
        return Math.exp(min + ( max - min ) * UniformRandom());
    }

    public double NormalRandom()
    {
        return r.nextGaussian();
    }

    public double NormalRandom(double m, double var)
    {
        return m + var * NormalRandom();
    }
    
    public double NormalRandom(double mu, double sigma, double min, double max)
    {
        double value = NormalRandom( mu, sigma );
        while( value <= min || value >= max )
            value = NormalRandom( mu, sigma );
        return value;
    }
    
    public double LogNormalRandom(double mu, double sigma, double min, double max)
    {
        double value = Math.exp( NormalRandom( mu, sigma) );
        while (value <= min || value >= max)
            value =  Math.exp( NormalRandom( mu, sigma) );
        return value;
    }
    
    public double Log10NormalRandom(double mu, double sigma, double min, double max)
    {
        double value = Log10NormalRandom( mu, sigma );
        while( value <= min || value >= max )
            value = Log10NormalRandom( mu, sigma );
        return value;
    }
    
    public double Log10NormalRandom(double mu, double sigma)
    {
        return Math.exp(Math.log( 10 ) * NormalRandom(mu, sigma));
    }

    public double NormalRestricted(double m, double var, double min, double max)
    {
        return restrict( m + var * NormalRandom(), min, max );
    }

    public double restrict(double val, double min, double max)
    {
        if( val < min )
            return min;
        if( val > max )
            return max;
        return val;
    }

    public boolean checkRandom(double probability)
    {
        if( probability <= 0 )
            return false;
        else if( probability >= 1 )
            return true;
        return r.nextDouble() < probability;
    }

    public double UniformRandom()
    {
        return r.nextDouble();
    }
    

//  Uniform","LogUniform","Normal","LogNormal","Log10Normal
    
    
    //    public void placeInBox(double[] box, CellDefinition cd, int number, Model model)
    //    {
    //        for( int i = 0; i < number; i++ )
    //        {
    //            double[] position = {0, 0, 0};
    //            position[0] = PhysiCellUtilities.UniformRandom( box[0], box[3] );
    //            position[1] = PhysiCellUtilities.UniformRandom( box[1], box[4] );
    //            position[2] = PhysiCellUtilities.UniformRandom( box[2], box[5] );
    //            Cell.createCell( cd, model, position );
    //        }
    //    }
    //
    //    public void place2D(Model model, String type, int number)
    //    {
    //        double[] box = model.m.mesh.boundingBox.clone();
    //        box[2] = 0.0;
    //        box[5] = 0.0;
    //        placeInBox( box, type, number, model );
    //    }
    //
    //    public void place(Model model, CellDefinition cd, int number)
    //    {
    //        Microenvironment m = model.getMicroenvironment();
    //        double[] box = m.mesh.boundingBox.clone();
    //        if( m.options.simulate2D )
    //        {
    //            box[2] = 0.0;
    //            box[5] = 0.0;
    //        }
    //        placeInBox( box, cd, number, model );
    //    }
    //
    //    public void place(Model model, String type, int number)
    //    {
    //        double[] box = model.m.mesh.boundingBox.clone();
    //        if( model.m.options.simulate2D )
    //        {
    //            box[2] = 0.0;
    //            box[5] = 0.0;
    //        }
    //        placeInBox( box, type, number, model );
    //    }
    //
    //    public void placeInBox(double[] box, String type, int number, Model model)
    //    {
    //        CellDefinition cd = model.getCellDefinition( type );
    //        for( int i = 0; i < number; i++ )
    //        {
    //            double[] position = {0, 0, 0};
    //            position[0] = PhysiCellUtilities.UniformRandom( box[0], box[3] );
    //            position[1] = PhysiCellUtilities.UniformRandom( box[1], box[4] );
    //            position[2] = PhysiCellUtilities.UniformRandom( box[2], box[5] );
    //            Cell.createCell( cd, model, position );
    //        }
    //    }
}
