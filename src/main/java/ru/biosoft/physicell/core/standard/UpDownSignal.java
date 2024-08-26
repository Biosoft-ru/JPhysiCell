package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.core.Model;

public class UpDownSignal
{
    public double up;
    public double down;
    public boolean noPromoters;
    public boolean noInhibitors;
    public double baseParameter;
    public double maxParameter;
    public double hillPower;
    public double halfMax;

    public UpDownSignal(Model model)
    {
        this( model.getParameterDouble( "hill_power" ), model.getParameterDouble( "half_max" ) );
    }

    public UpDownSignal(double hillPower, double halfMax)
    {
        up = 0.0;
        down = 0.0;

        baseParameter = 0.0;
        maxParameter = 1.0;

        noPromoters = true;
        noInhibitors = true;

        this.hillPower = hillPower;
        this.halfMax = halfMax;
    }

    public void addEffect(double factor, char factor_type)
    {
        // neutral signal 
        if( factor_type == 'N' || factor_type == 'n' )
        {
            return;
        }

        // promoter signal 
        if( factor_type == 'P' || factor_type == 'p' )
        {
            // up = sum of all (scaled) promoter signals 
            up += factor;
            noPromoters = false;
            return;
        }

        // inhibitor signal 
        if( factor_type == 'I' || factor_type == 'i' )
        {
            down += factor;
            noInhibitors = false;
            return;
        }
    }

    public void addEffect(double factor, String factor_type)
    {
        this.addEffect( factor, factor_type.charAt( 0 ) );
    }

    double computeEffectHill()
    {
        double denomConstant = Math.pow( halfMax, hillPower );
        double temp = Math.pow( up, hillPower );
        double up2 = noPromoters ? 0 : temp / ( denomConstant + temp );
        temp = Math.pow( down, hillPower );
        double down2 = noInhibitors ? 1.0 : denomConstant / ( denomConstant + temp );
        return ( baseParameter + ( maxParameter - baseParameter ) * up2 ) * down2; // return UP * DOWN;
    }

    public double computeEffect()
    {
        return computeEffectHill();
    }

    public void reset()
    {
        up = 0.0;
        down = 0.0;
        noPromoters = true;
        noInhibitors = true;
        baseParameter = 0.0;
        maxParameter = 1.0;
    }

    public String toString()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "up    : " + up + " (no promoters : " + noPromoters + ")\n" );
        sb.append( "down  : " + down + " (no inhibiters: " + noInhibitors + ")\n" );
        sb.append( "effect: " + computeEffect() + "\n" );
        return sb.toString();
    }

    public void display()
    {
        System.out.println( toString() );
    }

    public static void main(String ... args)
    {
        System.out.println( " testing ... " );
        UpDownSignal model = new UpDownSignal( 5, 0.1 );
        model.display();

        model.addEffect( 0.2, 'p' );
        model.display();

        model.addEffect( 0.2, 'p' );
        model.display();

        model.addEffect( 0.2, 'n' );
        model.display();

        model.addEffect( 0.2, 'p' );
        model.display();

        model.addEffect( 1, 'p' );
        model.display();

        model.addEffect( 0.1, 'i' );
        model.display();

        model.addEffect( 0.9, 'i' );
        model.display();

        model.addEffect( 0.9, 'i' );
        model.display();

        model.reset();
        model.display();

        model.addEffect( 1, 'i' );
        model.display();

        model.addEffect( 1, 'p' );
        model.display();

        model.addEffect( 10, 'p' );
        model.display();
    }
}