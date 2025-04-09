package ru.biosoft.physicell.core;

public class CellIntegrity implements Cloneable
{
    public static double TOLERANCE = 1E-8;
    public double damage = 0;
    public double damage_rate = 0;
    public double damage_repair_rate = 0;

    public void advanceDamage(double dt)
    {
        if( damage_rate > TOLERANCE || damage_repair_rate > TOLERANCE )
        {
            damage = ( damage + dt * damage_rate ) / ( 1 + dt * damage_repair_rate );
        }
    }
    
    @Override
    public CellIntegrity clone()
    {
        try
        {
            return (CellIntegrity)super.clone();
        }
        catch( CloneNotSupportedException e )
        {
            throw ( new InternalError( e ) );
        }
    }
}