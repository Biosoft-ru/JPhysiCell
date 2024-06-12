package ru.biosoft.physicell.core;

public class Variable implements Cloneable
{
    String name;
    public double value;
    String units;
    boolean conserved_quantity;

    public Variable()
    {
        name = "unnamed";
        units = "dimensionless";
        value = 0.0;
        conserved_quantity = false;
        return;
    }

    public String getName()
    {
        return name;
    }
    public void setName(String name)
    {
        this.name = name;
    }

    public void setValue(double value)
    {
        this.value = value;
    }
    public double getValue()
    {
        return value;
    }

    public boolean isConserved()
    {
        return conserved_quantity;
    }
    public void setConserved(boolean conserved)
    {
        this.conserved_quantity = conserved;
    }

    public String getUnits()
    {
        return units;
    }
    public void setUnits(String units)
    {
        this.units = units;
    }

    @Override
    public String toString()
    {
        return name + ": " + value + " " + units;
    }

    @Override
    public Variable clone()
    {
        try
        {
            return (Variable)super.clone();
        }
        catch( CloneNotSupportedException ex )
        {
            return null;
        }
    }
}