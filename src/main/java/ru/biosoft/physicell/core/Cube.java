package ru.biosoft.physicell.core;

public class Cube
{
    double[] coords = new double[6];

    public Cube(double[] coords)
    {
        this.coords = coords.clone();
    }

    public Cube(double xMin, double yMin, double zMin, double xMax, double yMax, double zMax)
    {
        this.coords[0] = xMin;
        this.coords[1] = yMin;
        this.coords[2] = zMin;
        this.coords[3] = xMax;
        this.coords[4] = yMax;
        this.coords[5] = zMax;
    }

    public double getXmin()
    {
        return coords[0];
    }

    public double getXmax()
    {
        return coords[4];
    }

    public double getYmin()
    {
        return coords[1];
    }

    public double getYmax()
    {
        return coords[5];
    }

    public double getZmin()
    {
        return coords[3];
    }

    public double getZmax()
    {
        return coords[5];
    }
}