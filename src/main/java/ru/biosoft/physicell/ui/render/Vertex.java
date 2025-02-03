package ru.biosoft.physicell.ui.render;

public class Vertex
{
    public double x;
    public double y;
    public double z;

    public Vertex(double x, double y, double z)
    {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public void offset(double x, double y, double z)
    {
        this.x += x;
        this.y += y;
        this.z += z;
    }

    public Vertex offset(Vertex shift)
    {
        this.x += shift.x;
        this.y += shift.y;
        this.z += shift.z;
        return this;
    }
    
    public Vertex minus(Vertex shift)
    {
        this.x -= shift.x;
        this.y -= shift.y;
        this.z -= shift.z;
        return this;
    }
    

    public String toString()
    {
        return "(" + x + " " + y + " " + z + ")";
    }

    public Vertex clone()
    {
        return new Vertex( x, y, z );
    }
}