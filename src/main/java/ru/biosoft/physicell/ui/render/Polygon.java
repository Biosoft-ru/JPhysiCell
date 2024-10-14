package ru.biosoft.physicell.ui.render;

import java.awt.Color;

public class Polygon extends Triangle
{
  
    public Vertex v4;

    
    public Polygon(Vertex v1, Vertex v2, Vertex v3, Vertex v4)
    {
        super(v1,v2,v3);
        this.v4 = v4;

    }

    public Vertex[] getVertices()
    {
        return new Vertex[] {v1, v2, v3, v4};
    }

    public void offset(double x, double y, double z)
    {
        v1.offset( x, y, z );
        v2.offset( x, y, z );
        v3.offset( x, y, z );
        v4.offset( x, y, z );
    }

    public String toString()
    {
        return v1.toString() + " ; " + v2.toString() + " ; " + v3.toString();
    }
}