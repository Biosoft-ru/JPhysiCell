package ru.biosoft.physicell.ui.render;

public class Triangle
{
    public Vertex v1;
    public Vertex v2;
    public Vertex v3;

    public Triangle(Vertex v1, Vertex v2, Vertex v3)
    {
        this.v1 = v1;
        this.v2 = v2;
        this.v3 = v3;;
    }

    public Vertex[] getVertices()
    {
        return new Vertex[] {v1, v2, v3};
    }

    public void offset(double x, double y, double z)
    {
        v1.offset( x, y, z );
        v2.offset( x, y, z );
        v3.offset( x, y, z );
    }

    public String toString()
    {
        return v1.toString() + " ; " + v2.toString() + " ; " + v3.toString();
    }
}