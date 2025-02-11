package ru.biosoft.physicell.ui.render;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Mesh implements Comparable<Mesh>
{
    private Color color;
    List<Triangle> triangles = new ArrayList<>();
    public Vertex center;
    public int depth;

    public Mesh(Vertex center)
    {
        this.center = center;
    }
    public Mesh()
    {
        this.center = new Vertex( 0, 0, 0 );
    }

    public Mesh(List<Triangle> triangles, Vertex center)
    {
        this.triangles = triangles;
        this.center = center;
    }

    public List<Triangle> getTriangles()
    {
        return triangles;
    }

    public Set<Vertex> getVertices()
    {
        Set<Vertex> result = new HashSet<>();
        for( Triangle triangle : triangles )
            for( Vertex vertex : triangle.getVertices() )
                result.add( vertex );
        return result;
    }

    public void add(Triangle triangle)
    {
        triangles.add( triangle );
    }

    public void offset(double x, double y, double z)
    {
        for( Triangle triangle : triangles )
            triangle.offset( x, y, z );
        center.offset( x, y, z );
    }

    public void section(Mesh polygon)
    {
        for( Triangle triangle : polygon.triangles )
        {
            section( triangle );
        }
    }

    public void section(Triangle triangle)
    {
    }
    
    public void setColor(Color color)
    {
        this.color = color;
    }
    
    public Color getColor()
    {
        return color;
    }
    
    @Override
    public int compareTo(Mesh o)
    {
        return -Integer.compare( depth , (int)((Mesh)o).depth);
    }
    
    public void setDepth(int depth)
    {
        this.depth = depth;
    }
}