package ru.biosoft.physicell.ui.render;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Mesh implements Comparable<Mesh>
{
    public static int SPHERE_TYPE = 0;
    public static int CIRCLE_TYPE = 1;
    
    private Color color;
    double radius;
    private int type = SPHERE_TYPE;
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
    
    public void setType(int type)
    {
        this.type = type;
    }
    
    public int getType()
    {
        return type;
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
    
    public void setColor(Color color)
    {
        this.color = color;
    }
    
    public Color getColor()
    {
        return color;
    }
    
    public void setRadius(double radius)
    {
        this.radius = radius;
    }
    
    public double getRadius()
    {
        return radius;
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