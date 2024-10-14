package ru.biosoft.physicell.ui.render;

public class Util
{
    public static double length(Vertex v)
    {
        return Math.sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
    }

    public static Vertex normalize(Vertex v)
    {
        double normalLength = length( v );
        v.x /= normalLength;
        v.y /= normalLength;
        v.z /= normalLength;
        return v;
    }

    public static void moveFrom(Vertex v, Vertex from, double distance)
    {
        v.x = from.x + ( v.x - from.x ) * distance;
        v.y = from.y + ( v.y - from.y ) * distance;
        v.z = from.z + ( v.z - from.z ) * distance;
    }

    public static double distance(Vertex v1, Vertex v2)
    {
        return Math.sqrt( ( v1.x - v2.x ) * ( v1.x - v2.x ) + ( v1.y - v2.y ) * ( v1.y - v2.y ) + ( v1.z - v2.z ) * ( v1.z - v2.z ) );
    }

    public static Vertex center(Vertex v1, Vertex v2)
    {
        return new Vertex( ( v1.x + v2.x ) / 2, ( v1.y + v2.y ) / 2, ( v1.z + v2.z ) / 2 );
    }
    
    public static double min(double x1, double x2, double x3)
    {
        return Math.min(Math.min(x1, x2), x3);
    }
    
    public static double max(double x1, double x2, double x3)
    {
        return Math.max(Math.max(x1, x2), x3);
    }
    
    public static Vertex multiply(Vertex v1, Vertex v2)
    {
        return new Vertex( v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x );
    }
    
    public static Vertex direction(Vertex v1, Vertex v2)
    {
        return new Vertex( v2.x - v1.x, v2.y - v1.y, v2.z - v1.z );
    }
    
    public static Vertex normal(Vertex v1, Vertex v2, Vertex v3)
    {
        Vertex ab = Util.direction( v1, v2 );
        Vertex ac = Util.direction( v1, v3 );
        Vertex norm = Util.multiply( ab, ac );
        return Util.normalize( norm );
    }
}