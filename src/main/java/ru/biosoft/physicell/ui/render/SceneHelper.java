package ru.biosoft.physicell.ui.render;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

public class SceneHelper
{

    public static Mesh createSphere(double[] position, double r, Color color)
    {
        return createSphere( position[0], position[1], position[2], r, color );
    }

    public static Mesh createSphere(double x, double y, double z, double r, Color color)
    {
        return createSphere( x, y, z, r, color, 4);
    }
    
    public static Mesh createSphere(double x, double y, double z, double r, Color color, int quality)
    {
        Mesh mesh = new Mesh();
        mesh.setColor( color );
        mesh.add( new Triangle( new Vertex( r, r, r ), new Vertex( -r, -r, r ), new Vertex( -r, r, -r ) ) );
        mesh.add( new Triangle( new Vertex( r, r, r ), new Vertex( -r, -r, r ), new Vertex( r, -r, -r ) ) );
        mesh.add( new Triangle( new Vertex( -r, r, -r ), new Vertex( r, -r, -r ), new Vertex( r, r, r ) ) );
        mesh.add( new Triangle( new Vertex( -r, r, -r ), new Vertex( r, -r, -r ), new Vertex( -r, -r, r ) ) );

        mesh.offset( x, y, z );

        for( int i = 0; i < quality; i++ )
            mesh = inflate( mesh, r );
        return mesh;
    }
    
    public static Mesh createCircle(double r, Vertex v, double d,  Color color)
    {
        Vertex center = new Vertex( v.x, v.y, d );
        Mesh mesh = new Mesh();
        mesh.setColor( color );
        mesh.setType( Mesh.CIRCLE_TYPE );
        mesh.center = center;
        Vertex prev = new Vertex( center.x, center.y + r, d );
        for( int i = 1; i <= 20; i++ )
        {
            double phi = i * Math.PI / 10;
            Vertex next = new Vertex( center.x + r * Math.sin( phi ), center.y + r * Math.cos( phi ), d );
            Triangle t = new Triangle( center.clone(), prev, next );
            mesh.add( t );
            prev = next.clone();
        }
        return mesh;
    }

    public static Mesh inflate(Mesh mesh, double radius)
    {
        Mesh result = new Mesh( mesh.center );
        result.setType( Mesh.SPHERE_TYPE );
        result.setColor( mesh.getColor() );
        result.setRadius( radius );
        for( Triangle t : mesh.getTriangles() )
        {
            Vertex c12 = Util.center( t.v1, t.v2 );
            Vertex c23 = Util.center( t.v2, t.v3 );
            Vertex c31 = Util.center( t.v3, t.v1 );

            result.add( new Triangle( t.v1, c12, c31 ) );
            result.add( new Triangle( t.v2, c12, c23 ) );
            result.add( new Triangle( t.v3, c23, c31 ) );
            result.add( new Triangle( c12, c23, c31 ) );
        }

        for( Vertex v : result.getVertices() )
            Util.moveFrom( v, mesh.center, radius / Util.distance( v, mesh.center ) );

        return result;
    }

    public static List<Vertex> createPositions(Vertex center, double sphereRadius, double cellRadius)
    {
        List<Vertex> result = new ArrayList<>();
        int xc = 0, zc = 0;
        double xSpacing = cellRadius * Math.sqrt( 3 );
        double ySpacing = cellRadius * 2;
        double zSpacing = cellRadius * Math.sqrt( 3 );

        for( double z = -sphereRadius; z < sphereRadius * 2; z += zSpacing, zc++ )
        {
            for( double x = -sphereRadius; x < sphereRadius * 2; x += xSpacing, xc++ )
            {
                for( double y = -sphereRadius; y < sphereRadius * 2; y += ySpacing )
                {
                    Vertex tempPoint = new Vertex( x + ( zc % 2 ) * 0.5 * cellRadius, y + ( xc % 2 ) * cellRadius, z );
                    tempPoint.offset( center );
                    if( Util.distance( tempPoint, center ) < sphereRadius )
                    {
                        result.add( tempPoint );
                    }
                }
            }
        }
        return result;
    }
}