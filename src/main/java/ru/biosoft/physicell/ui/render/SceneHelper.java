package ru.biosoft.physicell.ui.render;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

public class SceneHelper
{

    public static Mesh createSphere(double[] position, double r, Color color, Color insideColor)
    {
        return createSphere( position[0], position[1], position[2], r, color, insideColor );
    }

    public static Mesh createSphere(double x, double y, double z, double r, Color color, Color insideColor)
    {
        return createSphere( x, y, z, r, color, insideColor, 3);
    }
    
    public static Mesh createSphere(double x, double y, double z, double r, Color color, Color insideColor, int quality)
    {
        Mesh mesh = new Mesh();
        Vertex r111 = new Vertex( x + r, y + r, z + r );
        Vertex r001 = new Vertex( x - r, y - r, z + r );
        Vertex r100 = new Vertex( x + r, y - r, z - r );
        Vertex r010 = new Vertex( x - r, y + r, z - r );
        mesh.add( new Triangle( r111, r001, r010 ) );
        mesh.add( new Triangle( r111, r001, r100 ) );
        mesh.add( new Triangle( r010, r100, r111 ) );
        mesh.add( new Triangle( r010, r100, r001 ) );
        mesh.center = new Vertex( x, y, z );

        for( int i = 0; i < quality; i++ )
             inflate( mesh, r );
        
        for( Vertex v : mesh.getVertices() )
          Util.moveFrom( v, mesh.center, r / Util.distance( v, mesh.center ) );
            
        mesh.setType( Mesh.SPHERE_TYPE );
        mesh.setColor(color);
        mesh.setInsideColor( insideColor );
        mesh.setRadius( r );
        return mesh;
    }
    
    public static final int PLANE_XY = 0;
    public static final int PLANE_YZ = 1;
    public static final int PLANE_XZ = 2;
    
    public static Mesh createDisk(double r, Vertex v, double d, int plane, Color color)
    {
        Vertex center = getCenter( v, d, plane );
        Mesh mesh = new Mesh( center, Mesh.CIRCLE_TYPE, color );
        Vertex prev = getFirst( center, r, d, plane );
        for( int i = 1; i <= 20; i++ )
        {
            Vertex next = getNext( center, r, d, i * Math.PI / 10, plane );
            mesh.add( new Triangle( center.clone(), prev, next ) );
            prev = next.clone();
        }
        return mesh;
    }

    private static Vertex getNext(Vertex v, double r, double d, double phi, int plane)
    {
        switch( plane )
        {
            case PLANE_XY:
                return new Vertex( v.x + r * Math.sin( phi ), v.y + r * Math.cos( phi ), d );
            case PLANE_YZ:
                return new Vertex( d, v.y + r * Math.sin( phi ), v.z + r * Math.cos( phi ) );
            default:
                return new Vertex( v.x + r * Math.cos( phi ), d, v.z + r * Math.sin( phi ) );
        }
    }
    
    private static Vertex getFirst(Vertex v, double r, double d, int plane)
    {
        switch( plane )
        {
            case PLANE_XY:
                return new Vertex( v.x, v.y + r, d );
            case PLANE_YZ:
                return new Vertex( d, v.y, v.z + r );
            default:
                return new Vertex( v.x + r, d, v.z );
        }
    }

    private static Vertex getCenter(Vertex v, double d, int plane)
    {
        switch( plane )
        {
            case PLANE_XY:
                return new Vertex( v.x, v.y, d );
            case PLANE_YZ:
                return new Vertex( d, v.y, v.z );
            default:
                return new Vertex( v.x, d, v.z );
        }
    }

    public static void inflate(Mesh mesh, double radius)
    {
        List<Triangle> newTriangles = new ArrayList<>();
        for( Triangle t : mesh.getTriangles() )
        {
            Vertex c12 = Util.center( t.v1, t.v2 );
            Vertex c23 = Util.center( t.v2, t.v3 );
            Vertex c31 = Util.center( t.v3, t.v1 );

            newTriangles.add( new Triangle( t.v1, c12, c31 ) );
            newTriangles.add( new Triangle( t.v2, c12, c23 ) );
            newTriangles.add( new Triangle( t.v3, c23, c31 ) );
            t.v1 = c12;
            t.v2 = c23;
            t.v3 = c31;
        }
        newTriangles.forEach( t -> mesh.add( t ) );
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
    
    public static void addDisks(Scene scene, double d, int plane)
    {
        scene.clearLayer( plane );
        for( Mesh sphere : scene.getSpheres() )
        {
            double distance = getDistance( sphere.center, d, plane );
            if( Math.abs( distance ) < sphere.getRadius() )
            {
                //                Color meshColor = new Color( 2 * sphere.getColor().getRed(), 2 * sphere.getColor().getGreen(),
                //                        2 * sphere.getColor().getBlue() );
                double diskRadius = Math.sqrt( sphere.getRadius() * sphere.getRadius() - distance * distance ) - 1;
                Mesh disk = SceneHelper.createDisk( diskRadius, sphere.center, d, plane, sphere.getInsideColor() );
                scene.addDisk( disk, plane );
                //                double innerRadius = circleRadius / 2;
                //                Vertex innerCenter = new Vertex(sphere.center.x, sphere.center.y, sphere.center.z+3);
                //                Mesh circle2 = SceneHelper.createCircle( innerRadius, innerCenter, zCut, outer );
                //                result.add( circle2 );
            }
        }
    }

    public static double getDistance(Vertex v, double d, int plane)
    {
        switch( plane )
        {
            case SceneHelper.PLANE_XY:
                return v.z - d;
            case SceneHelper.PLANE_YZ:
                return v.x - d;
            default:
                return v.y - d;
        }
    }
}