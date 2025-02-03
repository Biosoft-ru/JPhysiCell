package ru.biosoft.physicell.ui.render;

import java.awt.Color;
import java.awt.image.BufferedImage;

public class Renderer3D
{

    private boolean cut = true;
    double xCutoff = 350;
    double yCutoff = 450;
    double zCutoff = 350;
    Vertex center;
    double[] zBuffer;
    Matrix3 hTransform;
    Matrix3 pTransform;
    Matrix3 transform;
    BufferedImage img;
    int width;
    int height;

    public Renderer3D(int width, int height, double h, double p)
    {
        zBuffer = new double[width * height];
        center = new Vertex( width / 2, width / 2, width / 2 );
        hTransform = new Matrix3( new double[] {Math.cos( h ), 0, -Math.sin( h ), 0, 1, 0, Math.sin( h ), 0, Math.cos( h )} );
        pTransform = new Matrix3( new double[] {1, 0, 0, 0, Math.cos( p ), Math.sin( p ), 0, -Math.sin( p ), Math.cos( p )} );

        transform = hTransform.multiply( pTransform );
        img = new BufferedImage( width, height, BufferedImage.TYPE_INT_ARGB );
    }

    public BufferedImage render(Scene scene)
    {  
        for( int q = 0; q < zBuffer.length; q++ )
            zBuffer[q] = Double.NEGATIVE_INFINITY;

        for( Mesh mesh : scene.getMeshes() )
            paintMesh( mesh, transform, zBuffer, img );

        //        System.out.println( "P=" + p + " H=" + h );
        return img;
    }

    public void setCut(boolean cut)
    {
        this.cut = cut;
    }

    public static Vertex sub(Vertex v1, Vertex v2)
    {
        Vertex result = v1.clone();
        return result.minus( v2 );
    }

    private static double dotProduct(Vertex v1, Vertex v2)
    {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    private static Vertex crossProduct(Vertex v1, Vertex v2)
    {
        Vertex result = new Vertex( 0, 0, 0 );
        result.x = v1.y * v2.z - v1.z * v2.y;
        result.y = v1.z * v2.x - v1.x * v2.z;
        result.z = v1.x * v2.y - v1.y * v2.x;
        return result;
    }

    public static boolean isBehind(Triangle triangle, Vertex p)
    {
        Vertex a = triangle.v1;
        Vertex b = triangle.v2;
        Vertex c = triangle.v3;

        Vertex ab = sub( b, a );// v2.minus( v1 );
        Vertex ac = sub( c, a );//v3.minus( v1 );

        Vertex normal = crossProduct( ab, ac );

        Vertex ap = sub( p, a );//p.minus( v1 );
        double dotProduct = dotProduct( normal, ap );
        return dotProduct < 0;
    }

    public static void filter(Mesh mesh)
    {
        for( Triangle t : mesh.triangles )
        {
            if( isBehind( t, mesh.center ) )
            {
                mesh.triangles.remove( t );
            }
        }
    }

    int removed = 0;
    public void filter(Scene scene)
    {
        removed = 0;
        for( Mesh mesh : scene.getMeshes() )
        {
            for( Triangle t : mesh.triangles )
            {
                if( isBehind( t, mesh.center ) )
                {
                    mesh.triangles.remove( t );
                    removed++;
                }
            }
        }
        System.out.println( "Removed: " + removed );
        System.out.println( "Remained: " + scene.getMehesCount() * 3 );
    }

    private void paintGrid(Matrix3 transform, BufferedImage img)
    {
        for( int i = 0; i < 500; i++ )
        {
            Vertex vertex = new Vertex( 0, 0, i );
            vertex = transform.transform( vertex );
            //        img.setRGB();
        }
    }

    Vertex triCenter;
    Vertex viewer;
    Vertex toCenter;
    Vertex toViewer;
    Vertex v1;
    Vertex v2;
    Vertex v3;
    Vertex norm;


    double b1;
    double b2;
    double b3;
    double depth;
    int zIndex;

    int minX;
    int maxX;
    int minY;
    int maxY;
    double dotP;
    double angleCos;

    double triangleArea;

    private void paintMesh(Mesh mesh, Matrix3 transform, double[] zBuffer, BufferedImage img)
    {
        if( cut && outOfBounds( mesh.center ) )
            return;

        //        int drawn =0;
        for( Triangle t : mesh.triangles )
        {
            triCenter = new Vertex( ( t.v1.x + t.v2.x + t.v3.x ) / 3, ( t.v1.y + t.v2.y + t.v3.y ) / 3, ( t.v1.z + t.v2.z + t.v3.z ) / 3 );

            viewer = new Vertex( mesh.center.x, mesh.center.y, 10000000 );

            toCenter = mesh.center.clone().minus( triCenter );
            toViewer = mesh.center.clone().minus( viewer );
            dotP = Renderer3D.dotProduct( toCenter, toViewer );
            if( dotP < 0 )
                continue;

            //            drawn++;
            v1 = transform.transform( ( t.v1.clone().minus( center ) ) ).offset( center );
            v2 = transform.transform( ( t.v2.clone().minus( center ) ) ).offset( center );
            v3 = transform.transform( ( t.v3.clone().minus( center ) ) ).offset( center );

            norm = Util.normal( v1, v2, v3 );

            angleCos = Math.abs( norm.z );

            minX = (int)Math.max( 0, Math.ceil( Util.min( v1.x, v2.x, v3.x ) ) );
            maxX = (int)Math.min( img.getWidth() - 1, Math.floor( Util.max( v1.x, v2.x, v3.x ) ) );
            minY = (int)Math.max( 0, Math.ceil( Util.min( v1.y, v2.y, v3.y ) ) );
            maxY = (int)Math.min( img.getHeight() - 1, Math.floor( Util.max( v1.y, v2.y, v3.y ) ) );

            triangleArea = ( v1.y - v3.y ) * ( v2.x - v3.x ) + ( v2.y - v3.y ) * ( v3.x - v1.x );

            for( int y = minY; y <= maxY; y++ )
            {
                for( int x = minX; x <= maxX; x++ )
                {
                    b1 = ( ( y - v3.y ) * ( v2.x - v3.x ) + ( v2.y - v3.y ) * ( v3.x - x ) ) / triangleArea;
                    b2 = ( ( y - v1.y ) * ( v3.x - v1.x ) + ( v3.y - v1.y ) * ( v1.x - x ) ) / triangleArea;
                    b3 = ( ( y - v2.y ) * ( v1.x - v2.x ) + ( v1.y - v2.y ) * ( v2.x - x ) ) / triangleArea;
                    if( b1 >= 0 && b1 <= 1 && b2 >= 0 && b2 <= 1 && b3 >= 0 && b3 <= 1 )
                    {
                        depth = b1 * v1.z + b2 * v2.z + b3 * v3.z;
                        zIndex = y * img.getWidth() + x;
                        if( zBuffer[zIndex] < depth )
                        {
                            img.setRGB( x, y, getShade( mesh.getColor(), angleCos ).getRGB() );
                            zBuffer[zIndex] = depth;
                        }
                    }
                }
            }
        }
        //                System.out.println( "Drawn: " + drawn );
    }

    public boolean outOfBounds(Vertex center)
    {
        return center.z > 750 && center.x > 750 && center.y < 750;
        //        return center.x > 500;
        //        return center.y < 500;
        //        return (center.x < -center.z || center.x > center.z) && center.y > 500;
    }

    public boolean outOfBoundsDefault(Vertex center)
    {
        return center.x > xCutoff && center.y < yCutoff && center.z > zCutoff && center.y > 500;
    }

    int red;
    int green;
    int blue;
    public Color getShade(Color color, double shade)
    {
        red = (int)Math.pow( Math.pow( color.getRed(), 2.4 ) * shade, 1 / 2.4 );
        green = (int)Math.pow( Math.pow( color.getGreen(), 2.4 ) * shade, 1 / 2.4 );
        blue = (int)Math.pow( Math.pow( color.getBlue(), 2.4 ) * shade, 1 / 2.4 );
        return new Color( red, green, blue );
    }
}
