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

    public BufferedImage render(Scene scene, double h, double p, int width, int height)
    {
        center = new Vertex( width / 2, width / 2, width / 2 );

        Matrix3 hTransform = new Matrix3( new double[] {Math.cos( h ), 0, -Math.sin( h ), 0, 1, 0, Math.sin( h ), 0, Math.cos( h )} );
        Matrix3 pTransform = new Matrix3( new double[] {1, 0, 0, 0, Math.cos( p ), Math.sin( p ), 0, -Math.sin( p ), Math.cos( p )} );

        Matrix3 transform = hTransform.multiply( pTransform );

        BufferedImage img = new BufferedImage( width, height, BufferedImage.TYPE_INT_ARGB );

        double[] zBuffer = new double[width * height];
        for( int q = 0; q < zBuffer.length; q++ )
            zBuffer[q] = Double.NEGATIVE_INFINITY;

        for( Mesh mesh : scene.getMeshes() )
            paintMesh( mesh, transform, zBuffer, img );

        System.out.println( "P=" +p +" H="+h);
        return img;
    }

    public void setCut(boolean cut)
    {
        this.cut = cut;
    }


    private void paintMesh(Mesh mesh, Matrix3 transform, double[] zBuffer, BufferedImage img)
    {
        if( cut && outOfBounds( mesh.center ) )
            return;

        for( Triangle t : mesh.triangles )
        {
            Vertex v1 = transform.transform( ( t.v1.clone().minus( center ) ) ).offset( center );
            Vertex v2 = transform.transform( ( t.v2.clone().minus( center ) ) ).offset( center );
            Vertex v3 = transform.transform( ( t.v3.clone().minus( center ) ) ).offset( center );

            Vertex norm = Util.normal( v1, v2, v3 );

            double angleCos = Math.abs( norm.z );

            int minX = (int)Math.max( 0, Math.ceil( Util.min( v1.x, v2.x, v3.x ) ) );
            int maxX = (int)Math.min( img.getWidth() - 1, Math.floor( Util.max( v1.x, v2.x, v3.x ) ) );
            int minY = (int)Math.max( 0, Math.ceil( Util.min( v1.y, v2.y, v3.y ) ) );
            int maxY = (int)Math.min( img.getHeight() - 1, Math.floor( Util.max( v1.y, v2.y, v3.y ) ) );

            double triangleArea = ( v1.y - v3.y ) * ( v2.x - v3.x ) + ( v2.y - v3.y ) * ( v3.x - v1.x );

            for( int y = minY; y <= maxY; y++ )
            {
                for( int x = minX; x <= maxX; x++ )
                {
                    double b1 = ( ( y - v3.y ) * ( v2.x - v3.x ) + ( v2.y - v3.y ) * ( v3.x - x ) ) / triangleArea;
                    double b2 = ( ( y - v1.y ) * ( v3.x - v1.x ) + ( v3.y - v1.y ) * ( v1.x - x ) ) / triangleArea;
                    double b3 = ( ( y - v2.y ) * ( v1.x - v2.x ) + ( v1.y - v2.y ) * ( v2.x - x ) ) / triangleArea;
                    if( b1 >= 0 && b1 <= 1 && b2 >= 0 && b2 <= 1 && b3 >= 0 && b3 <= 1 )
                    {
                        double depth = b1 * v1.z + b2 * v2.z + b3 * v3.z;
                        int zIndex = y * img.getWidth() + x;
                        if( zBuffer[zIndex] < depth )
                        {
                            img.setRGB( x, y, getShade( mesh.getColor(), angleCos ).getRGB() );
                            zBuffer[zIndex] = depth;
                        }
                    }
                }
            }
        }
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

    public static Color getShade(Color color, double shade)
    {
        double redLinear = Math.pow( color.getRed(), 2.4 ) * shade;
        double greenLinear = Math.pow( color.getGreen(), 2.4 ) * shade;
        double blueLinear = Math.pow( color.getBlue(), 2.4 ) * shade;

        int red = (int)Math.pow( redLinear, 1 / 2.4 );
        int green = (int)Math.pow( greenLinear, 1 / 2.4 );
        int blue = (int)Math.pow( blueLinear, 1 / 2.4 );

        return new Color( red, green, blue );
    }
}
