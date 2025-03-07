package ru.biosoft.physicell.ui.render;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.image.BufferedImage;

public class Renderer3D
{
    private boolean cut = true;
    private Vertex cutOff = new Vertex(500, 500, 500);
    private Vertex center;
    private double[] zBuffer;
    private Matrix3 hTransform;
    private Matrix3 pTransform;
    private Matrix3 transform;
    
    private BufferedImage img;
    int width;
    int height;

    public Renderer3D(int width, int height, double h, double p)
    {
        this.width = width;
        this.height = height;
        zBuffer = new double[width * height];
        center = new Vertex( width / 2, width / 2, width / 2 );
        hTransform = new Matrix3( new double[] {Math.cos( h ), 0, -Math.sin( h ), 0, 1, 0, Math.sin( h ), 0, Math.cos( h )} );
        pTransform = new Matrix3( new double[] {1, 0, 0, 0, Math.cos( p ), Math.sin( p ), 0, -Math.sin( p ), Math.cos( p )} );

        transform = hTransform.multiply( pTransform );
        img = new BufferedImage( width, height, BufferedImage.TYPE_INT_ARGB );
    }
    
    public BufferedImage render(Scene scene, double time)
    {  
        for( int q = 0; q < zBuffer.length; q++ )
            zBuffer[q] = Double.NEGATIVE_INFINITY;
        orderMeshes( scene );

        for( Mesh mesh : scene.getSpheres() )
            paintSphere( mesh );

        for( Mesh mesh : scene.getLayer( SceneHelper.PLANE_XY ) )
            paintDisk( mesh );

        for( Mesh mesh : scene.getLayer( SceneHelper.PLANE_YZ ) )
            paintDisk( mesh );

        for( Mesh mesh : scene.getLayer( SceneHelper.PLANE_XZ ) )
            paintDisk( mesh );

        drawText(scene, time, img.getGraphics() );
        drawLines(img.getGraphics());
        return img;
    }
    
    private void drawText(Scene scene, double time, Graphics g)
    {
        g.setFont( new Font( "TimesRoman", Font.PLAIN, 20 ) );
        g.setColor( Color.BLACK );
        g.drawString( "Time: " + time, 10, 40 );
        g.drawString( "Cells: " + scene.getSpheres().size() , 10, 70 );
    }

    private static int BLACK_RGB = Color.black.getRGB() ;
    
    private void drawLines(Graphics g)
    {

        for( int i = 750; i < width; i++ )
        {
            Vertex vx = rotate( new Vertex( i, 750, 750 ) );
            int x = (int)vx.x;// - 200;
            int y = (int)vx.y;// - 200;
            if( x >= 0 && x < width && y >= 0 && y <= height )
                img.setRGB( x, y, BLACK_RGB );

            Vertex vy = rotate( new Vertex( 750, 1500-i, 750 ) );
            x = (int)vy.x ;//- 200;
            y = (int)vy.y ;//- 200;
            if( x >= 0 && x < width && y >= 0 && y <= height )
                img.setRGB( x, y, BLACK_RGB );
            
            Vertex vz = rotate( new Vertex( 750, 750, i ) );
            x = (int)vz.x;// - 200;
            y = (int)vz.y;// - 200;
            if( x >= 0 && x < width && y >= 0 && y <= height )
                img.setRGB(x, y, BLACK_RGB );
        }
        paintTriangle(new Triangle(new Vertex(1500, 750, 750), new Vertex(1450, 740, 750), new Vertex(1450, 760, 750)), Color.black);
        paintTriangle(new Triangle(new Vertex(1500, 750, 750), new Vertex(1450, 750, 740), new Vertex(1450, 750, 760)), Color.black);
        paintTriangle(new Triangle(new Vertex(750, 0, 750), new Vertex(750, 50, 740), new Vertex(750, 50, 760)), Color.black);
        paintTriangle(new Triangle(new Vertex(750, 0, 750), new Vertex(740, 50, 750), new Vertex(760, 50, 750)), Color.black);
        paintTriangle(new Triangle(new Vertex(750, 750, 1500), new Vertex(740, 750, 1450), new Vertex(760, 750, 1450)), Color.black);
        paintTriangle(new Triangle(new Vertex(750, 750, 1500), new Vertex(750, 740, 1450), new Vertex(750, 760, 1450)), Color.black);
    }
    
    public void setIsCutOff(boolean cut)
    {
        this.cut = cut;
    }
    
    public void setCutOff(Vertex cut)
    {
        this.cutOff = cut;
    }

    private Triangle rotate(Triangle t)
    {
        return new Triangle( rotate( t.v1 ), rotate( t.v2 ), rotate( t.v3 ) );
    }
    
    private Vertex rotate(Vertex v)
    {
        return transform.transform( ( v.clone().minus( center ) ) ).offset( center );
    }
    
    private void orderMeshes(Scene scene)
    {
        for( Mesh mesh : scene.getSpheres() )
        {
            Vertex center = mesh.center;
            center = rotate( center );
            mesh.setDepth( (int)center.z );
        }
        scene.sortSpheres();
    }

    
    Vertex v1;
    Vertex v2;
    Vertex v3;
    
    double b1;
    double b2;
    double b3;
    double depth;
    int zIndex;
    int minX;
    int maxX;
    int minY;
    int maxY;
    double triangleArea;
    Triangle temp;
    Integer shade;
    
    int red;
    int green;
    int blue;
    
    private void paintSphere(Mesh mesh)
    {
        if( cut && outOfBounds( mesh.center, mesh.radius ) )
            return;

        Vertex meshCenter = rotate( mesh.center );
        for( Triangle t : mesh.triangles )
        {
//            if (outOfBounds( Util.center( t ) ))
//                continue;
            if (outOfBounds( t.v1 ) || outOfBounds( t.v2 ) ||outOfBounds( t.v3 ))
                continue;
            
            temp = rotate( t );

            if( Util.direction( Util.center( temp ), meshCenter ).z > 0 )
                continue;
            //                reatersizeB( temp, mesh.getColor(), false, true );

            reatersizeB( temp, mesh.getColor(), true, true );
        }
    }
    
    private void paintTriangle(Triangle t, Color c)
    {
        reatersizeB( rotate( t ), c, true, true );
    }
    
    private void paintDisk(Mesh mesh)
    {
        if( cut && outOfBounds( mesh.center ) )
            return;

        for( Triangle t : mesh.triangles )
        {
            if (outOfBounds( Util.center( t ) ))
                continue;

            reatersizeB( rotate( t ), mesh.getColor(), true, false );
        }
    }
    
    private void reatersizeB(Triangle t, Color c, boolean doShade, boolean reversePaint)
    {
        v1 = t.v1;
        v2 = t.v2;
        v3 = t.v3;
        depth = v1.z + v2.z + v3.z;    
        minX = (int)Math.max( 0, Math.ceil( Util.min( v1.x, v2.x, v3.x ) ) );
        maxX = (int)Math.min( img.getWidth() - 1, Math.floor( Util.max( v1.x, v2.x, v3.x ) ) );
        minY = (int)Math.max( 0, Math.ceil( Util.min( v1.y, v2.y, v3.y ) ) );
        maxY = (int)Math.min( img.getHeight() - 1, Math.floor( Util.max( v1.y, v2.y, v3.y ) ) );
        triangleArea = ( v1.y - v3.y ) * ( v2.x - v3.x ) + ( v2.y - v3.y ) * ( v3.x - v1.x );        
        shade = null;
        for( int y = minY; y <= maxY; y++ )
        {
            for( int x = minX; x <= maxX; x++ )
            {
                zIndex = y * img.getWidth() + x;
                if( zBuffer[zIndex] < depth && isInside( x, y ) )
                {
                    if( shade == null && doShade)
                        shade = getShade( c, t );
                    img.setRGB( x/*- 200*/, y/* -200*/, shade );
                    zBuffer[zIndex] = depth;
                }
            }
        }
    }
    
    public boolean isInside(double x, double y)
    {
        b1 = ( ( y - v3.y ) * ( v2.x - v3.x ) + ( v2.y - v3.y ) * ( v3.x - x ) ) / triangleArea;
        b2 = ( ( y - v1.y ) * ( v3.x - v1.x ) + ( v3.y - v1.y ) * ( v1.x - x ) ) / triangleArea;
        b3 = ( ( y - v2.y ) * ( v1.x - v2.x ) + ( v1.y - v2.y ) * ( v2.x - x ) ) / triangleArea;
        return b1 >= 0 && b1 <= 1 && b2 >= 0 && b2 <= 1 && b3 >= 0 && b3 <= 1;
    }

    public boolean outOfBounds(Vertex center)
    {
        return center.z > cutOff.z && center.x > cutOff.x && center.y < cutOff.y;
    }
    
    public boolean outOfBounds(Vertex center, double radius)
    {
        return center.z-radius > cutOff.z && center.x-radius > cutOff.x && center.y+radius < cutOff.y;
    }

    public Integer getShade(Color color, Triangle t)
    {
        return getShade( color, Math.abs( Util.normal( t.v1, t.v2, t.v3 ).z ) );
    }
    
    public Integer getShade(Color color, double shade)
    {
        red = (int)Math.pow( Math.pow( color.getRed(), 2.4 ) * shade, 1 / 2.4 );
        green = (int)Math.pow( Math.pow( color.getGreen(), 2.4 ) * shade, 1 / 2.4 );
        blue = (int)Math.pow( Math.pow( color.getBlue(), 2.4 ) * shade, 1 / 2.4 );
        return new Color( red, green, blue ).getRGB();
    }
}
