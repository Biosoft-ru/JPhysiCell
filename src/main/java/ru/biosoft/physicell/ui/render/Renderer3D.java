package ru.biosoft.physicell.ui.render;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.image.BufferedImage;

import ru.biosoft.physicell.ui.DensityState;
import ru.biosoft.physicell.ui.ModelData;
import ru.biosoft.physicell.ui.Visualizer2D.Section;

public class Renderer3D
{
    private static int BLACK_RGB = Color.black.getRGB();
    private boolean axes = true;
    private boolean statistics = true;
    private boolean density = true;
    private String substrate = "oxygen";
    private Point statisticsLocation = new Point( 10, 40 );
    private boolean cut = true;
    private Vertex cutOff = new Vertex( 500, 500, 500 );
    //    private Vertex center;
    private double[] zBuffer;
    private Matrix3 hTransform;
    private Matrix3 pTransform;
    private Matrix3 transform;
    private int xMax;
    private int yMax;
    private int zMax;
    private int xShift;
    private int yShift;
    private int zShift;
    private int xCells;
    private int yCells;
    private int zCells;
    private int dx;
    private int dy;
    private int dz;
    private Vertex shift;
    int width;
    int height;
    int depth;
    private Section sec = Section.Z;
    private int slice = 0;

    private DensityState densityState;

    public Renderer3D(int width, int height, double h, double p)
    {
        //        this( 0, 0, 0, width, height, height, h, p );
        shift = new Vertex( xShift, yShift, zShift );
        zBuffer = new double[width * height];
        //        center = new Vertex( width / 2, width / 2, width / 2 );
        hTransform = new Matrix3( new double[] {Math.cos( h ), 0, -Math.sin( h ), 0, 1, 0, Math.sin( h ), 0, Math.cos( h )} );
        pTransform = new Matrix3( new double[] {1, 0, 0, 0, Math.cos( p ), Math.sin( p ), 0, -Math.sin( p ), Math.cos( p )} );

        transform = hTransform.multiply( pTransform );
    }

    public Renderer3D(ModelData data, double h, double p)
    {
        this.width = (int)data.getXDim().getLength();
        this.height = (int)data.getYDim().getLength();
        this.depth = (int)data.getZDim().getLength();
        this.xShift = -(int)data.getXDim().getFrom();
        this.yShift = -(int)data.getYDim().getFrom();
        this.zShift = -(int)data.getZDim().getFrom();
        this.xMax = (int)data.getXDim().getTo();
        this.yMax = (int)data.getYDim().getTo();
        this.zMax = (int)data.getZDim().getTo();
        dx = (int)data.getXDim().getStep();
        dy = (int)data.getYDim().getStep();
        dz = (int)data.getZDim().getStep();
        xCells = (int)data.getXDim().getLength() / dx;
        yCells = (int)data.getYDim().getLength() / dx;
        zCells = (int)data.getZDim().getLength() / dx;
        shift = new Vertex( xShift, yShift, zShift );
        zBuffer = new double[width * height];
        //        center = new Vertex( width / 2, width / 2, width / 2 );
        hTransform = new Matrix3( new double[] {Math.cos( h ), 0, -Math.sin( h ), 0, 1, 0, Math.sin( h ), 0, Math.cos( h )} );
        pTransform = new Matrix3( new double[] {1, 0, 0, 0, Math.cos( p ), Math.sin( p ), 0, -Math.sin( p ), Math.cos( p )} );
        transform = hTransform.multiply( pTransform );
    }

    public void setDensityState(DensityState densityState)
    {
        this.densityState = densityState;
    }

    public void setSubstrate(String substrate)
    {
        this.substrate = substrate;
    }

    public BufferedImage render(Scene scene, double time)
    {
        BufferedImage img = new BufferedImage( width, height, BufferedImage.TYPE_INT_ARGB );
        
//        if (density && densityState != null)
//        {
//            this.addDensity( densityState.getDensity( substrate ), scene );
//            scene.getLayer( 4 ).forEach( m -> paintDisk( img, m ) );
//        }
        
        for( int q = 0; q < zBuffer.length; q++ )
            zBuffer[q] = Double.NEGATIVE_INFINITY;
        orderMeshes( scene );

        scene.getSpheres().parallelStream().forEach( m -> paintSphere( img, m ) );

        scene.getLayer( SceneHelper.PLANE_XY ).forEach( m -> paintDisk( img, m ) );
        scene.getLayer( SceneHelper.PLANE_YZ ).forEach( m -> paintDisk( img, m ) );
        scene.getLayer( SceneHelper.PLANE_XZ ).forEach( m -> paintDisk( img, m ) );

        if( statistics )
            drawText( scene, time, img.getGraphics() );

        if( axes )
            drawLines( img );
        return img;
    }

    private void drawText(Scene scene, double time, Graphics g)
    {
        g.setFont( new Font( "TimesRoman", Font.PLAIN, 20 ) );
        g.setColor( Color.BLACK );
        g.drawString( "Time: " + time, statisticsLocation.x, statisticsLocation.y );
        g.drawString( "Cells: " + scene.getSpheres().size(), statisticsLocation.x, statisticsLocation.y + 30 );
    }

    private boolean isValid(int x, int y)
    {
        return x >= 0 && x < width && y >= 0 && y < height;
    }


    private void drawLines(BufferedImage img)
    {
        int xLength = Math.max( xMax, 100 );
        int yLength = Math.max( yMax, 100 );
        int zLength = Math.max( zMax, 100 );
        int x = 0;
        int y = 0;
        Vertex v = null;
        for( int i = 0; i < xLength; i++ )
        {
            v = rotate( new Vertex( i, 0, 0 ) );
            x = (int)v.x;
            y = (int)v.y;
            if( isValid( x, y ) )
                img.setRGB( x, y, BLACK_RGB );
        }
        for( int i = 0; i > -yLength; i-- )
        {
            v = rotate( new Vertex( 0, i, 0 ) );
            x = (int)v.x;
            y = (int)v.y;
            if( isValid( x, y ) )
                img.setRGB( x, y, BLACK_RGB );
        }
        for( int i = 0; i < zLength; i++ )
        {
            v = rotate( new Vertex( 0, 0, i ) );
            x = (int)v.x;
            y = (int)v.y;
            if( isValid( x, y ) )
                img.setRGB( x, y, BLACK_RGB );
        }
        paintTriangle( img,
                new Triangle( new Vertex( xLength, 0, 0 ), new Vertex( xLength - 50, -10, 0 ), new Vertex( xLength - 50, 10, 0 ) ),
                Color.black );
        paintTriangle( img,
                new Triangle( new Vertex( xLength, 0, 0 ), new Vertex( xLength - 50, 0, -10 ), new Vertex( xLength - 50, 0, 10 ) ),
                Color.black );
        paintTriangle( img,
                new Triangle( new Vertex( 0, -yLength, 0 ), new Vertex( 0, -yLength + 50, -10 ), new Vertex( 0, -yLength + 50, 10 ) ),
                Color.black );
        paintTriangle( img,
                new Triangle( new Vertex( 0, -yLength, 0 ), new Vertex( -10, -yLength + 50, 0 ), new Vertex( 10, -yLength + 50, 0 ) ),
                Color.black );
        paintTriangle( img,
                new Triangle( new Vertex( 0, 0, zLength ), new Vertex( -10, 0, zLength - 50 ), new Vertex( 10, 0, zLength - 50 ) ),
                Color.black );
        paintTriangle( img,
                new Triangle( new Vertex( 0, 0, zLength ), new Vertex( 0, -10, zLength - 50 ), new Vertex( 0, 10, zLength - 50 ) ),
                Color.black );
    }
    
    public void setDrawDensity(boolean density)
    {
        this.density = density;
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

    //    private Vertex rotate(Vertex v)
    //    {
    //        return transform.transform( ( v.clone().minus( center ) ) ).offset( center );
    //    }

    private Vertex rotate(Vertex v)
    {
        return transform.transform( v.clone() ).offset( shift );
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


    //    Vertex v1;
    //    Vertex v2;
    //    Vertex v3;
    //    
    //    double b1;
    //    double b2;
    //    double b3;
    //    double depth;
    //    int zIndex;
    //    int minX;
    //    int maxX;
    //    int minY;
    //    int maxY;
    //    double triangleArea;
    //    Triangle temp;
    //    Integer shade;

    //    int red;
    //    int green;
    //    int blue;

    private void paintSphere(BufferedImage img, Mesh mesh)
    {
        if( cut && outOfBounds( mesh.center, mesh.radius ) )
            return;

        Vertex meshCenter = rotate( mesh.center );
        Triangle temp;
        for( Triangle t : mesh.triangles )
        {
            if( outOfBounds( t.v1 ) || outOfBounds( t.v2 ) || outOfBounds( t.v3 ) )
                continue;

            temp = rotate( t );

            if( Util.direction( Util.center( temp ), meshCenter ).z > 0 )
                continue;

            reatersizeB( img, temp, mesh.getColor(), true, true );
        }
    }

    private void paintTriangle(BufferedImage img, Triangle t, Color c)
    {
        reatersizeB( img, rotate( t ), c, true, true );
    }

    private void paintDisk(BufferedImage img, Mesh mesh)
    {
        if( cut && outOfBounds( mesh.center ) )
            return;

        for( Triangle t : mesh.triangles )
        {
            if( outOfBounds( Util.center( t ) ) )
                continue;

            reatersizeB( img, rotate( t ), mesh.getColor(), true, true );
        }
    }

    private void reatersizeB(BufferedImage img, Triangle t, Color c, boolean doShade, boolean reversePaint)
    {
        Vertex v1 = t.v1;
        Vertex v2 = t.v2;
        Vertex v3 = t.v3;
        double depth = v1.z + v2.z + v3.z;
        int minX = (int)Math.max( 0, Math.ceil( Util.min( v1.x, v2.x, v3.x ) ) );
        int maxX = (int)Math.min( img.getWidth() - 1, Math.floor( Util.max( v1.x, v2.x, v3.x ) ) );
        int minY = (int)Math.max( 0, Math.ceil( Util.min( v1.y, v2.y, v3.y ) ) );
        int maxY = (int)Math.min( img.getHeight() - 1, Math.floor( Util.max( v1.y, v2.y, v3.y ) ) );
        double triangleArea = ( v1.y - v3.y ) * ( v2.x - v3.x ) + ( v2.y - v3.y ) * ( v3.x - v1.x );
        Integer shade = null;
        for( int y = minY; y <= maxY; y++ )
        {
            for( int x = minX; x <= maxX; x++ )
            {
                int zIndex = y * img.getWidth() + x;
                if( !reversePaint || zBuffer[zIndex] < depth )
                {
                    if( isInside( v1, v2, v3, x, y, triangleArea ) )
                    {
                        if( shade == null && doShade )
                            shade = getShade( c, t );
                        //                    img.setRGB( x/*- 200*/, y/* -200*/, shade );
                        setRGB( img, x, y, shade );

                        if( reversePaint )
                            setBuffer( zIndex, depth );
                        //                    zBuffer[zIndex] = depth;
                    }
                }
            }
        }
    }

    synchronized public void setRGB(BufferedImage img, int x, int y, int shade)
    {
        img.setRGB( x/*- 200*/, y/* -200*/, shade );
    }

    synchronized public void setBuffer(int index, double value)
    {
        zBuffer[index] = value;
    }

    public boolean isInside(Vertex v1, Vertex v2, Vertex v3, double x, double y, double triangleArea)
    {
        double b1 = ( ( y - v3.y ) * ( v2.x - v3.x ) + ( v2.y - v3.y ) * ( v3.x - x ) ) / triangleArea;
        double b2 = ( ( y - v1.y ) * ( v3.x - v1.x ) + ( v3.y - v1.y ) * ( v1.x - x ) ) / triangleArea;
        double b3 = ( ( y - v2.y ) * ( v1.x - v2.x ) + ( v1.y - v2.y ) * ( v2.x - x ) ) / triangleArea;
        return b1 >= 0 && b1 <= 1 && b2 >= 0 && b2 <= 1 && b3 >= 0 && b3 <= 1;
    }

    public boolean outOfBounds(Vertex center)
    {
        return center.z > cutOff.z && center.x > cutOff.x && center.y < cutOff.y;
    }

    public boolean outOfBounds(Vertex center, double radius)
    {
        return center.z - radius > cutOff.z && center.x - radius > cutOff.x && center.y + radius < cutOff.y;
    }

    public Integer getShade(Color color, Triangle t)
    {
        return getShade( color, Math.abs( Util.normal( t.v1, t.v2, t.v3 ).z ) );
    }

    public Integer getShade(Color color, double shade)
    {
        int red = (int)Math.pow( Math.pow( color.getRed(), 2.4 ) * shade, 1 / 2.4 );
        int green = (int)Math.pow( Math.pow( color.getGreen(), 2.4 ) * shade, 1 / 2.4 );
        int blue = (int)Math.pow( Math.pow( color.getBlue(), 2.4 ) * shade, 1 / 2.4 );
        return new Color( red, green, blue ).getRGB();
    }

    public void setAxes(boolean axes)
    {
        this.axes = axes;
    }

    public void setStatistics(boolean statistics)
    {
        this.statistics = statistics;
    }

    public void setStatisticsLOcation(Point location)
    {
        this.statisticsLocation = location;
    }

    public void setAngle(int heading, int pitch)
    {
        double h = Math.toRadians( heading );
        double p = Math.toRadians( pitch );
        hTransform = new Matrix3( new double[] {Math.cos( h ), 0, -Math.sin( h ), 0, 1, 0, Math.sin( h ), 0, Math.cos( h )} );
        pTransform = new Matrix3( new double[] {1, 0, 0, 0, Math.cos( p ), Math.sin( p ), 0, -Math.sin( p ), Math.cos( p )} );
        transform = hTransform.multiply( pTransform );
    }

    private void addDensity(double[] densities, Scene scene)
    {
        int from1 = -this.xShift;
        int from2 = -this.yShift;
        int n1 = xCells;
        int n2 = yCells;
        int size1 = dx;
        int size2 = dy;
        int size3 = dz;
        int shift = this.zShift;
        double maxDensity = 1E-13;
        switch( sec )
        {
            case X:
                from1 = -this.yShift;
                from2 = -this.zShift;
                n1 = yCells;
                n2 = zCells;
                size1 = dy;
                size2 = dz;
                size3 = dx;
                shift = xShift;
                break;
            case Y:
                from1 = -this.xShift;
                from2 = -this.zShift;
                n1 = xCells;
                n2 = zCells;
                size1 = dx;
                size2 = dz;
                size3 = dy;
                shift = yShift;
                break;
            default:
                break;
        }

        int n = ( slice + shift ) / size3;//(int) ( ( options2D.getSlice() + shift ) / size3 );

        double actualMaxDensity = 0;
        for( int i = 0; i < n1; i++ )
        {
            for( int j = 0; j < n2; j++ )
            {
                int red;
                int index;
                switch( sec )
                {
                    case X:
                        index = n + n1 * i + n1 * n2 * j;
                        break;
                    case Y:
                        index = i + n * j + j * n1 * n2;
                        break;
                    default: //Z
                        index = i + n1 * j + n * n1 * n2;
                }
                double density = densities[index];
                if( density > actualMaxDensity )
                    actualMaxDensity = density;

                double ratio = ( density / maxDensity );
                ratio = Math.min( 1, ratio );
                red = (int) ( ( 1 - ratio ) * 255 );

                Mesh square = new Mesh();
                square.setColor( new Color( 255, red, red ) );
                square.add( new Triangle( new Vertex( from1 + i * size1, from2+j*size2, 0), new Vertex(from1+i*size1+size1, from2+j*size2, 0), new Vertex(from1+i*size1+size1, from2+j*size2+size2,0) ) );
                square.add( new Triangle( new Vertex(from1 + i*size1, from2+j*size2, 0), new Vertex(from1+i*size1, from2+j*size2+size2, 0), new Vertex(from1+i*size1+size1, j*size2+size2,0) ) );
                scene.addDisk( square, 4 );
//                scene.
                //                g.setColor( new Color( 255, red, red ) );
                //                g.fillRect( i * size1, j * size2, size1, size2 );
            }
        }
        if( actualMaxDensity > 0 )
        {
            maxDensity = actualMaxDensity;
            if( maxDensity < 1E-20 )
                maxDensity = 1E-20;
        }
    }
}
