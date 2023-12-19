package ru.biosoft.biofvm;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageWriter;
import javax.imageio.stream.ImageOutputStream;

import ru.biosoft.biofvm.cell.Cell;

public class Visualizer
{
    private boolean drawTitle = true;
    private boolean drawAgents = true;
    private boolean drawGrid = false;
    private boolean drawDensity = true;
    private boolean saveImage = true;
    private boolean saveGIF = true;

    private double maxDensity = 1;//6.06;

    private Map<String, Color> phaseColor = new HashMap<>();

    private ImageWriter writer = null;
    private ImageOutputStream ios = null;
    private Section sec;
    private double slice;
    private String folder;
    private String name;

    public Visualizer()
    {

    }

    public Visualizer(String folder, String name, Section sec, double slice)
    {
        this.folder = folder;
        this.name = name;
        this.sec = sec;
        this.slice = slice;
    }

    public void init() throws IOException
    {
        if( saveImage || saveGIF )
        {
            writer = ImageIO.getImageWritersByFormatName( "GIF" ).next();
            ios = ImageIO.createImageOutputStream( new File( folder + "/" + name + ".gif" ) );
            if( saveImage )
            {
                File images = new File( folder + "/" + name );
                images.mkdir();
            }
            writer.setOutput( ios );
            writer.prepareWriteSequence( null );
        }
    }

    public void saveResult(Microenvironment m, double t) throws Exception
    {
        BufferedImage image = draw( m, sec, slice, (int)t );
        if( saveImage )
            ImageIO.write( image, "PNG", new File( folder + "/" + name + "/Figure_" + (int)t + ".png" ) );
        if( saveGIF )
            writer.writeToSequence( new IIOImage( image, null, null ), writer.getDefaultWriteParam() );
    }

    public void finish() throws IOException
    {
        writer.endWriteSequence();
        writer.reset();
        writer.dispose();
        ios.flush();
        ios.close();
    }



    public enum Section
    {
        X, Y, Z
    }

    public BufferedImage draw(Microenvironment m, double slice, double time, String fileName) throws Exception
    {
        return draw( m, Section.Z, slice, time, fileName );
    }

    public BufferedImage draw(Microenvironment m, Section sec, double slice, double time) throws Exception
    {
        return draw( m, sec, slice, time, null );
    }

    public BufferedImage draw(Microenvironment m, Section sec, double slice, double time, String fileName) throws Exception
    {
        int xCells = m.mesh.x_coordinates.length;
        int yCells = m.mesh.y_coordinates.length;
        int width = (int) ( xCells * m.mesh.dx );
        int height = (int) ( yCells * m.mesh.dy );
        BufferedImage img = new BufferedImage( width, height, BufferedImage.TYPE_INT_RGB );
        Graphics g = img.getGraphics();
        g.setColor( Color.white );
        g.fillRect( 0, 0, width, height );
        if( drawGrid )
            drawGrid( xCells, yCells, (int)m.mesh.dx, (int)m.mesh.dy, g );
        if( drawDensity )
            drawDensity( m, slice, g );
        if( drawAgents )
            drawAgents( m, sec, slice, g );
        if( drawTitle )
            drawText( m, sec, time, g );
        if( fileName != null )
            ImageIO.write( img, "PNG", new File( fileName ) );
        return img;
    }

    private void drawText(Microenvironment m, Section sec, double time, Graphics g)
    {
        g.setFont( new Font( "TimesRoman", Font.PLAIN, 20 ) );
        g.setColor( Color.BLACK );
        g.drawString( "Time: " + time, 10, 40 );
        g.drawString( "Cells: " + m.getAgentsCount(), 10, 70 );
        if( sec == Section.X )
            g.drawString( "X = " + slice, 10, 100 );
        else if( sec == Section.Y )
            g.drawString( "Y = " + slice, 10, 100 );
        else
            g.drawString( "Z = " + slice, 10, 100 );
    }

    private void drawAgents(Microenvironment m, Section sec, double slice, Graphics g)
    {
        g.setColor( Color.black );
        for( BasicAgent agent : m.getAgents() )
        {
            String tag = agent.getTag();
            double[] position = agent.position;
            int x = (int)position[0];
            int y = (int)position[1];
            int z = (int)position[2];
            int c1 = x; //first coordinate
            int c2 = y; //second coordinate
            int c3 = z; //section coordinate;
            switch( sec )
            {
                case X:
                {
                    c1 = y;
                    c2 = z;
                    c3 = x;
                }
                case Y:
                {
                    c1 = x;
                    c2 = z;
                    c3 = y;
                }
                default:
                    break;
            }

            double radius = agent.getRadius(); //we consider agents to be spheres
            double d = Math.abs( c3 - slice );
            if( d > radius ) //it does not intersect slice;
                continue;
            int r = (int)Math.sqrt( radius * radius - d * d );
            if( agent instanceof Cell )
            {
                Cell cell = (Cell)agent;
                String phase = cell.phenotype.cycle.currentPhase().name;
                Color color = phaseColor.get( phase );
                if( cell.isOutOfDomain )
                    continue;

                if( color != null )
                {
                    g.setColor( color );
                    g.fillOval( c1 - r, c2 - r, 2 * r, 2 * r );
                }
                g.setColor( Color.black );
                g.drawOval( c1 - r, c2 - r, 2 * r, 2 * r );
            }
            else
            {
                g.setColor( Color.black );
                if( "Sink".equals( tag ) )
                    g.drawRect( c1 - r, c2 - r, 2 * r, 2 * r );
                else
                    g.drawOval( c1 - r, c2 - r, 2 * r, 2 * r );
            }
        }
    }

    private void drawDensity(Microenvironment m, double Slice, Graphics g)
    {
        int xCells = m.mesh.x_coordinates.length;
        int yCells = m.mesh.y_coordinates.length;
        int zCells = m.mesh.z_coordinates.length;

        int sizeX = (int)m.mesh.dx;
        int sizeY = (int)m.mesh.dy;
        int sizeZ = (int)m.mesh.dz;

        int n1 = xCells;
        int n2 = yCells;
        int size1 = sizeX;
        int size2 = sizeY;
        int size3 = sizeZ;

        switch( sec )
        {
            case X:
                n1 = yCells;
                n2 = zCells;
                size1 = sizeY;
                size2 = sizeZ;
                size3 = sizeX;
                break;
            case Y:
                n1 = xCells;
                n2 = zCells;
                size1 = sizeX;
                size2 = sizeZ;
                size3 = sizeY;
                break;
            default:
                break;
        }

        int n = (int) ( slice / size3 );

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
                red = (int) ( ( maxDensity - m.density[index][0] ) * 255 );
                g.setColor( new Color( 255, red, red ) );
                g.fillRect( i * size1, j * size2, size1, size2 );
            }
        }
    }

    private void drawDensity(int xNumber, int yNumber, int xSize, int ySize, double[][] p, int zCoord, Graphics g)
    {
        for( int i = 0; i < xNumber; i++ )
        {
            for( int j = 0; j < yNumber; j++ )
            {
                int offset = zCoord * xNumber * yNumber;
                int red = (int) ( ( maxDensity - p[offset + i + xNumber * j][0] ) * 255 );
                g.setColor( new Color( 255, red, red ) );
                g.fillRect( i * xSize, j * ySize, xSize, ySize );
            }
        }
    }

    private void drawGrid(int xNumber, int yNumber, int xSize, int ySize, Graphics g)
    {
        for( int i = 0; i < xNumber; i++ )
        {
            for( int j = 0; j < yNumber; j++ )
            {
                g.setColor( Color.black );
                g.drawRect( i * xSize, j * ySize, xSize, ySize );
            }
        }
    }

    public void setDrawTitle(boolean drawTitle)
    {
        this.drawGrid = drawTitle;
    }

    public void setDrawGrid(boolean drawGrid)
    {
        this.drawGrid = drawGrid;
    }

    public void setDrawAgents(boolean drawAgents)
    {
        this.drawGrid = drawAgents;
    }

    public void setDrawDensity(boolean drawDensity)
    {
        this.drawDensity = drawDensity;
    }

    public void setSaveImage(boolean saveImage)
    {
        this.saveImage = saveImage;
    }

    public void setSaveGIF(boolean saveGIF)
    {
        this.saveGIF = saveGIF;
    }

    public void setColorPhase(String phase, Color color)
    {
        this.phaseColor.put( phase, color );
    }
}