package ru.biosoft.physicell.ui;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageWriter;
import javax.imageio.stream.ImageOutputStream;

import ru.biosoft.physicell.biofvm.BasicAgent;
import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;

public class Visualizer
{
    private AgentVisualizer agentVisualizer = new AgentVisualizer();

    private boolean drawTitle = true;
    private boolean drawAgents = true;
    private boolean drawGrid = false;
    private boolean drawDensity = true;
    private boolean saveImage = true;
    private boolean saveGIF = true;

    private double maxDensity = 1E-13;//6.06;

    private ImageWriter writer = null;
    private ImageOutputStream ios = null;
    private Section sec;
    private double slice;
    private String folder;
    private String name;

    private int xShift = 0;
    private int yShift = 0;
    private int zShift = 0;
    private int substrateIndex = 0;

    public Visualizer()
    {

    }

    public Visualizer setStubstrateIndex(int i)
    {
        this.substrateIndex = i;
        return this;
    }


    public Visualizer(String folder, String name, Section sec, double slice)
    {
        this.folder = folder;
        this.name = name;
        this.sec = sec;
        this.slice = slice;
        setDefaultScheme();
    }

    public void setAgentVisualizer(AgentVisualizer agentVisualizer)
    {
        this.agentVisualizer = agentVisualizer;
    }

    private void setDefaultScheme()
    {
        setColorPhase( "Ki67-", Color.lightGray );
        setColorPhase( "Ki67+ (premitotic)", Color.green );
        setColorPhase( "Ki67+ (postmitotic)", new Color( 0, 128, 0 ) );
        setColorPhase( "Apoptotic", Color.red );
        setColorPhase( "Necrotic (swelling)", Color.magenta );
        setColorPhase( "Necrotic (lysed)", Color.pink );
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
        this.xShift = -(int) ( Math.floor( m.mesh.x_coordinates[0] - m.mesh.dx / 2.0 ) );
        this.yShift = -(int) ( Math.floor( m.mesh.y_coordinates[0] - m.mesh.dy / 2.0 ) );
        this.zShift = -(int) ( Math.floor( m.mesh.z_coordinates[0] - m.mesh.dz / 2.0 ) );

        int xCells = m.mesh.x_coordinates.length;
        int yCells = m.mesh.y_coordinates.length;
        int zCells = m.mesh.y_coordinates.length;
        int width = (int) ( xCells * m.mesh.dx );
        int height = (int) ( yCells * m.mesh.dy );
        switch( sec )
        {
            case X:
                width = (int) ( yCells * m.mesh.dy );
                height = (int) ( zCells * m.mesh.dz );
                break;
            case Y:
                width = (int) ( xCells * m.mesh.dx );
                height = (int) ( zCells * m.mesh.dz );
                break;
            default:
                break;
        }

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

    private void drawCoords(Microenvironment m, Section sec, Graphics g)
    {
        double[] c1 = null;
        double[] c2 = null;

        for( int i = 0; i < c1.length; i++ )
        {

        }
        for( int j = 0; j < c2.length; j++ )
        {

        }
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
            int c1 = x + xShift; //first coordinate
            int c2 = y + yShift; //second coordinate
            double d = Math.abs( z - slice ); //distance from slice;
            switch( sec )
            {
                case X:
                {
                    c1 = y + yShift;
                    c2 = z + zShift;
                    d = Math.abs( x - slice );
                    break;
                }
                case Y:
                {
                    c1 = x + xShift;
                    c2 = z + yShift;
                    d = Math.abs( y - slice );
                    break;
                }
                default:
                    break;
            }

            double radius = agent.getRadius(); //we consider agents to be spheres
            if( d > radius ) //it does not intersect slice;
                continue;
            int r = (int)Math.sqrt( radius * radius - d * d );
            if( agent instanceof Cell )
            {
                //                Cell cell = (Cell)agent;
                agentVisualizer.drawAgent( (Cell)agent, c1, c2, r, g );
                //                String phase = cell.phenotype.cycle.currentPhase().name;
                //                Color color = phaseColor.get( phase );
                //                if( cell.isOutOfDomain )
                //                    continue;
                //
                //                if( color != null )
                //                {
                //                    g.setColor( color );
                //                    g.fillOval( c1 - r, c2 - r, 2 * r, 2 * r );
                //                }
                //                else
                //                {
                //                    System.out.println( "Unknown phase " + phase );
                //                }
                //                g.setColor( Color.black );
                //                g.drawOval( c1 - r, c2 - r, 2 * r, 2 * r );
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

    private void drawDensity(Microenvironment m, double slice, Graphics g)
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
        int shift = this.zShift;
        switch( sec )
        {
            case X:
                n1 = yCells;
                n2 = zCells;
                size1 = sizeY;
                size2 = sizeZ;
                size3 = sizeX;
                shift = xShift;
                break;
            case Y:
                n1 = xCells;
                n2 = zCells;
                size1 = sizeX;
                size2 = sizeZ;
                size3 = sizeY;
                shift = yShift;
                break;
            default:
                break;
        }

        int n = (int) ( ( slice + shift ) / size3 );

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
                double density = m.density[index][substrateIndex];
                if( density > actualMaxDensity )
                    actualMaxDensity = density;

                double ratio = ( density / maxDensity );
                ratio = Math.min( 1, ratio );
                red = (int) ( ( 1 - ratio ) * 255 );

                g.setColor( new Color( 255, red, red ) );
                g.fillRect( i * size1, j * size2, size1, size2 );
            }
        }
        if( actualMaxDensity > 0 )
        {
            //            System.out.println( "Max density: " + actualMaxDensity );
            maxDensity = actualMaxDensity;
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

    public void setColorType(Integer cellType, Color color)
    {
        agentVisualizer.addTypeColor( cellType, color );
    }

    public void setColorPhase(String phase, Color color)
    {
        agentVisualizer.addPhaseColor( phase, color );
    }

    public void setMaxDensity(double density)
    {
        this.maxDensity = density;
    }
}