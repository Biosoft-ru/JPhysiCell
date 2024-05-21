package ru.biosoft.physicell.ui;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;

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

    private List<ResultGenerator> resultGenerators = new ArrayList<ResultGenerator>();

    private double maxDensity = 1E-13;//6.06;

    private Section sec;
    private double slice;
    private String folder;
    private String name;

    private int xShift = 0;
    private int yShift = 0;
    private int zShift = 0;
    private int substrateIndex = 0;

    public enum Section
    {
        X, Y, Z
    }

    public String getName()
    {
        return name;
    }

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

    public static Visualizer createWithGIF(String folder, String name, Section sec, double slice)
    {
        Visualizer result = new Visualizer( folder, name, sec, slice );
        result.setSaveImage( false );
        result.addResultGenerator( new GIFGenerator( folder, name + ".gif" ) );
        return result;
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

    public void addResultGenerator(ResultGenerator generator)
    {
        resultGenerators.add( generator );
    }

    public void init() throws IOException
    {
        if( saveImage )
        {
            File images = new File( folder + "/" + name );
            images.mkdir();
        }
        for( ResultGenerator generator : this.resultGenerators )
            generator.init();
    }

    public void saveResult(Microenvironment m, double t) throws IOException
    {
        BufferedImage image = draw( m, sec, slice, (int)t );
        if( saveImage )
            ImageIO.write( image, "PNG", new File( folder + "/" + name + "/Figure_" + (int)t + ".png" ) );
        update( image );
    }

    public void update(BufferedImage image) throws IOException
    {
        for( ResultGenerator generator : this.resultGenerators )
            generator.update( image );
    }

    public void finish() throws IOException
    {
        for( ResultGenerator generator : this.resultGenerators )
            generator.finish();
    }

    public BufferedImage getImage(Microenvironment m, double t) throws Exception
    {
        return draw( m, sec, slice, (int)t );
    }


    public BufferedImage draw(Microenvironment m, double slice, double time, String fileName) throws IOException
    {
        return draw( m, Section.Z, slice, time, fileName );
    }

    public BufferedImage draw(Microenvironment m, Section sec, double slice, double time) throws IOException
    {
        return draw( m, sec, slice, time, null );
    }

    public BufferedImage draw(Microenvironment m, Section sec, double slice, double time, String fileName) throws IOException
    {
        int extraWidth = 130;
        int textOffset = 10;
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

        width += extraWidth;
        BufferedImage img = new BufferedImage( width, height, BufferedImage.TYPE_INT_ARGB );
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
            drawText( m, sec, time, width - extraWidth + textOffset, g );
        if( fileName != null )
            ImageIO.write( img, "PNG", new File( fileName ) );
        return img;
    }

    private void drawText(Microenvironment m, Section sec, double time, int x, Graphics g)
    {
        g.setFont( new Font( "TimesRoman", Font.PLAIN, 20 ) );
        g.setColor( Color.BLACK );
        g.drawString( "Time: " + time, x, 40 );
        g.drawString( "Cells: " + m.getAgentsCount(), x, 70 );
        if( sec == Section.X )
            g.drawString( "X = " + slice, x, 100 );
        else if( sec == Section.Y )
            g.drawString( "Y = " + slice, x, 100 );
        else
            g.drawString( "Z = " + slice, x, 100 );
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
                Cell cell = (Cell)agent;
                double nuclearRadius = cell.phenotype.geometry.getNuclearRadius();
                if( d > nuclearRadius )
                    agentVisualizer.drawAgent( cell, c1, c2, r, g );
                else
                {
                    int nr = (int)Math.sqrt( nuclearRadius * nuclearRadius - d * d );
                    agentVisualizer.drawAgent( cell, c1, c2, r, nr, g );
                }
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
            if( maxDensity < 1E-20 )
                maxDensity = 1E-20;
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