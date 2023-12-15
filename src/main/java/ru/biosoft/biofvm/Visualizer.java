package ru.biosoft.biofvm;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;

import ru.biosoft.biofvm.cell.Cell;

public class Visualizer
{
    private boolean drawTitle = true;
    private boolean drawAgents = true;
    private boolean drawGrid = false;
    private boolean drawDensity = true;
    private double maxDensity = 6.06;

    private Map<String, Color> phaseColor = new HashMap<>();

    public void Visiualizer()
    {

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

    public void setColorPhase(String phase, Color color)
    {
        this.phaseColor.put( phase, color );
    }

    public BufferedImage draw(Microenvironment m, int zCoord, double time, String fileName) throws Exception
    {
        int xCells = m.mesh.x_coordinates.length;
        int yCells = m.mesh.y_coordinates.length;
        int width = (int) ( xCells * m.mesh.dx );
        int height = (int) ( yCells * m.mesh.dy );
        zCoord = (int) ( zCoord / m.mesh.dz );
        double[][] p = m.p_density_vectors;
        BufferedImage img = new BufferedImage( width, height, BufferedImage.TYPE_INT_RGB );
        Graphics g = img.getGraphics();
        g.setColor( Color.white );
        g.fillRect( 0, 0, width, height );
        if (drawGrid)
            drawGrid( xCells, yCells, (int)m.mesh.dx, (int)m.mesh.dy, g );
        if( drawDensity  )
            drawDensity( xCells, yCells, (int)m.mesh.dx, (int)m.mesh.dy, p, zCoord, g );
        if( drawAgents )
            drawAgents( m, g );
        if( drawTitle )
            drawText( m, time, g );
        ImageIO.write( img, "PNG", new File( fileName ) );
        return img;
    }

    private void drawText(Microenvironment m, double time, Graphics g)
    {
        g.setFont( new Font( "TimesRoman", Font.PLAIN, 20 ) );
        g.setColor( Color.BLACK );
        g.drawString( "Time: " + time, 10, 40 );
    }

    private void drawAgents(Microenvironment m, Graphics g)
    {
        List<BasicAgent> agents = BasicAgent.allBasicAgents;
        g.setColor( Color.black );
        for( BasicAgent agent : agents )
        {
            String tag = agent.getTag();
            double[] position = agent.position;
            int x = (int)position[0];
            int y = (int)position[1];
            int rX = (int) ( 2 * agent.getRadius() );
            int rY = (int) ( 2 * agent.getRadius() );
            if( agent instanceof Cell )
            {
                Cell cell = (Cell)agent;
                String phase = cell.phenotype.cycle.current_phase().name;
                Color color = phaseColor.get( phase );
                if( cell.is_out_of_domain )
                {
                    continue;
                }
                if( color != null )
                {
                    g.setColor( color );
                    g.fillOval( x - rX / 2, y - rY / 2, rX, rY );
                }
                g.setColor( Color.black );
                g.drawOval( x - rX / 2, y - rY / 2, rX, rY );
            }
            else
            {
                g.setColor( Color.black );
                if( "Sink".equals( tag ) )
                    g.drawRect( x - rX / 2, y - rY / 2, rX, rY );
                else
                    g.drawOval( x - rX / 2, y - rY / 2, rX, rY );
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
}