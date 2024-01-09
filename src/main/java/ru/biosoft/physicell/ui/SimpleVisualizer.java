package ru.biosoft.physicell.ui;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JPanel;

import ru.biosoft.physicell.biofvm.BasicAgent;
import ru.biosoft.physicell.biofvm.Microenvironment;

public class SimpleVisualizer
{

    public static void main(String ... ars)
    {
        //        show( null );
        try
        {
            draw( 1000, 1000, "C:/Users/Damag/BIOFVM/test.png" );
        }
        catch( Exception ex )
        {

        }
    }


    public static void readFile(String path)
    {
        try (BufferedReader br = new BufferedReader( new FileReader( new File( path ) ) ))
        {
            double[] p;
            String data = br.readLine();

            while( !"HEADER".equals( data ) )
                data = br.readLine();
            String[] vals = data.split( "\t" );
            int voxels = Integer.valueOf( vals[1] );
            p = new double[voxels];
            while( !"DATA".equals( data ) )
                data = br.readLine();

            int x;
            int y;
            int z;
            while( data != null && !data.isEmpty() )
            {
                vals = data.split( "\t" );
                x = Integer.valueOf( vals[0] );
                y = Integer.valueOf( vals[1] );
                z = Integer.valueOf( vals[2] );
                x = ( x - 5 ) / 10;
                y = ( y - 5 ) / 10;
                z = ( z - 5 ) / 10;
                //                p = Integer.valueOf( vals[4] );
            }
        }
        catch( Exception ex )
        {

        }
    }

    public static void show(Microenvironment m)
    {
        JFrame frame = new JFrame();
        MyPanel panel = new MyPanel();
        frame.setContentPane( panel );
        frame.pack();
        frame.setSize( 1000, 1000 );


        frame.setVisible( true );
        while( true )
            ;
    }

    public static class MyPanel extends JPanel
    {
        @Override
        public void paintComponent(Graphics g)
        {
            super.paintComponent( g );
            Graphics2D g2d = (Graphics2D)g;
            int xLength = 100;
            int yLength = 100;
            double[] p = new double[xLength * yLength];
            for( int i = 0; i < p.length; i++ )
                p[i] = Math.random();

            int frameWidth = 1000;
            int frameHeight = 1000;

            int actualWidth = frameWidth / xLength;
            int actualHeight = frameHeight / yLength;

            JPanel panel = new JPanel();
            JFrame frame = new JFrame();
            frame.add( panel );

            panel.setBounds( new Rectangle( 0, 0, frameWidth, frameHeight ) );
            panel.setVisible( true );
            drawDensity( xLength, yLength, actualWidth, actualHeight, p, 0, g );
        }
    }

    public static void draw(int width, int height, String fileName) throws Exception
    {
        BufferedImage img = new BufferedImage( width, height, BufferedImage.TYPE_INT_RGB );
        int actualWidth = width / 1000;
        int actualHeight = height / 1000;
        double[] p = new double[1000 * 1000];
        for( int i = 0; i < p.length; i++ )
            p[i] = Math.random();
        drawDensity( 1000, 1000, actualWidth, actualHeight, p, 0, img.getGraphics() );
        ImageIO.write( img, "PNG", new File( fileName ) );
    }

    public static BufferedImage draw(Microenvironment m, int zCoord, double time, String fileName) throws Exception
    {
        int xCells = m.mesh.x_coordinates.length;
        int yCells = m.mesh.y_coordinates.length;
        int width = (int) ( xCells * m.mesh.dx );
        int height = (int) ( yCells * m.mesh.dy );
        zCoord = (int) ( zCoord / m.mesh.dz );
        double[][] p = m.density;
        BufferedImage img = new BufferedImage( width, height, BufferedImage.TYPE_INT_RGB );
        Graphics g = img.getGraphics();
        drawDensity( xCells, yCells, (int)m.mesh.dx, (int)m.mesh.dy, p, zCoord, g );
        drawAgents( m, g );
        drawText( m, time, g );
        ImageIO.write( img, "PNG", new File( fileName ) );
		return img;
    }

    private static BufferedImage resizeImage(BufferedImage originalImage, int targetWidth, int targetHeight) throws IOException
    {
        Image resultingImage = originalImage.getScaledInstance( targetWidth, targetHeight, Image.SCALE_REPLICATE );
        BufferedImage resizedImage = new BufferedImage( targetWidth, targetHeight, BufferedImage.TYPE_INT_RGB );
        Graphics2D graphics2D = resizedImage.createGraphics();
        graphics2D.drawImage( resultingImage, 0, 0, targetWidth, targetHeight, null );
        graphics2D.dispose();
        return resizedImage;
    }

    private static void drawText(Microenvironment m, double time, Graphics g)
    {
        g.setFont( new Font( "TimesRoman", Font.PLAIN, 20 ) );
        g.setColor( Color.BLACK );
        g.drawString( "Time: " + time, 10, 40 );
    }

    private static void drawAgents(Microenvironment m, Graphics g)
    {
        g.setColor( Color.black );
        for( BasicAgent agent : m.getAgents() )
        {
            String tag = agent.getTag();
            double[] position = agent.position;
            int x = (int)position[0];
            int y = (int)position[1];
            int rX = (int) ( 2 * agent.getRadius() );
            int rY = (int) ( 2 * agent.getRadius() );
            if( "Source".equals( tag ) )
                g.drawOval( x - rX / 2, y - rY / 2, rX, rY );
            else
                g.drawRect( x - rX / 2, y - rY / 2, rX, rY );
        }
    }

    private static void drawDensity(int xNumber, int yNumber, int xSize, int ySize, double[][] p, int zCoord, Graphics g)
    {
        for( int i = 0; i < xNumber; i++ )
        {
            for( int j = 0; j < yNumber; j++ )
            {
                int offset = zCoord * xNumber * yNumber;
                int red = (int) ( ( 1 - p[offset + i + xNumber * j][0] ) * 255 );
                g.setColor( new Color( 255, red, red ) );
                g.fillRect( i * xSize, j * ySize, xSize, ySize );
                //                g.setColor( Color.black );
                //                g.drawRect( i * xSize, j * ySize, xSize, ySize );
            }
        }
    }

    private static void drawDensity(int xNumber, int yNumber, int xSize, int ySize, double[] p, int zCoord, Graphics g)
    {
        for( int i = 0; i < xNumber; i++ )
        {
            for( int j = 0; j < yNumber; j++ )
            {
                int offset = zCoord * xNumber * yNumber;
                int red = (int) ( p[offset + zCoord + i + xNumber * j] * 255 );
                g.setColor( new Color( 255, red, red ) );
                g.fillRect( i * xSize, j * ySize, xSize, ySize );
            }
        }
    }
}
