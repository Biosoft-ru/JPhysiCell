package ru.biosoft.physicell.ui.render;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;

import javax.swing.JPanel;
import javax.swing.JSlider;

public class RenderPanel extends JPanel
{
    private final JSlider headingSlider;
    private final JSlider pitchSlider;
    private Scene scene;
    public double[] x;

    Renderer3D renderer;// = new Renderer3D();



    BufferedImage img;

    public RenderPanel(JSlider headingSlider, JSlider pitchSlider)
    {
        this.headingSlider = headingSlider;
        this.pitchSlider = pitchSlider;
    }

    public void setScene(Scene scene)
    {
        this.scene = scene;
    }

    public BufferedImage getImage()
    {
        return img;
    }

    public void paintComponent(Graphics g)
    {
        double t0 = System.currentTimeMillis();
        Graphics2D g2 = (Graphics2D)g;
        g2.setColor( Color.BLACK );
        g2.fillRect( 0, 0, getWidth(), getHeight() );

        double heading = Math.toRadians( headingSlider.getValue() );
        double pitch = Math.toRadians( pitchSlider.getValue() );
        
        renderer= new Renderer3D( getWidth(), getHeight() , heading, pitch );
        renderer.setCut( true );
        img = renderer.render( scene);

        g2.drawImage( img, 0, 0, null );
        System.out.println( "ELAPSED: " + ( System.currentTimeMillis() - t0 ) / 1000 );
    }


}