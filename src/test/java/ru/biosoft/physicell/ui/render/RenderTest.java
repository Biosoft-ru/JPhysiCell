package ru.biosoft.physicell.ui.render;

import javax.swing.JFrame;
import javax.swing.JSlider;
import javax.swing.SwingConstants;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.util.List;


public class RenderTest
{

    public static void main(String[] args)
    {
        JFrame frame = new JFrame();
        Container pane = frame.getContentPane();
        pane.setLayout( new BorderLayout() );

        JSlider headingSlider = new JSlider( -180, 180, 0 );
        pane.add( headingSlider, BorderLayout.SOUTH );

        JSlider pitchSlider = new JSlider( SwingConstants.VERTICAL, -90, 90, 0 );
        pane.add( pitchSlider, BorderLayout.EAST );

        Scene scene = generateScene();
        System.out.println( "Number of cells " + scene.getMehesCount() );

        RenderPanel renderPanel = new RenderPanel( headingSlider, pitchSlider );
        renderPanel.setScene( scene );
        pane.add( renderPanel, BorderLayout.CENTER );

        headingSlider.addChangeListener( e -> renderPanel.repaint() );
        pitchSlider.addChangeListener( e -> renderPanel.repaint() );

        frame.setSize( 1000, 1000 );


//        headingSlider.setValue( -44 );
//        pitchSlider.setValue( -22 );
        headingSlider.setValue( 0 );
      pitchSlider.setValue( 0 );
        frame.setVisible( true );
    }

    public static Scene generateScene()
    {
        Scene result = new Scene();

//        List<Vertex> positions = SceneHelper.createPositions( new Vertex( 500, 500, 500 ), 200, 10 );
//        for( Vertex v : positions )
//        {
//            result.add( SceneHelper.createSphere( v.x, v.y, v.z, 5 + 6 * Math.random(), Color.gray ) );
//        }
//                result.add( SceneHelper.createSphere( 300, 400, -400, 20, Color.gray ) );
        //        result.add( SceneHelper.createSphere( 350, 400, -300, 10 ) );
        //        result.add( SceneHelper.createSphere( 400, 400, -200, 10 ) );
        //        result.add( SceneHelper.createSphere( 450, 400, -100, 10 ) );
//                result.add( SceneHelper.createSphere( 500, 400, 0, 10, Color.gray  ) );
        //        result.add( SceneHelper.createSphere( 550, 400, 100, 10 ) );
        //        result.add( SceneHelper.createSphere( 600, 400, 200, 10 ) );
        //        result.add( SceneHelper.createSphere( 650, 400, 300, 10 ) );
        //        result.add( SceneHelper.createSphere( 700, 400, 400, 10 ) );



        //        List<Vertex> positions = createPositions( 250, 18 );
        //        for( Vertex v : positions )
        //        {
        //            v.offset( 400, 400, 400 );
        //            meshes.add( createSphere( v.x, v.y, v.z, 5+10*Math.random() ) );
        //        }
        //        for( int i = 0; i < 100; i++ )
        //            meshes.add( createSphere( 1000 * Math.random(), 1000 * Math.random(), 1000 * Math.random(), 10 + Math.random() * 10 ) );

        result.add(  SceneHelper.createSphere( 300, 300, 0, 100, Color.gray ) );
        //        polygons.add( addSphere( 200, 200, 0, 10 ) );
        //        polygons.add( addSphere( 100, 200, 0, 10 ) );
        //        polygons.add( addSphere( 200, 100, 0, 10 ) );
        //        polygons.add( addSphere( 150, 150, 0, 10 ) );

        return result;
    }
}