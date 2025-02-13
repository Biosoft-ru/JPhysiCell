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

    private static double zCut = 480;
    private static int size = 1000;
    private static int quality = 3;
    public static void main(String[] args)
    {
        JFrame frame = new JFrame();
        Container pane = frame.getContentPane();
        pane.setLayout( new BorderLayout() );

        JSlider headingSlider = new JSlider( -180, 180, 0 );
        pane.add( headingSlider, BorderLayout.SOUTH );

        JSlider pitchSlider = new JSlider( SwingConstants.VERTICAL, -90, 90, 0 );
        pane.add( pitchSlider, BorderLayout.EAST );

        Scene scene = generateScene(zCut);
        System.out.println( "Number of cells " + scene.getSpheresCount() );

        RenderPanel renderPanel = new RenderPanel( headingSlider, pitchSlider );
        renderPanel.setZCut( zCut );
        renderPanel.setScene( scene );
        pane.add( renderPanel, BorderLayout.CENTER );

        headingSlider.addChangeListener( e -> renderPanel.repaint() );
        pitchSlider.addChangeListener( e -> renderPanel.repaint() );

        frame.setSize( size, size );

        headingSlider.setValue( 0 );
        pitchSlider.setValue( 0 );
        frame.setVisible( true );
    }

    public static Scene generateScene(double zCut)
    {
        Scene result = new Scene();

        List<Vertex> positions = SceneHelper.createPositions( new Vertex( 500, 500, 500 ), 450, 10);
//                List<Vertex> positions = new ArrayList<>();
//        positions.add( new Vertex( 500, 500, 409 ) );
        for( Vertex v : positions )
        {
            double radius = 10;////7 + 3 * Math.random();
            
            int oncoprotein = (int)(255*Math.random());
            Color outer = new Color( oncoprotein / 2, oncoprotein / 2, ( 255 - oncoprotein ) / 2 );
            Color inner = new Color( oncoprotein , oncoprotein , ( 255 - oncoprotein ) );
            Mesh sphere = SceneHelper.createSphere( v.x, v.y, v.z, radius, outer , quality);
            if( Math.abs( v.z - zCut ) < radius )
            {
                double circleRadius = Math.sqrt( radius * radius - ( v.z - zCut ) * ( v.z - zCut ) );
                Mesh circle = SceneHelper.createCircle( circleRadius, v, zCut, inner );
//                double innerRadius = circleRadius / 2;
//                Vertex innerCenter = new Vertex(sphere.center.x, sphere.center.y, sphere.center.z+3);
//                Mesh circle2 = SceneHelper.createCircle( innerRadius, innerCenter, zCut, outer );
                result.addCircle( circle );
//                result.add( circle2 );
            }
            result.addSphere( sphere );
        }


//        System.out.println( "Cells: "+positions.size() );
//                result.add( Scen/eHelper.createSphere( 500, 500, 500, 80, Color.gray ) );
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

//        result.add(  SceneHelper.createSphere( 500, 500, 500, 100, Color.gray ) );
//        result.add( SceneHelper.createCircle( 100, new Vertex(500,500,500), 4, Color.gray ) );
        //        polygons.add( addSphere( 200, 200, 0, 10 ) );
        //        polygons.add( addSphere( 100, 200, 0, 10 ) );
        //        polygons.add( addSphere( 200, 100, 0, 10 ) );
        //        polygons.add( addSphere( 150, 150, 0, 10 ) );

        return result;
    }
}