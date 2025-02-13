package ru.biosoft.physicell.ui.render;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.SwingConstants;
import javax.swing.plaf.basic.BasicSplitPaneUI.BasicVerticalLayoutManager;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.util.List;

public class RenderTest
{
    private static Vertex cutoff = new Vertex(500, 500, 500);
    
    private static int size = 1000;
    private static int quality = 3;
    public static void main(String[] args)
    {
        JFrame frame = new JFrame();
        Container pane = frame.getContentPane();
        pane.setLayout( new BorderLayout() );

        JSlider sectionXSlider = new JSlider( 0, 1000, 500 );
        JSlider sectionYSlider = new JSlider( 0, 1000, 500 );
        JSlider sectionZSlider = new JSlider( 0, 1000, 500 );
        
        JSlider headingSlider = new JSlider( -180, 180, 0 );   
        JSlider pitchSlider = new JSlider( SwingConstants.VERTICAL, -90, 90, 0 );
        Scene scene = generateScene();
        
        addDisks(scene, cutoff.z, SceneHelper.PLANE_XY);
        addDisks(scene, cutoff.x, SceneHelper.PLANE_YZ);
        addDisks(scene, cutoff.y, SceneHelper.PLANE_XZ);
        
        System.out.println( "Number of cells " + scene.getSpheresCount() );

        RenderPanel renderPanel = new RenderPanel( headingSlider, pitchSlider );
        renderPanel.setCutoff( cutoff );
        renderPanel.setScene( scene );
      
        JPanel leftPanel = new JPanel();
        leftPanel.setLayout(new GridBagLayout() );
        GridBagConstraints con = new GridBagConstraints();
        con.gridy = 0;
        con.gridx = 0;
        con.anchor = GridBagConstraints.NORTH;
        con.ipady = 10;
        leftPanel.add( sectionXSlider , con);
        con.gridy = 1;
        leftPanel.add( sectionYSlider , con);
        con.gridy = 2;
        leftPanel.add( sectionZSlider, con );
        
        sectionZSlider.addChangeListener( e -> {
            cutoff.z = sectionZSlider.getValue();
            renderPanel.setCutoff( cutoff );
            addDisks( scene, cutoff.z, SceneHelper.PLANE_XY );
            renderPanel.repaint();
        } );
        
        sectionXSlider.addChangeListener( e -> {
            cutoff.x = sectionXSlider.getValue();
            renderPanel.setCutoff( cutoff );
            addDisks( scene, cutoff.x, SceneHelper.PLANE_YZ );
            renderPanel.repaint();
        } );
        
        sectionYSlider.addChangeListener( e -> {
            cutoff.y = size-sectionYSlider.getValue();
            renderPanel.setCutoff( cutoff );
            addDisks( scene, cutoff.y, SceneHelper.PLANE_XZ );
            renderPanel.repaint();
        } );
        
        headingSlider.addChangeListener( e -> renderPanel.repaint() );
        pitchSlider.addChangeListener( e -> renderPanel.repaint() );

        pane.add( renderPanel, BorderLayout.CENTER );
        pane.add( pitchSlider, BorderLayout.EAST );
        pane.add( headingSlider, BorderLayout.SOUTH );
        
        pane.add( leftPanel, BorderLayout.LINE_START );
//        pane.add( sectionYSlider, BorderLayout.AFTER_LAST_LINE );
//        pane.add( sectionZSlider, BorderLayout.AFTER_LAST_LINE );

        frame.setSize( size, size );

        headingSlider.setValue( 0 );
        pitchSlider.setValue( 0 );
        frame.setVisible( true );
    }

    public static void addDisks(Scene scene, double d, int plane)
    {
        scene.clearLayer(plane);
        for( Mesh sphere : scene.getSpheres() )
        {
            double distance = getDistance(sphere.center , d, plane);
            if( Math.abs( distance ) < sphere.getRadius() )
            {
                Color meshColor = new Color( 2*sphere.getColor().getRed() , 2*sphere.getColor().getGreen() ,
                        2*sphere.getColor().getBlue() );
                double diskRadius = Math.sqrt( sphere.getRadius() * sphere.getRadius() - distance * distance ) - 1;
                Mesh disk = SceneHelper.createDisk( diskRadius, sphere.center, d , plane, meshColor);
                scene.addDisk( disk, plane );
                //                double innerRadius = circleRadius / 2;
                //                Vertex innerCenter = new Vertex(sphere.center.x, sphere.center.y, sphere.center.z+3);
                //                Mesh circle2 = SceneHelper.createCircle( innerRadius, innerCenter, zCut, outer );
                //                result.add( circle2 );
            }
        }
    }
    
    public static double getDistance(Vertex v, double d, int plane)
    {
        switch( plane )
        {
            case SceneHelper.PLANE_XY:
                return v.z - d;
            case SceneHelper.PLANE_YZ:
                return v.x - d;
            default:
                return v.y - d;
        }
    }
    
    public static Scene generateScene()
    {
        Scene result = new Scene();

        List<Vertex> positions = SceneHelper.createPositions( new Vertex( 500, 500, 500 ), 450, 30);
//                List<Vertex> positions = new ArrayList<>();
//        positions.add( new Vertex( 500, 500, 409 ) );
        for( Vertex v : positions )
        {
            double radius = 30;////7 + 3 * Math.random();
            
            int oncoprotein = (int)(255*Math.random());
            Color outer = new Color( oncoprotein / 2, oncoprotein / 2, ( 255 - oncoprotein ) / 2 );

            Mesh sphere = SceneHelper.createSphere( v.x, v.y, v.z, radius, outer , quality); 
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