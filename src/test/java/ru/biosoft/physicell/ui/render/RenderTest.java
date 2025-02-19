package ru.biosoft.physicell.ui.render;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.SwingConstants;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

public class RenderTest
{
    private static Vertex cutoff = new Vertex( 750, 750, 750 );

    private static int size = 1500;
    private static int quality = 2;
    private static boolean playing = false;
    private static File dir;
    private static RenderPanel renderPanel;
    private static Scene scene;
    private static Player player = new Player();
    private static int curTime;
    private static int step;
    private static TreeMap<Integer, File> files = new TreeMap<>();
    private static String folder = "C:/Users/Damag/Documents/Image text";
    private static JSlider timeSlider;

    private static void initDirectory(String path)
    {
        files.clear();
        dir = new File( path );
        for( File f : dir.listFiles() )
        {
            String name = f.getName();
            Integer time = Integer.parseInt( name.split( "_" )[1] );
            files.put( time, f );
        }
        curTime = files.navigableKeySet().first();
        step = files.navigableKeySet().higher( curTime );

        timeSlider.setMaximum( files.navigableKeySet().last() );
        timeSlider.setValue( curTime );
    }

    private static Scene readScene(File f)
    {
        Scene scene = new Scene();
        try (BufferedReader br = new BufferedReader( new FileReader( f ) ))
        {
            String line = br.readLine();
            line = br.readLine();
            while( line != null )
            {
                String[] parts = line.split( "\t" );
                double x = Double.parseDouble( parts[0] );
                double y = Double.parseDouble( parts[1] );
                double z = Double.parseDouble( parts[2] );
                double outerRadius = Double.parseDouble( parts[3] );
                //                double innerRadius = Double.parseDouble( parts[4] );
                Color outerColor = decodeColor( parts[5] );
                Color innerColor = decodeColor( parts[7] );

                Mesh mesh = SceneHelper.createSphere( x, y, z, outerRadius, outerColor, innerColor, quality );
                scene.addSphere( mesh );
                line = br.readLine();
            }
        }
        catch( Exception ex )
        {
            ex.printStackTrace();
        }
        return scene;
    }

    public static Color decodeColor(String s)
    {
        s = s.substring( 1, s.length() - 1 );
        String[] parts = s.split( "," );
        int r = Integer.parseInt( parts[0].trim() );
        int g = Integer.parseInt( parts[1].trim() );
        int b = Integer.parseInt( parts[2].trim() );
        return new Color( r, g, b );
    }


    private static void setScene(Scene newScene)
    {
        scene = newScene;
        SceneHelper.addDisks( scene, cutoff.z, SceneHelper.PLANE_XY );
        SceneHelper.addDisks( scene, cutoff.x, SceneHelper.PLANE_YZ );
        SceneHelper.addDisks( scene, cutoff.y, SceneHelper.PLANE_XZ );
        renderPanel.setScene( scene, curTime);
        renderPanel.repaint();
        System.out.println( "Number of cells " + scene.getSpheresCount() );
    }

    public static void main(String[] args)
    {
        JFrame frame = new JFrame();
        Container pane = frame.getContentPane();
        pane.setLayout( new BorderLayout() );

        JButton browseButton = new JButton( "Load File" );
        JButton playButton = new JButton( "Play" );
        JButton pauseButton = new JButton( "Pause" );
        timeSlider = new JSlider( 0, 100, 0 );
        JSlider sectionXSlider = new JSlider( 0, 1500, 750 );
        JSlider sectionYSlider = new JSlider( 0, 1500, 750 );
        JSlider sectionZSlider = new JSlider( 0, 1500, 750 );
        JSlider headingSlider = new JSlider( -180, 180, 0 );
        JSlider pitchSlider = new JSlider( SwingConstants.VERTICAL, -90, 90, 0 );
        
        renderPanel = new RenderPanel( headingSlider, pitchSlider );
        renderPanel.setSize( 1500, 1500 );
        renderPanel.setPreferredSize(new Dimension( 1500, 1500 ));
        renderPanel.setCutoff( cutoff );

        JPanel leftPanel = new JPanel();
        leftPanel.setLayout( new GridBagLayout() );
        GridBagConstraints con = new GridBagConstraints();
        con.gridy = 0;
        con.gridx = 0;
        con.gridwidth = 2;
        con.anchor = GridBagConstraints.NORTH;
        con.ipady = 10;
        leftPanel.add( browseButton, con );
        con.gridy++;
        leftPanel.add( sectionXSlider, con );
        con.gridy++;
        leftPanel.add( sectionYSlider, con );
        con.gridy++;
        leftPanel.add( sectionZSlider, con );
        con.gridy++;
        con.gridwidth = 1;
        con.weightx = 1;
        leftPanel.add( playButton, con );
        con.gridx++;
        leftPanel.add( pauseButton, con );

        timeSlider.addChangeListener( e -> {

        } );

        sectionZSlider.addChangeListener( e -> {
            cutoff.z = sectionZSlider.getValue();
            renderPanel.setCutoff( cutoff );
            SceneHelper.addDisks( scene, cutoff.z, SceneHelper.PLANE_XY );
            renderPanel.repaint();
        } );

        sectionXSlider.addChangeListener( e -> {
            cutoff.x = sectionXSlider.getValue();
            renderPanel.setCutoff( cutoff );
            SceneHelper.addDisks( scene, cutoff.x, SceneHelper.PLANE_YZ );
            renderPanel.repaint();
        } );

        sectionYSlider.addChangeListener( e -> {
            cutoff.y = size - sectionYSlider.getValue();
            renderPanel.setCutoff( cutoff );
            SceneHelper.addDisks( scene, cutoff.y, SceneHelper.PLANE_XZ );
            renderPanel.repaint();
        } );

        timeSlider.addChangeListener( e -> {
            setTime(timeSlider.getValue());
        });
        
        pauseButton.addActionListener( new ActionListener()
        {
            @Override
            public void actionPerformed(ActionEvent e)
            {
                playing = false;
            }
        } );

        playButton.addActionListener( new ActionListener()
        {
            @Override
            public void actionPerformed(ActionEvent e)
            {
                if (playing)
                    return;
                player = new Player();
                playing = true;
                player.start();
            }
        } );

        browseButton.addActionListener( new ActionListener()
        {
            @Override
            public void actionPerformed(ActionEvent e)
            {
                JFileChooser jchooser = new JFileChooser();
                jchooser.showDialog( frame, "Open" );
                File selectedFile = jchooser.getSelectedFile();
                if( selectedFile != null )
                {
                    setScene( readScene( selectedFile ) );
                }
            }
        } );

        initDirectory(folder);
      
        //        setScene( generateScene() );

        headingSlider.addChangeListener( e -> renderPanel.repaint() );
        pitchSlider.addChangeListener( e -> renderPanel.repaint() );

        JScrollPane scrollPane = new JScrollPane(renderPanel);
        scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        scrollPane.setBounds(50, 30, 500, 500);
        
        pane.add( scrollPane, BorderLayout.CENTER );
        pane.add( pitchSlider, BorderLayout.EAST );
        pane.add( headingSlider, BorderLayout.SOUTH );
        pane.add( timeSlider, BorderLayout.NORTH );
        pane.add( leftPanel, BorderLayout.LINE_START );
        //        pane.add( sectionYSlider, BorderLayout.AFTER_LAST_LINE );
        //        pane.add( sectionZSlider, BorderLayout.AFTER_LAST_LINE );

        frame.setSize( 1500, 1000 );

        headingSlider.setValue( 0 );
        pitchSlider.setValue( 0 );
        frame.setVisible( true );
    }

    public static class Player extends Thread
    {
        @Override
        public void run()
        {
            while( playing )
            {
                if( !files.containsKey( curTime ) )
                {
                    playing = false;
                    return;
                }
                timeSlider.setValue( curTime + step );
            }
        }
    }
    
    private static void setTime(int time)
    {
        curTime = files.navigableKeySet().floor( time );
//        setScene(generateScene());
        setScene( readScene( files.get( curTime  ) ) );
    }

    public static Scene generateScene()
    {
        Scene result = new Scene();

        List<Vertex> positions = SceneHelper.createPositions( new Vertex( 500, 500, 500 ), 250, 10 );
//                        List<Vertex> positions = new ArrayList<>();
//                positions.add( new Vertex( 500, 500, 500 ) );
        for( Vertex v : positions )
        {
            double radius = 10;////7 + 3 * Math.random();

            int oncoprotein = (int) ( 255 * Math.random() );
            Color outer = new Color( oncoprotein / 2, oncoprotein / 2, ( 255 - oncoprotein ) / 2 );
            Color inner = new Color( oncoprotein, oncoprotein, ( 255 - oncoprotein ) );
            Mesh sphere = SceneHelper.createSphere( v.x, v.y, v.z, radius, outer, inner, quality );
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