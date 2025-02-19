package ru.biosoft.physicell.ui.render;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import javax.imageio.ImageIO;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.ui.AgentColorer;
import ru.biosoft.physicell.ui.AgentColorerBW;
import ru.biosoft.physicell.ui.ResultGenerator;
import ru.biosoft.physicell.ui.Visualizer;

public class Visualizer3D implements Visualizer
{
    private boolean saveImage = true;

    private List<ResultGenerator> resultGenerators = new ArrayList<ResultGenerator>();
    private AgentColorer colorer = new AgentColorerBW();
    private Microenvironment m;
    private String folder;
    private String name;

    double xCutoff = 350;
    double yCutoff = 450;
    double zCutoff = 350;
    
    double p = -0.3839;
    double h = -0.7679;
    
    public String getName()
    {
        return name;
    }

    public Visualizer3D(String folder, String name, Microenvironment m)
    {
        this.folder = folder;
        this.name = name;
        this.m = m;
        render = new Renderer3D( (int)m.mesh.boundingBox[3], (int)m.mesh.boundingBox[4], h, p);
        scene  = new Scene();
    }

    public void addResultGenerator(ResultGenerator generator)
    {
        resultGenerators.add( generator );
    }
    public void clearResultGenerators()
    {
        resultGenerators.clear();
    }

    public Iterable<ResultGenerator> getGenerators()
    {
        return resultGenerators;
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
        BufferedImage image = draw( (int)t );
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
    
    Renderer3D render;
    Scene scene;
    
    private BufferedImage draw(double time) throws IOException
    {
        scene.clear();
        
        Set<Cell> cells = m.getAgents( Cell.class );
        for( Cell cell : cells )
        {
            Color c = colorer.findColors( cell )[1];
            if( c == null )
                c = Color.gray;
            Mesh sphere = SceneHelper.createSphere( cell.position, cell.getRadius(), c, null );
//            Renderer3D.filter( sphere );
            scene.addSphere( sphere );
        }

        return render.render( scene, time );
    }

    @Override
    public BufferedImage getImage(Microenvironment m, double t) throws Exception
    {
        return draw( (int)t );
    }

    public Visualizer setAgentColorer(AgentColorer colorer)
    {
        this.colorer = colorer;
        return this;
    }
}
