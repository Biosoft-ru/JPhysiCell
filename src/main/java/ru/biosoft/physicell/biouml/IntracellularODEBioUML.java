package ru.biosoft.physicell.biouml;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.CopyOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;

import biouml.model.Diagram;
import biouml.plugins.simulation.InfiniteSpan;
import biouml.plugins.simulation.Model;
import biouml.plugins.simulation.Simulator;
import biouml.plugins.simulation.java.EventLoopSimulator;
import biouml.plugins.simulation.java.JavaSimulationEngine;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Intracellular;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.ode.IntracellularODE;
import ru.biosoft.util.TempFiles;

public class IntracellularODEBioUML extends IntracellularODE
{
    double next_librr_run = 0;
    protected Map<String, Integer> variableIndex = null;
    private Diagram diagram;
    private Model model;
    private Simulator solver = new EventLoopSimulator();
    private double dt;
    private String sbmlPath;
    private String[] inputs;
    private String[] outputs;
    private JavaSimulationEngine engine = new JavaSimulationEngine();

    public IntracellularODEBioUML(ru.biosoft.physicell.core.Model model, CellDefinition cd) throws Exception
    {
        super( model, cd );
    }

    public static String copySRC() throws IOException
    {
        InputStream is = IntracellularODEBioUML.class.getResourceAsStream( "src.jar" );
        Path path = TempFiles.dir( "JPhysicell_temp" ).toPath().resolve( "src.jar" );
        Files.copy( is, path, new CopyOption[0] );
        return path.toString();
    }

    public void setDiagram(Diagram diagram) throws Exception
    {
        this.diagram = diagram;
        engine.setDiagram( diagram );
        String path = copySRC();
        engine.setClassPath( path );
        engine.setLogLevel( Level.OFF );
        model = engine.createModel();
        variableIndex = engine.getShortNameMapping();
        model.init();
    }
    public Diagram getDiagram()
    {
        return diagram;
    }

    public String getSBMLPath()
    {
        return sbmlPath;
    }

    public void setModel(Model model)
    {
        this.model = model;
    }

    public void setSolver(Simulator simulator)
    {
        this.solver = simulator;
    }

    @Override
    public boolean need_update(double curTime)
    {
        return curTime >= this.next_librr_run;
    }

    @Override
    public void step() throws Exception
    {
        if( solver.getProfile().getTime() == 0 )
            solver.init( model, model.getInitialValues(), new InfiniteSpan( dt ), null, null );
        else
            solver.setInitialValues( model.getCurrentValues() );
        solver.doStep();
    }

    public void setVarIndexes(Map<String, Integer> mapping)
    {
        this.variableIndex = new HashMap<>( mapping );
    }

    public String getVariableName(String name)
    {
        return phenotypeSpecies.get( name );
    }

    @Override
    public double getParameterValue(String name) throws Exception
    {
        int index = variableIndex.get( name );
        return model.getCurrentValues()[index];
    }

    @Override
    public void setParameterValue(String name, double value) throws Exception
    {
        int index = variableIndex.get( name );
        double[] vals = model.getCurrentValues();
        vals[index] = value;
        model.setCurrentValues( vals );
    }

    @Override
    public void update(Cell cell, Phenotype phenotype, double dt)
    {
        // TODO Auto-generated method stub

    }

    @Override
    public void inherit(Cell cell)
    {
        // TODO Auto-generated method stub

    }

    @Override
    public void setDT(double dt)
    {
        this.dt = dt;
    }

    @Override
    public void addPhenotypeSpecies(String code, String species)
    {
        this.phenotypeSpecies.put( code, species );
    }

    @Override
    public Intracellular clone()
    {
        try
        {
            IntracellularODEBioUML result = (IntracellularODEBioUML)super.clone();
            result.phenotypeSpecies = new HashMap<>( phenotypeSpecies );
            result.diagram = diagram;
            if( model != null )
                result.model = model.clone();
            if( solver != null )
                result.solver = solver.getClass().newInstance();
            result.dt = dt;
            result.solver.init( result.model, result.model.getInitialValues(), new InfiniteSpan( dt ), null, null );
            return result;

        }
        catch( Exception e )
        {
            throw ( new InternalError( e ) );
        }
    }

    @Override
    public String[] getInputs()
    {
        return inputs;
    }

    @Override
    public String[] getOutputs()
    {
        return this.phenotypeSpecies.keySet().toArray( String[]::new );
    }

    @Override
    public void setInputs(String[] inputs)
    {
        this.inputs = inputs;
    }

    @Override
    public void setOutputs(String[] outputs)
    {
        this.outputs = outputs;
    }
}