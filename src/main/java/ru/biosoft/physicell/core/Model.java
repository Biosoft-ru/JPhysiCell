package ru.biosoft.physicell.core;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.ui.Visualizer;
import ru.biosoft.physicell.ui.Visualizer.Section;

public class Model
{
    private List<Visualizer> visualizers = new ArrayList<Visualizer>();
    private String logFile;

    private Microenvironment m;
    private Map<String, String> parameters = new HashMap<>();
    private double tMax;
    private double diffusion_dt;
    private double mechanics_dt = 0.1;
    private double phenotype_dt = 6.0;
    private double full_save_interval;
    private String resultFolder;
    private boolean enableFullSaves;
    private boolean hasEvents = false;
    private List<Event> events = new ArrayList<>();

    public Iterable<Visualizer> getVisualizers()
    {
        return visualizers;
    }

    public void setResultFolder(String folder)
    {
        this.resultFolder = folder;
    }

    public void addEvent(Event event)
    {
        this.events.add( event );
    }

    public void setLogFile(String path)
    {
        this.logFile = path;
    }

    public Visualizer addVisualizer(int zSlice, String name)
    {
        Visualizer visualizer = new Visualizer( resultFolder, name, Section.Z, zSlice );
        //        visualizer.setDrawDensity( false );
        visualizer.setSaveImage( false );
        //        visualizer.setColorPhase( "Ki67-", Color.lightGray );
        //        visualizer.setColorPhase( "Ki67+ (premitotic)", Color.green );
        //        visualizer.setColorPhase( "Ki67+ (postmitotic)", new Color( 0, 128, 0 ) );
        //        visualizer.setColorPhase( "Apoptotic", Color.red );
        this.visualizers.add( visualizer );
        return visualizer;
    }

    public Model()
    {
        m = new Microenvironment();
    }

    public Model(Microenvironment m)
    {
        this.m = m;
    }

    public Microenvironment getMicroenvironment()
    {
        return m;
    }

    public void createContainer(double voxelSize)
    {
        CellContainer.createCellContainer( m, voxelSize );
    }

    public void simulate() throws Exception
    {
        double curTime = 0;
        double next_full_save_time = 0;
        boolean enable_legacy_saves = false;
        File reportFile;

        for( Visualizer listener : visualizers )
            listener.init();

        double startTime = System.currentTimeMillis();
        boolean hasEvents = !events.isEmpty();
        boolean eventsFired = false;
        try
        {
            while( curTime < tMax + 0.1 * diffusion_dt )
            {
                if( hasEvents )
                {
                    eventsFired = false;
                    Set<Event> executedEvens = new HashSet<>();
                    for( Event event : events )
                    {
                        if( curTime > event.executionTime - 0.01 * diffusion_dt )
                        {
                            event.execute( this );
                            executedEvens.add( event );
                        }
                    }
                    events.removeAll( executedEvens ); //events are one-time things
                    hasEvents = !events.isEmpty();
                }

                // save data if it's time. 
                if( Math.abs( curTime - next_full_save_time ) < 0.01 * diffusion_dt || eventsFired )
                {
                    if( enableFullSaves )
                    {
                        for( Visualizer listener : visualizers )
                            listener.saveResult( m, curTime );

                        String info = PhysiCellUtilities.getCurrentTime() + "\tElapsed\t"
                                + ( System.currentTimeMillis() - startTime ) / 1000 + "\tTime:\t" + (int)Math.round( curTime ) + "\tCells\t"
                                + m.getAgentsCount();

                        if( logFile != null )
                        {
                            try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( logFile ), true ) ))
                            {
                                bw.append( info );
                                bw.append( "\n" );
                            }
                        }
                        System.out.println( info );
                    }
                    writeReport( resultFolder + "/step_" + curTime + ".txt", "name\ti1\ti2\tp\tpressure\tphase\telapsed\n" );
                    for( Cell cell : m.getAgents( Cell.class ) )
                    {
                        String report = cell.typeName + "\t" + cell.get_current_mechanics_voxel_index() + "\t" + cell.currentVoxelIndex
                                + "\t"
                                + cell.nearest_density_vector()[0] + "\t" + cell.state.simplePressure + "\t"
                                + cell.phenotype.cycle.currentPhase().name + "\t" + cell.phenotype.cycle.data.elapsedTimePhase + "\t"
                                + cell.isOutOfDomain + "\n";
                    //                        if( cell.phenotype.cycle.code == 5 )
                    //                            report = cell.typeName + "\t" + cell.phenotype.cycle.currentPhase().name + "\t"
                    //                                    + cell.phenotype.cycle.transition_rate( 0, 0 ) + "\t" + cell.state.simplePressure + "\n";
                    //                        else
                    //                            report = cell.typeName + "\t" + cell.phenotype.cycle.currentPhase().name + "\t0.0\t"
                    //                                    + cell.state.simplePressure + "\n";
                    //
                    writeReport( resultFolder + "/step_" + curTime + ".txt", report );
                }
                    next_full_save_time += full_save_interval;
                }
                m.simulate_diffusion_decay( diffusion_dt );
                ( (CellContainer)m.agentContainer ).updateAllCells( m, curTime, phenotype_dt, mechanics_dt, diffusion_dt );
                curTime += diffusion_dt;
            }

            for( Visualizer listener : visualizers )
                listener.finish();
        }
        catch( Exception ex )
        {
            ex.printStackTrace();
        }
    }

    private void writeReport(String name, String str) throws Exception
    {
        File f = new File( name );
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( f, true ) ))
        {
            bw.append( str );
        }
    }


    public double[] getAverageTransition(Model model)
    {
        double sum = 0;
        int number = 0;
        for( Cell cell : model.getMicroenvironment().getAgents( Cell.class ) )
        {
            if( cell.phenotype.cycle.code != 5 )
                continue;
            number++;
            //            double prot = SignalBehavior.get_single_signal( cell, "custom:oncoprotein" );
            double rate = cell.phenotype.cycle.transition_rate( 0, 0 );
            sum += rate;
        }
        sum /= number;
        return new double[] {sum, number};
    }

    public void addParameter(String name, String val)
    {
        this.parameters.put( name, val );
    }

    public String getParameter(String name)
    {
        return parameters.get( name );
    }

    public int getParameterInt(String name)
    {
        return Integer.parseInt( parameters.get( name ) );
    }

    public double getParameterDouble(String name)
    {
        return Double.parseDouble( parameters.get( name ) );
    }

    public boolean getParameterBoolean(String name)
    {
        return Boolean.parseBoolean( parameters.get( name ) );
    }

    public void setTMax(double tMax)
    {
        this.tMax = tMax;
    }

    public void setDiffusionDt(double diffusion_dt)
    {
        this.diffusion_dt = diffusion_dt;
    }

    public void setMechanicsDt(double mechanics_dt)
    {
        this.mechanics_dt = mechanics_dt;
    }

    public void setPhenotypeDt(double phenotype_dt)
    {
        this.phenotype_dt = phenotype_dt;
    }

    public void setSaveInterval(double full_save_interval)
    {
        this.full_save_interval = full_save_interval;
    }

    public void setEnableFullSaves(boolean enable)
    {
        enableFullSaves = enable;
    }

    public boolean isEnableFullSaves()
    {
        return enableFullSaves;
    }

    public static abstract class Event
    {
        public double executionTime;
        public boolean executed = false;
        public abstract void execute(Model model) throws Exception;

        public Event(double executionTime)
        {
            this.executionTime = executionTime;
        }
    }

    public void saveResult()
    {

    }
}
