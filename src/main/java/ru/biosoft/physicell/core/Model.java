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
    protected Microenvironment m;
    private Map<String, String> parameters = new HashMap<>();

    private double tMax;
    protected double startTime;
    protected double curTime = 0;
    protected double diffusion_dt;
    protected double mechanics_dt = 0.1;
    protected double phenotype_dt = 6.0;
    protected double intracellular_dt = 0.01;
    protected double next_intracellular_update = intracellular_dt;

    private boolean hasEvents = false;
    private List<Event> events = new ArrayList<>();

    private boolean saveFull = true;
    private double saveFullInterval;
    private double saveFullNext = 0;
    private boolean saveDensity = true;
    private boolean saveReport = true;
    private boolean saveImg = true;
    private double saveImgNext = 0;
    private double saveImgInterval;

    private String resultFolder;
    private String densityFolder = null;
    private String cellDataFolder = null;
    private String modelFile = null;

    private String initialPath = null;
    private boolean rulesEnabled = false;
    private String rulesPath = null;

    public Iterable<Visualizer> getVisualizers()
    {
        return visualizers;
    }

    public void setResultFolder(String folder)
    {
        this.resultFolder = folder;
        this.densityFolder = folder + "/density";
        this.cellDataFolder = folder + "/cells";
        this.logFile = folder + "/log.txt";
        this.modelFile = folder + "/model.txt";
    }

    public void setRulesPath(String path)
    {
        rulesEnabled = true;
        rulesPath = path;
    }

    public void addEvent(Event event)
    {
        this.events.add( event );
        this.hasEvents = true;
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

    public void createContainer()
    {
        createContainer( 30 );

    }

    public void createContainer(double voxelSize)
    {
        CellContainer container = CellContainer.createCellContainer( m, voxelSize );
        container.setRulesEnabled( this.rulesEnabled );
    }

    public void setWriteDensity(boolean writeDensity)
    {
        this.saveDensity = writeDensity;
    }

    public void init() throws Exception
    {
        curTime = 0;
        saveFullNext = 0;

        File f = new File( resultFolder );
        if( f.exists() )
            deleteDirectory( f );
        if( saveDensity )
            new File( densityFolder ).mkdirs();
        if( saveReport )
            new File( cellDataFolder ).mkdirs();
        f.mkdirs();

        writeReport( modelFile, display() );

        for( Visualizer listener : visualizers )
            listener.init();

        startTime = System.currentTimeMillis();
        hasEvents = !events.isEmpty();
    }

    public static void deleteDirectory(File directoryToBeDeleted)
    {
        File[] allContents = directoryToBeDeleted.listFiles();
        if( allContents != null )
        {
            for( File file : allContents )
            {
                deleteDirectory( file );
            }
        }
        directoryToBeDeleted.delete();
    }

    public void simulate() throws Exception
    {
        while( curTime < tMax + 0.1 * diffusion_dt )
        {
            boolean eventsFired = executeEvents( curTime );
            if( saveFull && ( Math.abs( curTime - saveFullNext ) < 0.01 * diffusion_dt || eventsFired ) )
            {
                saveFull();
                saveFullNext += saveFullInterval;
            }
            if( saveImg && ( Math.abs( curTime - saveImgNext ) < 0.01 * diffusion_dt || eventsFired ) )
            {
                saveImg();
                saveImgNext += saveImgInterval;
            }
            doStep();
        }

        for( Visualizer listener : visualizers )
            listener.finish();
    }

    protected void doStep() throws Exception
    {
        m.simulateDiffusionDecay( diffusion_dt );
        ( (CellContainer)m.agentContainer ).updateAllCells( m, curTime, phenotype_dt, mechanics_dt, diffusion_dt );
        if( curTime >= next_intracellular_update )
        {
            updateIntracellular();
            next_intracellular_update += intracellular_dt;
        }
        curTime += diffusion_dt;
        m.time = curTime;
    }

    private void saveImg() throws Exception
    {
        for( Visualizer listener : visualizers )
            listener.saveResult( m, curTime );
    }

    private void saveFull() throws Exception
    {
        if( logFile != null )
        {
            try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( logFile ), true ) ))
            {
                bw.append( getLogInfo() );
                bw.append( "\n" );
            }
        }
        if( saveDensity )
        {
            getMicroenvironment().writeDensity( densityFolder + "/Density " + Math.round( curTime ) + ".txt" );
        }
        if( saveReport )
        {
            writeReport( cellDataFolder + "/Cells " + Math.round( curTime ) + ".txt", this.getReportHeader() );
            for( Cell cell : getMicroenvironment().getAgents( Cell.class ) )
            {
                writeReport( cellDataFolder + "/Cells " + Math.round( curTime ) + ".txt", this.getReport( cell ) );
            }
        }
        System.out.println( getLogInfo() );
    }

    private boolean executeEvents(double curTime) throws Exception
    {
        boolean eventsFired = false;
        if( hasEvents )
        {
            Set<Event> executedEvens = new HashSet<>();
            for( Event event : events )
            {
                if( curTime > event.executionTime - 0.01 * diffusion_dt )
                {
                    event.execute( this );
                    executedEvens.add( event );
                    eventsFired = true;
                }
            }
            events.removeAll( executedEvens ); //events are one-time things
            hasEvents = !events.isEmpty();
        }
        return eventsFired;
    }

    private void writeReport(String name, String str) throws Exception
    {
        File f = new File( name );
        if( !f.exists() )
            f.createNewFile();
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( f, true ) ))
        {
            bw.append( str );
        }
    }

    public void printMaxInterferon(Model model)
    {
        double sum = 0;
        for( int i = 0; i < model.getMicroenvironment().numberVoxels(); i++ )
            sum += model.getMicroenvironment().getDensity( i )[0];
        System.out.println( sum );
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

    public void setDiffusionDt(double dt)
    {
        this.diffusion_dt = dt;
    }

    public void setMechanicsDt(double dt)
    {
        this.mechanics_dt = dt;
    }

    public void setPhenotypeDt(double dt)
    {
        this.phenotype_dt = dt;
    }

    public void setSaveFullInterval(double interval)
    {
        this.saveFullInterval = interval;
    }

    public void setSaveImgInterval(double interval)
    {
        this.saveImgInterval = interval;
    }

    public void setSaveFull(boolean enable)
    {
        saveFull = enable;
    }

    public void setSaveImg(boolean enable)
    {
        saveImg = enable;
    }

    public boolean isEnableFullSaves()
    {
        return saveFull;
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

    public String display()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "================================" );
        sb.append( "\nSimulation Options" );
        sb.append( "\n================================" );
        sb.append( "\n\tMaximum Time: " + tMax );
        sb.append( "\tSave interval " + saveFullInterval );
        sb.append( "\tSeed " + PhysiCellUtilities.getSeed() );
        sb.append( "\n\n" + getMicroenvironment().display() );
        sb.append( "\n" );
        sb.append( "\nCell Types: ( " + CellDefinition.getDefinitionsCount() + " total)" );
        sb.append( "\n--------------------------------" );
        for( int i = 0; i < CellDefinition.getDefinitionsCount(); i++ )
            sb.append( "\n\t" + i + ". " + CellDefinition.getCellDefinitionByIndex( i ).name + " # "
                    + calcCells( CellDefinition.getCellDefinitionByIndex( i ) ) );

        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
            sb.append( "\n\n" + cd.display() );
        return sb.toString();
    }

    public int calcCells(CellDefinition cd)
    {
        int result = 0;
        for( Cell cell : m.getAgents( Cell.class ) )
            if( cell.typeName.equals( cd.name ) )
                result++;
        return result;
    }

    public void setInitialPath(String path)
    {
        this.initialPath = path;
    }

    public String getInitialPath()
    {
        return initialPath;
    }

    public String getReport(Cell cell) throws Exception
    {
        return "\n" + cell.ID + "\t" + cell.position[0] + "\t" + cell.position[1] + "\t" + cell.position[2] + "\t"
                + cell.phenotype.cycle.currentPhase().name + "\t" + cell.phenotype.cycle.data.elapsedTimePhase;
    }

    public String getReportHeader()
    {
        return "ID\tX\tY\tZ\tCycle\tElapsed";
    }

    public void updateIntracellular() throws Exception
    {

    }

    public String getLogInfo() throws Exception
    {
        return PhysiCellUtilities.getCurrentTime() + "\tElapsed\t" + ( System.currentTimeMillis() - startTime ) / 1000 + "\tTime:\t"
                + (int)Math.round( curTime ) + "\tCells\t" + m.getAgentsCount();
    }
}