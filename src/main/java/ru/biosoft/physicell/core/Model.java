package ru.biosoft.physicell.core;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.standard.StandardModels;
import ru.biosoft.physicell.ui.AgentColorer;
import ru.biosoft.physicell.ui.Visualizer;
import ru.biosoft.physicell.ui.Visualizer2D;
import ru.biosoft.physicell.ui.Visualizer2D.Section;
import ru.biosoft.physicell.xml.ModelReader.ExternalFile;

public class Model
{
    protected Rules rules = new Rules();

    public SignalBehavior signals;

    //Random number generator for the model
    protected RandomGenerator rng;

    //Generates report for all cells at each time point of simulation
    protected ReportGenerator reportGenerator = new ReportGenerator();

    protected GlobalReportGenerator globalReportGenerator = new GlobalReportGenerator();

    private List<CellDefinition> cellDefinitions = new ArrayList<>();
    private Map<Integer, Integer> typeToIndex = new HashMap<>();
    private Map<String, CellDefinition> cellDefinitionNames = new HashMap<>();
    private Map<Integer, CellDefinition> cellDefinitionTypes = new HashMap<>();

    private List<Visualizer> visualizers = new ArrayList<Visualizer>();
    private String logFile;
    protected Microenvironment m;
    private Map<String, UserParameter> parameters = new HashMap<>();

    private double tMax;
    protected double startTime;
    protected double curTime = 0;
    protected double diffusionStep;
    protected double mechanicsStep = 0.1;
    protected double phenotypeStep = 6.0;
    protected double intracellularStep = 0.01;
    protected double nextIntracellularUpdate = 0;

    public boolean disableAutomatedSpringAdhesions = false;

    private boolean hasEvents = false;
    private List<Event> events = new ArrayList<>();

    private boolean verbose = true;
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

    private boolean rulesEnabled = false;

    public ExternalFile initialInfo = null;
    public ExternalFile reportInfo = null;
    public ExternalFile visualizerInfo = null;

    public static double tDiffusion = 0;

    public void setSeed(long seed)
    {
        this.rng.setSeed( seed );
    }

    public long getSeed()
    {
        return rng.getSeed();
    }

    public Rules getRules()
    {
        return rules;
    }

    public Iterable<Visualizer> getVisualizers()
    {
        return visualizers;
    }

    public SignalBehavior getSignals()
    {
        return signals;
    }

    public RandomGenerator getRNG()
    {
        return rng;
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
        //        rulesPath = path;
    }

    public void addEvent(Event event)
    {
        this.events.add( event );
        this.hasEvents = true;
    }

    public Visualizer addVisualizer(Visualizer visualizer)
    {
        this.visualizers.add( visualizer );
        return visualizer;
    }

    public Visualizer2D addGIFVisualizer(int zSlice, String name)
    {
        Visualizer2D visualizer = Visualizer2D.createWithGIF( resultFolder, name, Section.Z, zSlice );
        this.visualizers.add( visualizer );
        return visualizer;
    }

    public Model()
    {
        m = new Microenvironment();
        this.signals = new SignalBehavior();
        this.rng = new RandomGenerator();
    }

    public Model(Microenvironment m)
    {
        this.m = m;
        this.signals = new SignalBehavior();
        this.rng = new RandomGenerator();
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

    public void createContainer(double voxelSize, String name) throws Exception
    {
        CellContainer container = CellContainer.createCellContainer( m, name, voxelSize );
        container.setRulesEnabled( this.rulesEnabled );
    }

    public void setWriteDensity(boolean writeDensity)
    {
        this.saveDensity = writeDensity;
    }


    public void initFiles() throws Exception
    {
        if( resultFolder == null )
            return;
        File f = new File( resultFolder );
        if( f.exists() )
            deleteDirectory( f );
        if( saveDensity )
            new File( densityFolder ).mkdirs();
        if( saveReport )
            new File( cellDataFolder ).mkdirs();
        f.mkdirs();

        writeReport( modelFile, display() );
    }
    
    public void setupInitial() throws Exception
    {
        //Nothing by default
    }
    
    public void init(boolean outputFiles) throws Exception
    {
        curTime = 0;
        saveFullNext = 0;

        if( outputFiles )
            initFiles();

        for( Visualizer listener : visualizers )
            listener.init();

        startTime = System.currentTimeMillis();
        hasEvents = !events.isEmpty();

        signals.setupDictionaries( this );
        
        for (CellDefinition cd: getCellDefinitions())
            cd.distribution.apply(cd, this);
    }

    public void init() throws Exception
    {
        init( true );
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
        while( curTime < tMax + 0.1 * diffusionStep )
        {
            boolean eventsFired = executeEvents( curTime );
            if( ( Math.abs( curTime - saveFullNext ) < 0.01 * diffusionStep || eventsFired ) )
            {
                if( saveFull )
                    saveFull();
                if( verbose )
                    System.out.println( getLogInfo() );
                saveFullNext += saveFullInterval;
            }
            if( saveImg && ( Math.abs( curTime - saveImgNext ) < 0.01 * diffusionStep || eventsFired ) )
            {
                saveImg();
                saveImgNext += saveImgInterval;
            }

            doStep();
        }

        for( Visualizer listener : visualizers )
            listener.finish();
    }

    public double getCurrentTime()
    {
        return curTime;
    }
    
    public void stepBeforeCells()
    {
        //do nothing by default
    }
    
    public void stepAfterCells() throws Exception
    {
        //do nothing by default
    }

    public void doStep() throws Exception
    {
        double tDiff = System.nanoTime();
        m.simulateDiffusionDecay( diffusionStep );
        tDiff = System.nanoTime() - tDiff;
        tDiffusion += tDiff;
        stepBeforeCells();
        ( (CellContainer)m.agentContainer ).updateAllCells( this, curTime, phenotypeStep, mechanicsStep, diffusionStep );
        stepAfterCells();
        updateIntracellular();
        curTime += diffusionStep;
        m.time = curTime;
    }

    public void setRulesEnabled(boolean rulesEnabled)
    {
        this.rulesEnabled = rulesEnabled;
        if( m.agentContainer instanceof CellContainer )
            ( (CellContainer)m.agentContainer ).setRulesEnabled( rulesEnabled );
    }

    public void updateIntracellular() throws Exception
    {
        if( curTime >= nextIntracellularUpdate )
        {
            m.getAgents( Cell.class ).parallelStream().filter( cell -> !cell.isOutOfDomain && !cell.phenotype.death.dead )
                    .forEach( cell -> {
                        try
                        {
                            Intracellular intra = cell.phenotype.intracellular;
                            if( intra != null )
                            {
                                intra.updateIntracellularParameters( m, cell.phenotype );
                                intra.step();
                                intra.updatePhenotypeParameters( m, cell.phenotype );
                            }
                        }
                        catch( Exception ex )
                        {
                            ex.printStackTrace();
                        }
                    } );
            nextIntracellularUpdate += intracellularStep;
        }
    }

    private void saveImg() throws Exception
    {
        for( Visualizer listener : visualizers )
            listener.saveResult( m, curTime );
    }

    public void saveFull() throws Exception
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
        //        System.out.println( getLogInfo() );
    }

    public boolean executeEvents(double curTime) throws Exception
    {
        boolean eventsFired = false;
        if( hasEvents )
        {
            Set<Event> executedEvens = new HashSet<>();
            for( Event event : events )
            {
                if( curTime > event.executionTime - 0.01 * diffusionStep )
                {
                    event.execute();
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


    public void addParameter(String name, String val, String description)
    {
        this.parameters.put( name, new UserParameter( name, description, val ) );
    }

    public Set<String> getParameters()
    {
        return parameters.keySet();
    }

    public UserParameter getParameter(String name)
    {
        return parameters.get( name );
    }

    public String getParameterString(String name)
    {
        return parameters.get( name ).getValue();
    }

    public int getParameterInt(String name)
    {
        return Integer.parseInt( parameters.get( name ).getValue() );
    }

    public double getParameterDouble(String name)
    {
        return Double.parseDouble( parameters.get( name ).getValue() );
    }

    public boolean getParameterBoolean(String name)
    {
        return Boolean.parseBoolean( parameters.get( name ).getValue() );
    }

    public void setTMax(double tMax)
    {
        this.tMax = tMax;
    }
    public double getTMax()
    {
        return tMax;
    }

    public void setDiffusionDt(double dt)
    {
        this.diffusionStep = dt;
    }
    public double getDiffusionDt()
    {
        return diffusionStep;
    }

    public void setMechanicsDt(double dt)
    {
        this.mechanicsStep = dt;
    }

    public double getMechanicsDt()
    {
        return mechanicsStep;
    }

    public void setPhenotypeDt(double dt)
    {
        this.phenotypeStep = dt;
    }

    public double getPhenotypeDt()
    {
        return phenotypeStep;
    }

    public void setSaveFullInterval(double interval)
    {
        this.saveFullInterval = interval;
    }

    public double getSaveFullInterval()
    {
        return saveFullInterval;
    }

    public void setSaveImgInterval(double interval)
    {
        this.saveImgInterval = interval;
    }

    public double getSaveImgInterval()
    {
        return saveImgInterval;
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
        protected Model model;
        public boolean executed = false;
        public abstract void execute() throws Exception;

        public Event(Model model)
        {
            this.model = model;
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
        sb.append( "\tSeed " + getSeed() );
        sb.append( "\tCell update: " + ( (CellContainer)m.agentContainer ).getName() );
        sb.append( "\tDiffusion: " + m.getSolver().getName() );
        sb.append( "\n\n" + getMicroenvironment().display() );
        sb.append( "\n" );
        sb.append( "\nCell Types: ( " + getDefinitionsCount() + " total)" );
        sb.append( "\n--------------------------------" );
        for( int i = 0; i < getDefinitionsCount(); i++ )
            sb.append( "\n\t" + i + ". " + getCellDefinitionByIndex( i ).name + " # " + calcCells( getCellDefinitionByIndex( i ) ) );

        for( CellDefinition cd : getCellDefinitions() )
            sb.append( "\n\n" + cd.display() );
        sb.append( "\n\n================================" );
        sb.append( "\nGlobal parameters" );
        sb.append( "\n================================" );
        for( Entry<String, UserParameter> e : parameters.entrySet() )
            sb.append( "\n\t" + e.getKey() + " " + e.getValue().getValue() );

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

    public ExternalFile getVisualizerInfo()
    {
        return visualizerInfo;
    }
    public void setVisualizerInfo(ExternalFile visualizerInfo)
    {
        this.visualizerInfo = visualizerInfo;
    }

    public ExternalFile getReportInfo()
    {
        return reportInfo;
    }
    public void setReportInfo(ExternalFile reportInfo)
    {
        this.reportInfo = reportInfo;
    }

    public void setInitialInfo(ExternalFile info)
    {
        this.initialInfo = info;
    }
    public ExternalFile getInitialInfo()
    {
        return initialInfo;
    }

    public void setReportGenerator(ReportGenerator reportGenerator)
    {
        this.reportGenerator = reportGenerator;
    }
    public ReportGenerator getReportGenerator()
    {
        return reportGenerator;
    }

    public void setReportGenerator(GlobalReportGenerator globalReportGenerator)
    {
        this.globalReportGenerator = globalReportGenerator;
    }
    public GlobalReportGenerator getGlobalReportGenerator()
    {
        return globalReportGenerator;
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

    public String getLogInfo() throws Exception
    {
        return PhysiCellUtilities.getCurrentTime() + "\tElapsed\t" + ( System.currentTimeMillis() - startTime ) / 1000 + "\tTime:\t"
                + (int)Math.round( curTime ) + "\tCells\t" + m.getAgentsCount();
    }

    public void setVerbose(boolean verbose)
    {
        this.verbose = verbose;
    }

    public int getDefinitionsCount()
    {
        return cellDefinitions.size();
    }

    public void registerCellDefinition(CellDefinition cd)
    {
        cellDefinitionNames.put( cd.name, cd );
        cellDefinitionTypes.put( cd.type, cd );
        typeToIndex.put( cd.type, cellDefinitions.size() );
        cellDefinitions.add( cd );
        sync();
    }

    public void clearCellDefinitions()
    {
        cellDefinitions.clear();
        cellDefinitionNames.clear();
        cellDefinitionTypes.clear();
        typeToIndex.clear();
        StandardModels.defaults = null;
        sync();
    }

    public Iterable<CellDefinition> getCellDefinitions()
    {
        return cellDefinitions;
    }

    public CellDefinition getCellDefinitionByIndex(int index)
    {
        return cellDefinitions.get( index );
    }

    public int getCellDefinitionIndex(int type)
    {
        return typeToIndex.get( type );
    }

    public CellDefinition getCellDefinition(int type)
    {
        return cellDefinitionTypes.get( type );
    }

    public Set<String> getCellDefinitionNames()
    {
        return cellDefinitionNames.keySet();
    }

    public CellDefinition getCellDefinition(String name)
    {
        return cellDefinitionNames.get( name );
    }

    private void sync()
    {
        for( CellDefinition cd : cellDefinitions )
        {
            cd.phenotype.sync( this );
        }
    }

    public int findCellDefinitionIndex(String name)
    {
        CellDefinition cd = cellDefinitionNames.get( name );
        if( cd != null )
            return cd.type;
        return -1;
    }
    
    public AgentColorer getDefaultColorer()
    {
        return null;
    }
}