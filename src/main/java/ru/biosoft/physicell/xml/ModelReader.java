package ru.biosoft.physicell.xml;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.MicroenvironmentOptions;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.CellCSVReader;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellInteractions;
import ru.biosoft.physicell.core.CellTransformations;
import ru.biosoft.physicell.core.CycleModel;
import ru.biosoft.physicell.core.Death;
import ru.biosoft.physicell.core.DeathParameters;
import ru.biosoft.physicell.core.Mechanics;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Motility;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.PhysiCellSettings;
import ru.biosoft.physicell.core.Rules;
import ru.biosoft.physicell.core.Secretion;
import ru.biosoft.physicell.core.Volume;
import ru.biosoft.physicell.core.standard.AdvancedChemotaxis;
import ru.biosoft.physicell.core.standard.AdvancedChemotaxisNormalized;
import ru.biosoft.physicell.core.standard.Chemotaxis;
import ru.biosoft.physicell.core.standard.DomainEdgeAvoidance;
import ru.biosoft.physicell.core.standard.StandardModels;


public class ModelReader extends ModelReaderSupport
{
    private IntracellularReader intracellularReader = null;
    private FunctionsReader functionsReader = null;
    private static final String OPTIONS = "options";
    private Map<String, File> additionalFiles = new HashMap<>();

    private boolean readFromJAR = false;
    //folder from which we read file
    private Path filePath = null;

    public Model read(InputStream is) throws Exception
    {
        return this.read( is, null );
    }

    public Model read(InputStream is, Class<?> clazz) throws Exception
    {
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        DocumentBuilder builder = factory.newDocumentBuilder();
        Document doc = builder.parse( is );
        return read( doc, clazz );
    }

    public Model read(File f, Class<?> clazz) throws Exception
    {
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        DocumentBuilder builder = factory.newDocumentBuilder();
        Document doc = builder.parse( f );
        filePath = f.toPath();
        return read( doc, clazz );
    }

    public Model read(Document doc, Class<?> clazz) throws Exception
    {
        NodeList nodes = doc.getChildNodes();

        Element physicell = findElement( nodes, PHYSICELL_ELEMENT );
        if( physicell == null )
            throw new Exception( "Physicell base element not found" );

        Model model = clazz != null ? (Model)clazz.getDeclaredConstructor().newInstance() : null;
        Microenvironment m = model.getMicroenvironment();
        readDomain( physicell, m );
        readOverall( physicell, model );
        readOptions( physicell );
        readSave( physicell, model );
        readMicroenvironmentSetup( physicell, m );
        readCellDefinitions( physicell, m, model );
        readPhenotypes( physicell, model );
        readInitialConditions( physicell, model );
        readUserParameters( physicell, model );
        readRules( physicell, model );
        return model;
    }

    public void setAdditionalFiles(Map<String, File> additionalFiles)
    {
        this.additionalFiles = additionalFiles;
    }

    public void setReadFromJAR(boolean readFromJAR)
    {
        this.readFromJAR = readFromJAR;
    }

    private void readRules(Element physicell, Model model) throws Exception
    {
        Element rulesElement = findElement( physicell, "cell_rules" );
        if( rulesElement == null )
            return;
        Element rulesetsElement = findElement( rulesElement, "rulesets" );
        for( Element rulesetElement : findAllElements( rulesetsElement, "ruleset" ) )
        {
            boolean enabled = getBoolAttr( rulesetElement, "enabled" );
            if( !enabled )
                continue;
//            String format = getAttr( rulesetElement, "format" );
//            String version = getAttr( rulesetElement, "version" );
//            String protocol = getAttr( rulesetElement, "CBHG" );
            Element folderElement = findElement( rulesetElement, "folder" );
            String folder = getVal( folderElement );
            Element filenameElement = findElement( rulesetElement, "filename" );
            String filename = getVal( filenameElement );

            InputStream is = null;

            if( folder.startsWith( "./" ) )
                folder = folder.substring( 2 );
            String relativePath = folder + "/" + filename;
            if( additionalFiles.containsKey( relativePath ) )
            {
                is = new FileInputStream( additionalFiles.get( relativePath ) );
            }
            else if( filePath != null )
            {
                //Path rulesPath = filePath.getParent().resolve( folder, filename );
                //is = new FileInputStream( rulesPath.toFile() );
                is = new FileInputStream( new File( folder, filename ) );
            }

            model.getSignals().setupDictionaries( model );
            Rules.setupRules( model );
            Rules.parseCSVRules2( model, is );
            model.setRulesPath( folder + "/" + filename );
        }
    }

    private void readDomain(Element physicell, Microenvironment m)
    {
        Element domainElement = findElement( physicell, DOMAIN_ELEMENT );
        double minX = 0;
        double maxX = 100;
        double minY = 0;
        double maxY = 100;
        double minZ = 0;
        double maxZ = 100;
        double dx = 10;
        double dy = 10;
        double dz = 10;
        boolean use2D = false;
        for( Element el : getAllElements( domainElement ) )
        {
            switch( el.getTagName() )
            {
                case X_MIN:
                    minX = getDoubleVal( el );
                    break;
                case X_MAX:
                    maxX = getDoubleVal( el );
                    break;
                case Y_MIN:
                    minY = getDoubleVal( el );
                    break;
                case Y_MAX:
                    maxY = getDoubleVal( el );
                    break;
                case Z_MIN:
                    minZ = getDoubleVal( el );
                    break;
                case Z_MAX:
                    maxZ = getDoubleVal( el );
                    break;
                case DX:
                    dx = getDoubleVal( el );
                    break;
                case DY:
                    dy = getDoubleVal( el );
                    break;
                case DZ:
                    dz = getDoubleVal( el );
                    break;
                case USE_2D:
                    use2D = getBoolVal( el );
                    break;
            }
        }
        m.resizeSpace( minX, maxX, minY, maxY, minZ, maxZ, dx, dy, dz );
        m.options.simulate2D = use2D;
        if( m.options.simulate2D )
        {
            minZ = -0.5 * dz;
            maxZ = 0.5 * dz;
        }
        m.options.X_range = new double[] {minX, maxX};
        m.options.Y_range = new double[] {minY, maxY};
        m.options.Z_range = new double[] {minZ, maxZ};

        m.options.dx = dx;
        m.options.dy = dy;
        m.options.dz = dz;
    }

    public void readOverall(Element physicell, Model model)
    {
        Microenvironment m = model.getMicroenvironment();
        double maxTime = 100;
        //        String maxTimeUnits = "min";
        String timeUnits = "min";
        String spaceUnits = "micron";
        //        String diffusionUnits = "min";
        //        String mechanicsUnits = "min";
        //        String phenotypeUnits = "min";
        double diffusionStep = 0.01;
        double mechanicsStep = 0.1;
        double phenotypeSteps = 6;
        Element overallElement = findElement( physicell, OVERALL );
        for( Element el : getAllElements( overallElement ) )
        {
            switch( el.getTagName() )
            {
                case MAX_TIME:
                    maxTime = getDoubleVal( el );
                    //                    maxTimeUnits = el.getAttribute( UNITS );
                    break;
                case TIME_UNITS:
                    timeUnits = getVal( el );
                    break;
                case SPACE_UNITS:
                    spaceUnits = getVal( el );
                    break;
                case DT_DIFFUSION:
                    //                    mechanicsUnits = el.getAttribute( UNITS );
                    diffusionStep = getDoubleVal( el );
                    break;
                case DT_MECHANICS:
                    //                    diffusionUnits = el.getAttribute( UNITS );
                    mechanicsStep = getDoubleVal( el );
                    break;
                case DT_PHENOTYPE:
                    //                    phenotypeUnits = el.getAttribute( UNITS );
                    phenotypeSteps = getDoubleVal( el );
                    break;

            }
        }
        m.timeUnits = timeUnits;
        m.spatialUnits = spaceUnits;
        model.setDiffusionDt( diffusionStep );
        model.setPhenotypeDt( phenotypeSteps );
        model.setMechanicsDt( mechanicsStep );
        model.setTMax( maxTime );
    }


    public void readSave(Element physicell, Model model)
    {
        Element saveElement = findElement( physicell, "save" );
        for( Element el : getAllElements( saveElement ) )
        {
            switch( el.getTagName() )
            {
                case "folder":
//                    String folderName = getVal( el );
                    break;
                case "full_data":
                    Element intervalElement = findElement( el, "interval" );
                    Element enabledElement = findElement( el, "enable" );
                    if( enabledElement != null )
                    {
                        model.setSaveFull( getBoolVal( enabledElement ) );
                    }
                    if( intervalElement != null && model.isEnableFullSaves() )
                    {
                        double interval = getDoubleVal( intervalElement );
                        //                        String intervalUnits = getAttr( intervalElement, "units" );
                        model.setSaveFullInterval( interval );
                    }

                    break;
                case "SVG":
                    intervalElement = findElement( el, "interval" );
                    enabledElement = findElement( el, "enable" );
                    if( enabledElement != null )
                    {
                        model.setSaveImg( getBoolVal( enabledElement ) );
                    }
                    if( intervalElement != null )
                    {
                        double interval = getDoubleVal( intervalElement );
                        //                        String intervalUnits = getAttr( intervalElement, UNITS );
                        model.setSaveImgInterval( interval );
                    }

                    break;
                case "legacy_data":
                    enabledElement = findElement( el, "" );
                    if( enabledElement != null )
                    {
                        //                        enabled = getBoolVal( enabledElement );
                    }
                    break;
                case "report":
                    boolean enabled = getBoolAttr( el, "enabled" );
                    if( enabled )
                    {
                        String folder = getVal( findElement( el, "folder" ) );
                        String filename = getVal( findElement( el, "filename" ) );
                        String format = getAttr( el, "format" );
                        String path = folder + "/" + filename;
                        model.setReportInfo( new ExternalFile( format, path ) );
                    }
                    break;
                case "visualizer":
                    enabled = getBoolAttr( el, "enabled" );//findElement( el, "enable" );
                    if( enabled )
                    {
                        String folder = getVal( findElement( el, "folder" ) );
                        String filename = getVal( findElement( el, "filename" ) );
                        String format = getAttr( el, "format" );
                        String path = folder + "/" + filename;
                        model.setVisualizerInfo( new ExternalFile( format, path ) );
                    }
                    break;
            }
        }
    }

    public void readOptions(Element physicell)
    {
        Element optionsElement = findElement( physicell, OPTIONS );
        for( Element el : getAllElements( optionsElement ) )
        {
            switch( el.getTagName() )
            {
                //                case "legacy_random_points_on_sphere_in_divide":
                //                    boolean legacy_random_points_on_sphere_in_divide = getBoolVal( el );
                //                    PhysiCellSettings.cell_division_orientation = LegacyRandomOnUnitSphere;
                //                    break;
                case "virtual_wall_at_domain_edge":
                {
                    try
                    {
                        boolean virtual_wall_at_domain_edge = getBoolVal( el );
                        if( virtual_wall_at_domain_edge )
                            StandardModels.getDefaultCellDefinition().functions.membraneInteraction = new DomainEdgeAvoidance();
                        break;
                    }
                    catch( Exception ex )
                    {

                    }
                }
                case "disable_automated_spring_adhesions":
                    PhysiCellSettings.disable_automated_spring_adhesions = getBoolVal( el );
                    break;
            }
        }
    }

    public void readMicroenvironmentSetup(Element physicell, Microenvironment m)
    {
        Element microenvironmentSetupElement = findElement( physicell, "microenvironment_setup" );
        List<Element> variableElements = findAllElements( microenvironmentSetupElement, "variable" );

        boolean activated_Dirichlet_boundary_detected = false;

        int varSize = variableElements.size();
        double[] initialValues = new double[varSize];
        double[] Dirichlet_condition = new double[varSize];
        boolean[] Dirichlet_activation = new boolean[varSize];

        boolean[] Dirichlet_xmin = new boolean[varSize];
        boolean[] Dirichlet_xmax = new boolean[varSize];
        boolean[] Dirichlet_ymin = new boolean[varSize];
        boolean[] Dirichlet_ymax = new boolean[varSize];
        boolean[] Dirichlet_zmin = new boolean[varSize];
        boolean[] Dirichlet_zmax = new boolean[varSize];

        double[] Dirichlet_xmin_values = new double[varSize];
        double[] Dirichlet_xmax_values = new double[varSize];
        double[] Dirichlet_ymin_values = new double[varSize];
        double[] Dirichlet_ymax_values = new double[varSize];
        double[] Dirichlet_zmin_values = new double[varSize];
        double[] Dirichlet_zmax_values = new double[varSize];
        boolean[] Dirichlet_all = new boolean[varSize];

        for( int i = 0; i < varSize; i++ )
        {
            Element variableElement = variableElements.get( i );
            String name = getAttr( variableElement, "name" );
            String units = getAttr( variableElement, UNITS );
            int id = getIntAttr( variableElement, "ID" );
            Element physicalElement = findElement( variableElement, "physical_parameter_set" );
            Element diffusionElement = findElement( physicalElement, "diffusion_coefficient" );
            //            String diffusionUnits = getAttr( diffusionElement, UNITS );
            double diffusionValue = getDoubleVal( diffusionElement );
            Element decayElement = findElement( physicalElement, "decay_rate" );
            //            String decayUnits = getAttr( diffusionElement, UNITS );
            double decayValue = getDoubleVal( decayElement );

            if( id == 0 )
                m.setDensity( id, name, units, diffusionValue, decayValue );
            else
                m.addDensity( name, units, diffusionValue, decayValue );

            Element initialConditionElement = findElement( variableElement, "initial_condition" );
            //            String initialUnits = initialConditionElement.getAttribute( UNITS );
            initialValues[i] = getDoubleVal( initialConditionElement );

            Element DirichletBoundaryElement = findElement( variableElement, "Dirichlet_boundary_condition" );
            //            String dirichletUnits = getAttr( DirichletBoundaryElement, UNITS );
            boolean dirichletEnabled = getBoolAttr( DirichletBoundaryElement, "enabled" );
            double dirichletValue = getDoubleVal( DirichletBoundaryElement );
            Dirichlet_condition[i] = dirichletValue;
            Dirichlet_activation[i] = dirichletEnabled;
            Dirichlet_all[i] = dirichletEnabled;

            if( dirichletEnabled )
                activated_Dirichlet_boundary_detected = true;

            Dirichlet_xmin[i] = dirichletEnabled;
            Dirichlet_xmax[i] = dirichletEnabled;
            Dirichlet_ymin[i] = dirichletEnabled;
            Dirichlet_ymax[i] = dirichletEnabled;
            Dirichlet_zmin[i] = dirichletEnabled;
            Dirichlet_zmax[i] = dirichletEnabled;

            Dirichlet_xmin_values[i] = dirichletValue;
            Dirichlet_xmax_values[i] = dirichletValue;
            Dirichlet_ymin_values[i] = dirichletValue;
            Dirichlet_ymax_values[i] = dirichletValue;
            Dirichlet_zmin_values[i] = dirichletValue;
            Dirichlet_zmax_values[i] = dirichletValue;

            Element dirichletOptionsElement = findElement( variableElement, "Dirichlet_options" );
            if( dirichletOptionsElement != null )
            {
                for( Element el : findAllElements( dirichletOptionsElement, "boundary_value" ) )
                {
                    boolean enabled = getBoolAttr( el, "enabled" );
                    if( !enabled )
                        Dirichlet_all[i] = false;
                    switch( getAttr( el, "ID" ) )
                    {
                        case "xmin":
                            Dirichlet_xmin[i] = enabled;
                            if( enabled )
                                Dirichlet_xmin_values[i] = getDoubleVal( el );
                            break;
                        case "xmax":
                            Dirichlet_xmax[i] = enabled;
                            if( enabled )
                                Dirichlet_xmax_values[i] = getDoubleVal( el );
                            break;
                        case "ymin":
                            Dirichlet_ymin[i] = enabled;
                            if( enabled )
                                Dirichlet_ymin_values[i] = getDoubleVal( el );
                            break;
                        case "ymax":
                            Dirichlet_ymax[i] = enabled;
                            if( enabled )
                                Dirichlet_ymax_values[i] = getDoubleVal( el );
                            break;
                        case "zmin":
                            Dirichlet_zmin[i] = enabled;
                            if( enabled )
                                Dirichlet_zmin_values[i] = getDoubleVal( el );
                            break;
                        case "zmax":
                            Dirichlet_zmax[i] = enabled;
                            if( enabled )
                                Dirichlet_zmax_values[i] = getDoubleVal( el );
                            break;
                    }
                }
            }
        }
        MicroenvironmentOptions options = m.getOptions();
        options.outer_Dirichlet_conditions = false;

        options.Dirichlet_condition_vector = Dirichlet_condition;
        options.Dirichlet_activation_vector = Dirichlet_activation;
        options.initial_condition_vector = initialValues;

        options.Dirichlet_all = Dirichlet_all;

        options.Dirichlet_xmin = Dirichlet_xmin;
        options.Dirichlet_xmax = Dirichlet_xmax;
        options.Dirichlet_ymin = Dirichlet_ymin;
        options.Dirichlet_ymax = Dirichlet_ymax;
        options.Dirichlet_zmin = Dirichlet_zmin;
        options.Dirichlet_zmax = Dirichlet_zmax;

        options.Dirichlet_xmin_values = Dirichlet_xmin_values;
        options.Dirichlet_xmax_values = Dirichlet_xmax_values;
        options.Dirichlet_ymin_values = Dirichlet_ymin_values;
        options.Dirichlet_ymax_values = Dirichlet_ymax_values;
        options.Dirichlet_zmin_values = Dirichlet_zmin_values;
        options.Dirichlet_zmax_values = Dirichlet_zmax_values;

        for( int i = 0; i < m.numberVoxels(); i++ )
            m.density[i] = options.initial_condition_vector.clone();

        // if any of the substrates have outer Dirichlet conditions enables, then set the outer_Dirichlet_conditions = true;       
        if( activated_Dirichlet_boundary_detected )
            options.outer_Dirichlet_conditions = true;

        Element optionsElement = findElement( microenvironmentSetupElement, "options" );
        if( options != null )
        {
            Element calculate_gradientsElement = findElement( optionsElement, "calculate_gradients" );
            if( calculate_gradientsElement != null )
                options.calculate_gradients = getBoolVal( calculate_gradientsElement );

            Element trackSubstrate = findElement( optionsElement, "track_internalized_substrates_in_each_agent" );
            if( trackSubstrate != null )
                options.track_internalized_substrates_in_each_agent = getBoolVal( trackSubstrate );
        }
        Microenvironment.initialize( m );
    }

    public void readCellDefinitions(Element physicell, Microenvironment m, Model model) throws Exception
    {
        Element cellDefinitionsElement = findElement( physicell, "cell_definitions" );
        if( cellDefinitionsElement == null )
            return;
        for( Element cdElement : findAllElements( cellDefinitionsElement, "cell_definition" ) )
        {
            CellDefinition cd;
            String name = getAttr( cdElement, "name" );
            Integer ID = getIntAttr( cdElement, "ID" );
            int id = ID == null ? -1 : ID;
            boolean defaultDefinition = name.equals( "default" );
            if( defaultDefinition ) //TODO: check ID=0 
                cd = StandardModels.createFromDefault( name, id, m );
            else
                cd = new CellDefinition( m, id, name );

            CellDefinition parent = null;
            String parentType = getAttr( cdElement, "parent_type" );

            if( parentType.isEmpty() && !defaultDefinition )
            {
                parent = StandardModels.createFromDefault( "default", -1, m );
                cd = parent.clone( name, id, m );
            }
            Element functionsElement = findElement( cdElement, "functions" );
            if( functionsElement != null )
                readFunctions( functionsElement, cd );

            model.registerCellDefinition( cd );
        }
    }

    public void readPhenotypes(Element physicell, Model model) throws Exception
    {
        Microenvironment m = model.getMicroenvironment();
        Element cellDefinitionsElement = findElement( physicell, "cell_definitions" );
        if( cellDefinitionsElement == null )
            return;
        for( Element cdElement : findAllElements( cellDefinitionsElement, "cell_definition" ) )
        {
            String name = getAttr( cdElement, "name" );
            CellDefinition cd = model.getCellDefinition( name );
            Phenotype p = cd.phenotype;

            CellDefinition parent = null;
            String parentType = getAttr( cdElement, "parent_type" );
            if( !parentType.isEmpty() )
                parent = model.getCellDefinition( parentType );

            if( parent != null )
            {
                CellDefinition.copy( parent, cd );
                p = cd.phenotype;
            }

            p.sync( m );
            Element phenotypeElement = findElement( cdElement, "phenotype" );
            for( Element el : getAllElements( phenotypeElement ) )
            {
                switch( el.getTagName() )
                {
                    case "cycle":
                        readCycle( el, cd );
                        break;
                    case "death":
                        readDeath( el, p );
                        break;
                    case "volume":
                        readVolume( el, p );
                        break;
                    case "mechanics":
                        readMechanics( el, p, model );
                        break;
                    case "motility":
                        readMotility( el, cd, m );
                        break;
                    case "secretion":
                        readSecretion( el, p, m );
                        break;
                    case "cell_interactions":
                        readCellInteractions( el, cd, model );
                        break;
                    case "cell_transformations":
                        readCellTransformations( el, cd, model );
                        break;
                    case "intracellular":
                        readIntracellular( el, model, cd );
                        break;
                }
            }
            Element customDataEl = findElement( cdElement, "custom_data" );
            readCustomData( customDataEl, cd );
        }
    }

    public void readIntracellular(Element el, Model model, CellDefinition cd) throws Exception
    {
        if( intracellularReader == null )
            throw new Exception( "No intracellular reader set" );
        intracellularReader.readIntracellular( filePath, el, model, cd );
    }

    public void setIntracellularReader(IntracellularReader reader)
    {
        this.intracellularReader = reader;
    }

    private void readCycle(Element el, CellDefinition cd) throws Exception
    {
        Phenotype p = cd.phenotype;
        CycleModel model = p.cycle;
        String name = getAttr( el, "name" );
        Integer code = getIntAttr( el, "code" );
        if( code == null )
            code = PhysiCellConstants.find_cycle_model_code( name );
        if( code < 0 )
            throw new Exception( "Error. Unable to identify cycle model\n" + "node.attribute(\"name\").value()\n" + "("
                    + getIntAttr( el, "code" ) + ")" );
        model.addPhase( code, name );

        if( getIntAttr( el, "code" ) != null )
        {
            if( code == PhysiCellConstants.advanced_Ki67_cycle_model )
            {
                p.cycle = StandardModels.Ki67_advanced.clone();
            }
            else if( code == PhysiCellConstants.basic_Ki67_cycle_model )
            {
                p.cycle = StandardModels.Ki67_basic.clone();
            }
            else if( code == PhysiCellConstants.flow_cytometry_cycle_model )
            {
                p.cycle = StandardModels.flow_cytometry_cycle_model.clone();
            }
            else if( code == PhysiCellConstants.live_apoptotic_cycle_model ) // ?
            {
                p.cycle = StandardModels.live.clone(); // ?
                System.out.println( "Warning: live_apoptotic_cycle_model not directly supported.\n"
                        + "         Substituting live cells model. Set death rates=0." );
            }
            else if( code == PhysiCellConstants.total_cells_cycle_model )
            {
                p.cycle = StandardModels.live.clone();
                System.out.println( "Warning: total_cells_cycle_model not directly supported.\n"
                        + "         Substituting live cells model. Set death rates=0." );
            }
            else if( code == PhysiCellConstants.live_cells_cycle_model )
            {
                p.cycle = StandardModels.live.clone();
            }
            else if( code == PhysiCellConstants.flow_cytometry_separated_cycle_model )
            {
                p.cycle = StandardModels.flow_cytometry_separated_cycle_model.clone();
            }
            else if( code == PhysiCellConstants.cycling_quiescent_model )
            {
                p.cycle = StandardModels.cycling_quiescent.clone();
            }
            else
            {
                throw new Exception( "Warning: Unknown cycle model " );
            }
            //            cd.phenotype.cycle.sync_to_cycle_model( cd.functions.cycle_model );
        }

        Element transitionElement = findElement( el, "phase_transition_rates" );
        //        String units = getAttr( transitionElement, "units" );
        if( transitionElement != null )
        {
            for( Element rate : findAllElements( transitionElement, "rate" ) )
            {
                try
                {
                    int startIndex = getIntAttr( rate, "start_index" );
                    int endIndex = getIntAttr( rate, "end_index" );
                    boolean fixedDuration = getBoolAttr( rate, "fixed_duration" );
                    double rateValue = getDoubleVal( rate );
                    p.cycle.setBasicTransitionRate( startIndex, endIndex, rateValue );
                    p.cycle.phase_link( startIndex, endIndex ).fixedDuration = fixedDuration;
                }
                catch( Exception ex )
                {
                    ex.printStackTrace();//TODO: change to something nice
                }
            }
        }
        Element durationElement = findElement( el, "phase_durations" );
        if( durationElement != null )
        {
            for( Element duration : findAllElements( durationElement, "duration" ) )
            {
                int start = getIntAttr( duration, "index" );
                boolean fixed = getBoolAttr( duration, "fixed_duration" );
                double value = getDoubleVal( duration );
                p.cycle.data.setExitRate( start, 1.0 / ( value + 1e-16 ) );
                p.cycle.phaseLinks.get( start ).get( 0 ).fixedDuration = fixed;
            }
        }
    }

    private void readDeath(Element el, Phenotype p) throws Exception
    {
        Death death = p.death;
        for( Element modelElement : findAllElements( el, "model" ) )
        {
            Integer code = getIntAttr( modelElement, "code" );
            String name = getAttr( modelElement, "name" );

            if( code == null )
                code = PhysiCellConstants.find_cycle_model_code( name );
            if( code < 0 )
                throw new Exception( "Can not identify death model " + name );

            int death_index = death.findDeathModelIndex( code );

            boolean alreadyExists = false;
            if( death.models.get( death_index ).code == code )
            {
                alreadyExists = true;
            }

            DeathParameters parameters = new DeathParameters();
            Element deathRateElement = findElement( modelElement, "death_rate" );
            double rateValue = getDoubleVal( deathRateElement );
            parameters.time_units = getAttr( deathRateElement, "units" );
            Element durationsElement = findElement( modelElement, "phase_durations" );
            if( durationsElement != null )
            {
                //            String durationUnits = getAttr( durationsElement, "units" );
                for( Element durationElement : findAllElements( durationsElement, "duration" ) )
                {
                    int index = getIntAttr( durationElement, "index" );
                    boolean fixedDuration = getBoolAttr( durationElement, "fixed_duration" );
                    double durationValue = getDoubleVal( durationElement );
                    p.death.models.get( death_index ).data.setExitRate( index, 1.0 / ( durationValue + 1e-16 ) );
                    p.death.models.get( death_index ).phaseLinks.get( index ).get( 0 ).fixedDuration = fixedDuration;
                }
            }

            Element transitionsElement = findElement( modelElement, "phase_transition_rates" );
            if( transitionsElement != null )
            {
                //            String durationUnits = getAttr( durationsElement, "units" );
                for( Element transitionElement : findAllElements( transitionsElement, "rate" ) )
                {
                    int start = getIntAttr( transitionElement, "start_index" );
                    int end = getIntAttr( transitionElement, "end_index" );
                    boolean fixed = getBoolAttr( transitionElement, "fixed_duration" );
                    double value = getDoubleVal( transitionElement );
                    p.death.models.get( death_index ).setBasicTransitionRate( start, end, value );
                    p.death.models.get( death_index ).phase_link( start, end ).fixedDuration = fixed;
                }
            }

            Element parametersElement = findElement( modelElement, "parameters" );
            if( parametersElement != null )
            {
                for( Element paramEl : getAllElements( parametersElement ) )
                {
                    switch( paramEl.getTagName() )
                    {
                        case "unlysed_fluid_change_rate":
                            parameters.unlysed_fluid_change_rate = getDoubleVal( paramEl );
                            //                        String unlysed_fluid_change_rate_units = getAttr( paramEl, "units" ); //TODO
                            break;
                        case "lysed_fluid_change_rate":
                            parameters.lysed_fluid_change_rate = getDoubleVal( paramEl );
                            //                        String lysed_fluid_change_rate_units = getAttr( paramEl, "units" );
                            break;
                        case "cytoplasmic_biomass_change_rate":
                            parameters.cytoplasmic_biomass_change_rate = getDoubleVal( paramEl );
                            //                        String cytoplasmic_biomass_change_rate_units = getAttr( paramEl, "units" );
                            break;
                        case "nuclear_biomass_change_rate":
                            parameters.nuclear_biomass_change_rate = getDoubleVal( paramEl );
                            //                        String nuclear_biomass_change_rate_units = getAttr( paramEl, "units" );
                            break;
                        case "calcification_rate":
                            parameters.calcification_rate = getDoubleVal( paramEl );
                            //                        String calcification_rate_units = getAttr( paramEl, "units" );
                            break;
                        case "relative_rupture_volume":
                            parameters.relative_rupture_volume = getDoubleVal( paramEl );
                            //                        String relative_rupture_volume_units = getAttr( paramEl, "units" );
                            break;
                    }
                }
            }
            if( code == PhysiCellConstants.apoptosis_death_model )
            {
                if( !alreadyExists )
                {
                    death.addDeathModel( rateValue, StandardModels.apoptosis, parameters );
                    death_index = death.findDeathModelIndex( code );
                }
                else
                {
                    death.parameters.set( death_index, parameters );
                    death.rates.set( death_index, rateValue );
                }
            }
            else if( code == PhysiCellConstants.necrosis_death_model ) // set necrosis parameters 
            {
                if( !alreadyExists )
                {
                    death.addDeathModel( rateValue, StandardModels.necrosis, parameters );
                    death_index = death.findDeathModelIndex( code );
                }
                else
                {
                    death.parameters.set( death_index, parameters );
                    death.rates.set( death_index, rateValue );
                }
            }
            else if( code == PhysiCellConstants.autophagy_death_model )
            {
                System.out.println( "Warning: autophagy_death_model not yet supported.<br>         Skipping this model." );
            }
            else
            {
                throw new Exception( "Warning: Unknown death model " );
            }
        }
    }

    private void readVolume(Element el, Phenotype p)
    {
        Volume volume = p.volume;
        for( Element paramElement : getAllElements( el ) )
        {
            switch( paramElement.getTagName() )
            {
                case "total":
                    volume.total = getDoubleVal( paramElement );
                    //                    String total_units = getAttr( paramElement, "units" );
                    break;
                case "fluid_fraction":
                    volume.fluid_fraction = getDoubleVal( paramElement );
                    //                    String fluid_units = getAttr( paramElement, "units" );
                    break;
                case "nuclear":
                    volume.nuclear = getDoubleVal( paramElement );
                    //                    String nuclear_units = getAttr( paramElement, "units" );
                    break;
                case "fluid_change_rate":
                    volume.fluid_change_rate = getDoubleVal( paramElement );
                    //                    String fluid_change_rate_units = getAttr( paramElement, "units" );
                    break;
                case "cytoplasmic_biomass_change_rate":
                    volume.cytoplasmic_biomass_change_rate = getDoubleVal( paramElement );
                    //                    String cytoplasmic_biomass_change_rate_units = getAttr( paramElement, "units" );
                    break;
                case "nuclear_biomass_change_rate":
                    volume.nuclear_biomass_change_rate = getDoubleVal( paramElement );
                    //                    String nuclear_biomass_change_rate_units = getAttr( paramElement, "units" );
                    break;
                case "calcified_fraction":
                    volume.calcified_fraction = getDoubleVal( paramElement );
                    //                    String calcified_fraction_units = getAttr( paramElement, "units" );
                    break;
                case "calcification_rate":
                    volume.calcification_rate = getDoubleVal( paramElement );
                    //                    String calcification_rate_units = getAttr( paramElement, "units" );
                    break;
                case "relative_rupture_volume":
                    volume.relative_rupture_volume = getDoubleVal( paramElement );
                    //                    String relative_rupture_volume_units = getAttr( paramElement, "units" );
                    break;
            }
        }
        volume.adjust();
        p.geometry.update( null, p, 0.0 );
    }

    private void readMechanics(Element el, Phenotype p, Model m) throws Exception
    {
        Mechanics mechanics = p.mechanics;
        for( Element paramElement : getAllElements( el ) )
        {
            switch( paramElement.getTagName() )
            {
                case "cell_cell_adhesion_strength":
                    mechanics.cellCellAdhesionStrength = getDoubleVal( paramElement );
                    //                    String cell_cell_adhesion_strength_units = getAttr( paramElement, "units" );
                    break;
                case "cell_cell_repulsion_strength":
                    mechanics.cellCellRepulsionStrength = getDoubleVal( paramElement );
                    //                    String cell_cell_repulsion_strength_units = getAttr( paramElement, "units" );
                    break;
                case "relative_maximum_adhesion_distance":
                    mechanics.relMaxAdhesionDistance = getDoubleVal( paramElement );
                    //                    String relative_maximum_adhesion_distance_units = getAttr( paramElement, "units" );
                    break;
                case "cell_BM_adhesion_strength":
                    mechanics.cellBMAdhesionStrength = getDoubleVal( paramElement );
                    //                    String cell_BM_adhesion_strength_units = getAttr( paramElement, "units" );
                    break;
                case "cell_BM_repulsion_strength":
                    mechanics.cellBMRepulsionStrength = getDoubleVal( paramElement );
                    //                    String cell_BM_repulsion_strength_units = getAttr( paramElement, "units" );
                    break;
                case "attachment_elastic_constant":
                    mechanics.attachmentElasticConstant = getDoubleVal( paramElement );
                    //                    String attachment_elastic_constant_units = getAttr( paramElement, "units" );
                    break;
                case "attachment_rate":
                    mechanics.attachmentRate = getDoubleVal( paramElement );
                    //                    String attachment_rate_units = getAttr( paramElement, "units" );
                    break;
                case "detachment_rate":
                    mechanics.detachmentRate = getDoubleVal( paramElement );
                    //                    String detachment_rate_units = getAttr( paramElement, "units" );
                    break;
                case "cell_adhesion_affinities":
                    for( Element adhesionElement : findAllElements( paramElement, "cell_adhesion_affinity" ) )
                    {
                        String target = getAttr( adhesionElement, "name" );
                        double value = getDoubleVal( adhesionElement );
                        int ind = m.findCellDefinitionIndex( target );
                        if( ind > -1 )
                            mechanics.cellAdhesionAffinities[ind] = value;
                        else
                            throw new Exception( "Unknown Cell Definition " + target );
                        //                                                { std::cout << "what?!?" << std::endl; }
                    }
                    break;
                case OPTIONS:
                    for( Element optionElement : getAllElements( paramElement ) )
                    {
                        switch( optionElement.getTagName() )
                        {
                            case "set_relative_equilibrium_distance":
                                //                                String set_relative_equilibrium_distance_units = getAttr( optionElement, "units" );
                                boolean set_relative_equilibrium_distance_units_isEnabled = getBoolAttr( optionElement, "enabled" );
                                double set_relative_equilibrium_distance_units_value = getDoubleVal( optionElement );
                                if( set_relative_equilibrium_distance_units_isEnabled )
                                    mechanics.setRelEquilibriumDistance( set_relative_equilibrium_distance_units_value );
                                break;
                            case "set_absolute_equilibrium_distance":
                                //                                String set_absolute_equilibrium_distance_units = getAttr( optionElement, "units" );
                                boolean set_absolute_equilibrium_distance_units_isEnabled = getBoolAttr( optionElement, "enabled" );
                                double set_absolute_equilibrium_distance_value = getDoubleVal( optionElement );
                                if( set_absolute_equilibrium_distance_units_isEnabled )
                                    mechanics.setAbsEquilibriumDistance( p, set_absolute_equilibrium_distance_value );
                                break;
                        }

                    }
            }
        }
    }

    private void readMotility(Element el, CellDefinition cd, Microenvironment m)
    {
        Phenotype p = cd.phenotype;
        Motility motility = p.motility;

        for( Element child : getAllElements( el ) )
        {
            switch( child.getTagName() )
            {
                case "speed":
                    motility.migrationSpeed = getDoubleVal( child );
                    //                    String speed_units = getAttr( child, "units" );
                    break;
                case "persistence_time":
                    motility.persistenceTime = getDoubleVal( child );
                    //                    String persistence_time_units = getAttr( child, "units" );
                    break;
                case "migration_bias":
                    motility.migrationBias = getDoubleVal( child );
                    //                    String migration_bias_units = getAttr( child, "units" );
                    break;
                case OPTIONS:
                {
                    Element enabledElement = findElement( child, "enabled" );
                    if( enabledElement != null )
                        motility.isMotile = getBoolVal( enabledElement );

                    Element use2dElement = findElement( child, "use_2D" );
                    if( use2dElement != null )
                        motility.restrictTo2D = getBoolVal( use2dElement );

                    if( m.options.simulate2D && !motility.restrictTo2D )
                    {
                        System.out.println( "Note: Overriding to set cell motility for " + cd.name
                                + " to 2D based on microenvironment domain settings ... " );
                        motility.restrictTo2D = true;
                        break;
                    }
                    Element chemotaxisElement = findElement( child, "chemotaxis" );
                    if( chemotaxisElement != null )
                    {
                        boolean enabled = getBoolVal( findElement( chemotaxisElement, "enabled" ) );
                        if( enabled )
                        {
                            cd.functions.updateMigration = new Chemotaxis();

                            // search for the right chemo index               
                            String substrate_name = getVal( findElement( chemotaxisElement, "substrate" ) );

                            motility.chemotaxisIndex = m.findDensityIndex( substrate_name );
                            if( motility.chemotaxisIndex < 0 )
                            {
                                System.out.println( "Error: parsing phenotype:motility:options:chemotaxis:  invalid substrate" );
                                System.out.println( "Substrate " + substrate_name + " was not found in the microenvironment." );
                            }
                            String actual_name = m.densityNames[motility.chemotaxisIndex];
                            if( !substrate_name.equals( actual_name ) )
                            {
                                System.out.println( "Error: attempted to set chemotaxis to \"" + substrate_name
                                        + "\", which was not found in the microenvironment."
                                        + " Please double-check your substrate name in the config file." );
                            }
                            motility.chemotaxisDirection = getIntVal( findElement( chemotaxisElement, "direction" ) );
                        }
                    }
                    Element advancedChemotaxisElement = findElement( child, "advanced_chemotaxis" );
                    if( advancedChemotaxisElement != null )
                    {
                        if( getBoolVal( findElement( advancedChemotaxisElement, "enabled" ) ) )
                        {
                            //                            cd.functions.update_migration_bias = new StandardModels.chemotaxis_function();
                            if( cd.functions.updateMigration instanceof Chemotaxis )
                            {
                                System.out.println( "Warning: when processing motility for " + cd.name + " cells: \n"
                                        + "\tBoth chemotaxis and advanced_chemotaxis are enabled.\n"
                                        + "\tThe settings for advanced_chemotaxis override those of chemotaxis." );
                            }

                            if( getBoolAttr( advancedChemotaxisElement, "normalize_each_gradient" ) )
                                cd.functions.updateMigration = new AdvancedChemotaxisNormalized();
                            else
                                cd.functions.updateMigration = new AdvancedChemotaxis();

                            Element sensitivityEl = findElement( advancedChemotaxisElement, "chemotactic_sensitivities" );
                            if( sensitivityEl != null )
                            {
                                for( Element sensEl : findAllElements( sensitivityEl, "chemotactic_sensitivity" ) )
                                {
                                    String substrate_name = getAttr( sensEl, "substrate" );
                                    int index = m.findDensityIndex( substrate_name );
                                    String actual_name = "";
                                    if( index > -1 )
                                    {
                                        actual_name = m.densityNames[index];
                                    }

                                    // error check 
                                    if( !substrate_name.equals( actual_name ) )
                                    {
                                        System.out.println( "Warning: when processing advanced chemotaxis for " + cd.name + " cells: "
                                                + "\tInvalid substrate " + substrate_name + " specified."
                                                + "\tIgnoring this invalid substrate in the chemotaxis function .. " );
                                    }
                                    else
                                    {
                                        cd.phenotype.motility.chemotacticSensitivities[index] = getDoubleVal( sensEl );
                                    }
                                }
                            }
                            else
                            {
                                System.out.println( "Warning: when processing motility for " + cd.name + " cells: "
                                        + "\tAdvanced chemotaxis requries chemotactic_sensitivities."
                                        + "\tBut you have none. Your migration bias will be the zero vector." );
                            }
                        }
                    }
                }
            }
        }
        // display summary for diagnostic help 
        //        if( cd.functions.update_migration_bias instanceof Chemotaxis && motility.isMotile )
        //        {
        //            System.out.println( "Cells of type " + cd.name + " use standard chemotaxis: \n" + "\t d_bias (before normalization) = "
        //                    + motility.chemotaxisDirection + " * grad(" + m.density_names[motility.chemotaxisIndex] + ")" );
        //        }

        //        if( cd.functions.update_migration_bias instanceof AdvancedChemotaxis && motility.isMotile )
        //        {
        //            int number_of_substrates = m.density_names.length;

        //            System.out.println( "Cells of type " + cd.name + " use advanced chemotaxis: \n" + "\t d_bias (before normalization) = "
        //                    + motility.chemotacticSensitivities[0] + " * grad(" + m.density_names[0] + ")" );

        //            for( int n = 1; n < number_of_substrates; n++ )
        //                System.out.println( motility.chemotacticSensitivities[n] + " * grad(" + m.density_names[n] + ")" );
        //        }
        //
        //        if( cd.functions.update_migration_bias instanceof AdvancedChemotaxisNormalized && motility.isMotile )
        //        {
        //            int number_of_substrates = m.density_names.length;
        //
        //            //            System.out.println( "Cells of type " + cd.name + " use normalized advanced chemotaxis: \n"
        //            //                    + "\t d_bias (before normalization) = " + motility.chemotacticSensitivities[0] + " * grad(" + m.density_names[0] + ")"
        //            //                    + " / ||grad(" + m.density_names[0] + ")||" );
        //
        //            for( int n = 1; n < number_of_substrates; n++ )
        //            {
        //                System.out.println( motility.chemotacticSensitivities[n] + " * grad(" + m.density_names[n] + ")" + " / ||grad("
        //                        + m.density_names[n] + ")||" );
        //            }
        //        }
    }

    private void readSecretion(Element el, Phenotype p, Microenvironment m)
    {
        Secretion secretion = p.secretion;
        for( Element substrateElement : findAllElements( el, "substrate" ) )
        {
            String name = getAttr( substrateElement, "name" );
            int index = m.findDensityIndex( name );
            String actual_name = m.densityNames[index];
            if( !name.equals( actual_name ) )
            {
                System.out.println( "Error: attempted to set secretion/uptake/export for \"" + name
                        + "\", which was not found in the microenvironment.\n"
                        + "       Please double-check your substrate name in the config file." );
            }
            for( Element child : getAllElements( substrateElement ) )
            {
                switch( child.getTagName() )
                {
                    case "secretion_rate":
                        secretion.secretionRates[index] = getDoubleVal( child );
                        //                        String secreion_rate_units = getAttr( child, "units" );
                        break;
                    case "secretion_target":
                        secretion.saturationDensities[index] = getDoubleVal( child );
                        //                        String secretion_target_units = getAttr( child, "units" );
                        break;
                    case "uptake_rate":
                        secretion.uptakeRates[index] = getDoubleVal( child );
                        //                        String uptake_rate_units = getAttr( child, "units" );
                        break;
                    case "net_export_rate":
                        secretion.netExportRates[index] = getDoubleVal( child );
                        //                        String net_export_rate_units = getAttr( child, "units" );
                        break;
                }
            }
        }
    }

    private void readCellInteractions(Element el, CellDefinition cd, Model m)
    {
        CellInteractions cellInteractions = cd.phenotype.cellInteractions;
        Element dprElement = findElement( el, "dead_phagocytosis_rate" );
        if( dprElement != null )
        {
            cellInteractions.deadPhagocytosisRate = getDoubleVal( dprElement );
        }
        Element lprsElement = findElement( el, "live_phagocytosis_rates" );
        if( lprsElement != null )
        {
            for( Element lprElement : findAllElements( lprsElement, "phagocytosis_rate" ) )
            {
                // get the name of the target cell type
                String target_name = getAttr( lprElement, "name" );
                // now find its index 
                int index = m.findCellDefinitionIndex( target_name );
                // safety first! 
                if( index >= 0 )
                {
                    // if the target is found, set the appropriate rate 
                    cellInteractions.livePhagocytosisRates[index] = getDoubleVal( lprElement );
                    //                    String units = getAttr( lprElement, "units" );
                }
                else
                {
                    System.out.println( "Warning: When processing the " + cd.name + " cell definition: \n" + "\tCould not find cell type "
                            + target_name + " for phagocytosis.\n" + "\tIgnoring this live phagocytosis rate!" );
                }
            }
        }
        Element attackRatesElement = findElement( el, "attack_rates" );
        if( attackRatesElement != null )
        {
            for( Element arElement : findAllElements( attackRatesElement, "attack_rate" ) )
            {
                // get the name of the target cell type
                String target_name = getAttr( arElement, "name" );
                //            // now find its index 
                int index = m.findCellDefinitionIndex( target_name );
                if( index >= 0 )
                {
                    cellInteractions.attackRates[index] = getDoubleVal( arElement );
                }
                else
                {
                    System.out.println( "Warning: When processing the " + cd.name + " cell definition: \n" + "\tCould not find cell type "
                            + target_name + " for cell attack.\n" + "\tIgnoring this cell attack rate!" );
                }
            }
        }
        Element fusionRatesElement = findElement( el, "fusion_rates" );
        if( fusionRatesElement != null )
        {
            for( Element drElement : findAllElements( fusionRatesElement, "fusion_rate" ) )
            {
                // get the name of the target cell type
                String target_name = getAttr( drElement, "name" );
                //            // now find its index 
                int index = m.findCellDefinitionIndex( target_name );
                if( index >= 0 )
                {
                    cellInteractions.fusionRates[index] = getDoubleVal( drElement );
                }
                else
                {
                    System.out.println( "Warning: When processing the " + cd.name + " cell definition: \n" + "\tCould not find cell type "
                            + target_name + " for cell fusion.\n" + "\tIgnoring this cell fusion rate!" );
                }
            }

            Element damageRateElement = findElement( el, "damage_rate" );
            if( damageRateElement != null )
            {
                cellInteractions.damageRate = getDoubleVal( damageRateElement );
                //                String units = getAttr( damageRateElement, "units" );
            }
        }
    }

    private void readCellTransformations(Element el, CellDefinition cd, Model m)
    {
        CellTransformations transformations = cd.phenotype.cellTransformations;
        Element ratesElement = findElement( el, "transformation_rates" );
        for( Element rateElement : findAllElements( ratesElement, "transformation_rate" ) )
        {
            String name = getAttr( rateElement, "name" );
            int index = m.findCellDefinitionIndex( name );
            if( index >= 0 )
            {
                double rate = getDoubleVal( rateElement );
                //              String units = getAttr( rateElement, "units" );
                if( name == cd.name && rate > 1e-16 )
                {
                    System.out.println( "Warning: When processing the " + cd.name + " cell definition:\n" + "\tTransformation from "
                            + cd.name + " to " + name + " is not allowed.\n" + "\tIgnoring this cell transformation rate!" );
                }
                else
                    transformations.transformationRates[index] = rate;
            }
            else
            {
                System.out.println( "Warning: When processing the " + cd.name + " cell definition:\n" + "\tCould not find cell type " + name
                        + " for cell transformation.\n" + "\tIgnoring this cell transformation rate!" );

            }
        }
    }

    private void readCustomData(Element el, CellDefinition cd) throws Exception
    {
        for( Element child : getAllElements( el ) )
        {
            CellDefinition defaults = StandardModels.getDefaultCellDefinition();
            // name of teh custom data 
            String name = child.getTagName();

            // units 
            String units = getAttr( child, "units" );
            double[] values; // ( length, 0.0 ); 

            // conserved quantity 
//            boolean conserved = getBoolAttr( child, "conserved" );

            // get value(s)
            String str_values = getVal( child );
            values = VectorUtil.csv_to_vector( str_values );

            // add variable if cell defaults  
            // if the custom data is not yet found, add it 
            // first, try scalar 
            if( values.length == 1 )
            {
                // find the variable 
                int n = cd.custom_data.findVariableIndex( name );
                // if it exists, overwrite 
                if( n > -1 )
                {
                    cd.custom_data.variables.get( n ).value = values[0];
                }
                // otherwise, add 
                else
                {
                    cd.custom_data.addVariable( name, units, values[0] );
                    defaults.custom_data.addVariable( name, units, values[0] );
                }

                n = cd.custom_data.findVariableIndex( name );
            }
            // or vector 
            else
            {
                // find the variable 
                int n = cd.custom_data.findVectorVariableIndex( name );
                // if it exists, overwrite 
                if( n > -1 )
                {
                    cd.custom_data.vectorVariables.get( n ).value = values;
                }
                // otherwise, add 
                else
                {
                    cd.custom_data.addVectorVariable( name, units, values );
                }

                n = cd.custom_data.findVectorVariableIndex( name );
            }
        }
    }

    public void readUserParameters(Element el, Model model)
    {
        Element parametersElement = findElement( el, "user_parameters" );
        if( parametersElement != null )
        {
            for( Element element : getAllElements( parametersElement ) )
            {
                if( element.hasChildNodes() )
                    model.addParameter( element.getTagName(), getVal( element ), getAttr( element, "description" ) );
            }
        }
    }

    public void readInitialConditions(Element el, Model m) throws Exception
    {
        Element parametersElement = findElement( el, "initial_conditions" );
        if( parametersElement != null )
        {
            Element positionsElement = findElement( parametersElement, "cell_positions" );
            String type = getAttr( positionsElement, "format" );
            boolean enabled = getBoolAttr( positionsElement, "enabled" );
            if( !enabled )
                return;
            String folder = getVal( findElement( positionsElement, "folder" ) );
            String fileName = getVal( findElement( positionsElement, "filename" ) );
            String inputFilename = folder + "/" + fileName;
            m.setInitialInfo( new ExternalFile( type, inputFilename ) );

            if( readFromJAR )
            {
                CellCSVReader.load_cells_csv( m.getClass().getResourceAsStream( inputFilename ), m );
            }
        }
    }

    public static class ExternalFile
    {
        public String format;
        public String path;

        public ExternalFile(String format, String path)
        {
            this.format = format;
            this.path = path;
        }
    }

    public void setFunctionsReader(FunctionsReader reader)
    {
        this.functionsReader = reader;
    }

    public void readFunctions(Element el, CellDefinition cd)
    {
        if( functionsReader != null )
            functionsReader.readFunctions( el, cd );
    }
}